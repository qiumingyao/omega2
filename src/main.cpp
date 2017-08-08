//============================================================================
// Name        : main.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Main code
//============================================================================

#include "Config.h"
#include "DataSet.h"
#include "OverlapGraph.h"
#include "logcpp/log.h"

int OverlapGraph::s_nReads_in_goodEdges = 0;
int OverlapGraph::s_nGoodEdges = 0;
TLogLevel loglevel = logINFO;                   /* verbosity level of logging */
string outputFilenamePrefix = "omega2";

int main(int argc, char **argv) {

	// Parse command line options:
	if(!Config::setConfig(argc, argv)){
		cerr << "Error: wrong configurations" << endl;
		return false;
	}

	CLOCKSTART;
	vector<string> readFilenameList 	= Config::getReadFilenames();
	vector<string> edgeFilenameList 	= Config::getEdgeFilenames();
	outputFilenamePrefix 			= Config::getOutputFilenamePrefix();
	UINT64 minOvl 				= Config::getMinOvl();
	UINT64 maxSubs				= Config::getMaxSubs();
	UINT64 maxEdits				= Config::getMaxEdits();
	string samFileName 			= Config::pacbioSamFile;
	loglevel 				= FILELog::ReportingLevel();

	FILE_LOG(logINFO) << "Log level is " << FILELog::ReportingLevel() << ":\t" << Config::getLogLevel() << endl;
	FILE_LOG(logINFO) << "Minimum overlap length is: " << minOvl << endl;
	FILE_LOG(logINFO) << "File(s) including reads: ";
	if(loglevel > 1){
		for(vector<std::string>::iterator it = readFilenameList.begin(); it!=readFilenameList.end(); ++it)
			FILE_LOG(logINFO) << *it << "\t";
	}
	FILE_LOG(logINFO) << endl;
	FILE_LOG(logINFO) << "File(s) including edges: ";
	if(loglevel > 1){
		for(vector<std::string>::iterator it = edgeFilenameList.begin(); it!=edgeFilenameList.end(); ++it)
			FILE_LOG(logINFO) << *it;
	}
	FILE_LOG(logINFO) << endl;
	FILE_LOG(logINFO) << "Output file names' prefix is: " << outputFilenamePrefix << endl;
	FILE_LOG(logINFO) << "Maximum number of substitutions allowed in the overlap is: " << maxSubs << endl;
	FILE_LOG(logINFO) << "Maximum number of edits allowed in the overlap is: " << maxEdits << endl;
	FILE_LOG(logINFO) << "Maximum read count in dead-end edge is: " << minReadsCountInEdgeToBeNotDeadEnd << endl;
	FILE_LOG(logINFO) << "Maximum edge length in dead-end edge is: " << minEdgeLengthToBeNotDeadEnd << endl;
	FILE_LOG(logINFO) << "Minimum read count in edges with flow is: " << minReadsCountInEdgeToBe1MinFlow << endl;
	FILE_LOG(logINFO) << "Minimum edge length of edges with flow is: " << minEdgeLengthToBe1MinFlow << endl;
	FILE_LOG(logINFO) << "Minimum edge length for edges to be reported is: " << minContigLengthTobeReported << endl;
	FILE_LOG(logINFO) << "Minimum overlap length difference for branches to clip: " << minOvlDiffToClip << endl;
	FILE_LOG(logINFO) << "Minimum fold difference to consider branches to be short: " << minFoldToBeShortBranch << endl;
	if (samFileName.length() > 0){
		FILE_LOG(logINFO) << "Sam file used for long read threading is: " << samFileName << endl;
	}
	
//	FILE_LOG(logINFO) << "Report isolated reads to contigs? " << std::boolalpha << reportUnoverlappedReadsToContigs << endl;


	OverlapGraph *overlapGraph = new OverlapGraph(edgeFilenameList, readFilenameList, 
			minOvl, maxSubs, maxEdits);
	overlapGraph->simplifyGraph();

	if(loglevel > 3){
		std::string graph_file = outputFilenamePrefix+"_before_flow_graph.gdl";
		ofstream g_out(graph_file.c_str());
		g_out << *overlapGraph;
		g_out.close();
	}
	// Flow analysis
	overlapGraph->calculateFlowStream();
	overlapGraph->removeAllEdgesWithoutFlow();
	overlapGraph->simplifyGraph();

	std::string graph_file = outputFilenamePrefix+"_graph.gdl";
	ofstream g_out(graph_file.c_str());
	g_out << *overlapGraph;
	g_out.close();

	std::string contig_file = outputFilenamePrefix+"_contigs.fasta";
	ofstream f_out(contig_file.c_str());
	overlapGraph->printContigs(f_out, readFilenameList);
	f_out.close();

	// Pacbio read threading
	if(samFileName.length() > 0){
		ifstream inputSamStream(samFileName.c_str());
		UINT64 totalEdgesThreaded = overlapGraph->readsThreading(inputSamStream);
		inputSamStream.close();
		cout << "Total number of edge merging is " << totalEdgesThreaded << endl;
		overlapGraph->simplifyGraph();


		std::string graph_file = outputFilenamePrefix+"_threading_graph.gdl";
		ofstream g_out(graph_file.c_str());
		g_out << *overlapGraph;
		g_out.close();

		std::string contig_file = outputFilenamePrefix+"_threading_contigs.fasta";
		ofstream f_out(contig_file.c_str());
		overlapGraph->printContigs(f_out, readFilenameList);
		f_out.close();
	}

	delete overlapGraph;

	CLOCKSTOP;
	return 0;
}
