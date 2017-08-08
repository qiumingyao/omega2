//============================================================================
// Name        : Config.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Config cpp file
//============================================================================

#include "Config.h"


//=============================================================================
// Default constructor
//=============================================================================
Config::Config() {
}


//=============================================================================
// Default destructor
//=============================================================================
Config::~Config() {
}

//=============================================================================
//Here is to create and initialize the static members
//=============================================================================
vector<string> Config::readFilenamesList;
vector<string> Config::edgeFilenamesList;
string Config::outFilenamePrefix = "new_omega_out";
string Config::pacbioSamFile = "";

// global variables with default value
unsigned int minReadsCountInEdgeToBeNotDeadEnd = 10;
unsigned int minEdgeLengthToBeNotDeadEnd = 500;
unsigned int minReadsCountInEdgeToBe1MinFlow = 20;
unsigned int minEdgeLengthToBe1MinFlow = 1000;
unsigned int minContigLengthTobeReported = 500;
unsigned int minOvlDiffToClip = 50;
unsigned int minFoldToBeShortBranch = 5;
bool reportUnoverlappedReadsToContigs = false;
UINT64 Config::minOvl(0);
UINT64 Config::maxSubs(0);
UINT64 Config::maxEdits(0);

//=============================================================================
// print help usage
//=============================================================================
void Config::printHelp()
{
	cout << endl
		<< "  Usage:" << endl
		<< "    new_omega [OPTION]...<PARAM>..." << endl
		<< endl
		<< "  <PARAM>" << std::endl
		<< "    -f\t contained read reduction read filename(s) (comma separated fasta/fastq)" << endl
		<< "    -e\t ovelapped edge property graph filename(s) (comma separated edge list)" << endl
		<< endl
		<< "  [OPTION]" << std::endl
		<< "    -h/--help\t only print out the help contents" << endl
		<< "    -ovl\t minimum overlap length (default: 0, use all overlap found in edge property graph files)" << endl
		<< "    -o\t\t all output filename prefix (default: new_omega_out)" << endl
		<< "    -log\t verbosity level of log messages: ERROR, WARNING, INFO (default: INFO)" << endl
		<< "    -sam\t sam file with Illumina reads mapped to pacbio reads (default: none)" << endl
		<< "    -sub\t maximum allowed substitutions in overlap (default: 0)" << endl
		<< "    -edit\t maximum allowed edits in overlap (default: 0)" << endl
		<< "    -mcd\t Minimum reads Count in an edge regarding to be not a Dead-end (default: 10)" << endl
		<< "    -mld\t Minimum edge Length regarding to be not a Dead-end (default: 500)" << endl
		<< "    -mcf\t Minimum reads Count in an edge regarding to be 1 minimum Flow (default: 20)" << endl
		<< "    -mlf\t Minimum edge Length regarding to be 1 minimum Flow (default: 1000)" << endl
		<< "    -mlr\t Minimum reads Length to be Reported (default: 200)" << endl
		<< "    -mdc\t Minimum overlap length difference to clip branches (default: 50)" << endl
		<< "    -mfs\t Minimum fold difference to consider branches to be short (default: 5)" << endl
//		<< "    -rur\t Report Unoverlapped Reads to contigs (default: false)" << endl
		<< endl;
}

bool Config::setConfig(int argc, char **argv)
{
	vector<string> argumentsList;
	cout << "PRINTING ARGUMENTS" << endl;
	for(int i = 0; i < argc; i++)
	{
		cout << argv[i] << ' ';
	}
	cout << endl;
	while(argc--)
		argumentsList.push_back(*argv++);

	if(argumentsList.size() == 1)
	{
		Config::printHelp();
		return false;
	}

	for(UINT64 i = 1; i <= argumentsList.size()-1; i++)
	{
		// -h/--help: only print out the help contents
		if (argumentsList[i] == "-h" || argumentsList[i] == "--help")
		{
			Config::printHelp();
			exit(0);
		}
		// -f: contained read reduction read filename(s) (comma separated fasta/fastq)
		else if (argumentsList[i] == "-f") {
			string inputFilenames=argumentsList[++i];
			stringstream ss(inputFilenames);
			string item;

			while (getline(ss, item, ','))
			{
				Config::readFilenamesList.push_back(item);
			}
		}
		// -e: ovelapped edge property graph filename(s) (comma separated edge list)
		else if (argumentsList[i] == "-e") {
			string inputFilenames=argumentsList[++i];
			stringstream ss(inputFilenames);
			string item;

			while (getline(ss, item, ','))
			{
				Config::edgeFilenamesList.push_back(item);
			}
		}
		// -ovl: minimum overlap length (default: 0)
		else if (argumentsList[i] == "-ovl") {
			Config::minOvl = stoi(argumentsList[++i]);
		}
		else if (argumentsList[i] == "-sub") {
			Config::maxSubs = stoi(argumentsList[++i]);
		}
		else if (argumentsList[i] == "-edit") {
			Config::maxEdits = stoi(argumentsList[++i]);
		}
		// -log: verbosity level of log messages
		else if (argumentsList[i] == "-log") {
			FILELog::ReportingLevel() = FILELog::FromString(argumentsList[++i]);
		}
		else if (argumentsList[i] == "-sam")
		{
			Config::pacbioSamFile = argumentsList[++i];
		}
		// -o: all output filename prefix (default: new_omega_out)
		else if (argumentsList[i] == "-o") {
			Config::outFilenamePrefix = argumentsList[++i];
		}
		// -mcd: Minimum reads Count in an edge regarding to be not a Dead-end (default: 10)
		else if (argumentsList[i] == "-mcd") {
			istringstream iss_a(argumentsList[++i]);
			iss_a >> minReadsCountInEdgeToBeNotDeadEnd;
		}
		// -mld: Minimum edge Length regarding to be not a Dead-end (default: 500)
		else if (argumentsList[i] == "-mld") {
			istringstream iss_b(argumentsList[++i]);
			iss_b >> minEdgeLengthToBeNotDeadEnd;
		}
		// -mcf: Minimum reads Count in an edge regarding to be 1 minimum Flow (default: 20)
		else if (argumentsList[i] == "-mcf") {
			istringstream iss_c(argumentsList[++i]);
			iss_c >> minReadsCountInEdgeToBe1MinFlow;
		}
		// -mlf: Minimum edge Length regarding to be 1 minimum Flow (default: 1000)
		else if (argumentsList[i] == "-mlf") {
			istringstream iss_d(argumentsList[++i]);
			iss_d >> minEdgeLengthToBe1MinFlow;
		}
		// -mlr: Minimum reads Length to be Reported (default: 200)
		else if (argumentsList[i] == "-mlr") {
			istringstream iss_d(argumentsList[++i]);
			iss_d >> minContigLengthTobeReported;
		}
		// -mdc: Minimum overlap length difference to clip branches (default: 50)
		else if (argumentsList[i] == "-mdc") {
			istringstream iss_d(argumentsList[++i]);
			iss_d >> minOvlDiffToClip;
		}
		// -mfs: Minimum fold difference to consider branches to be short (default: 5)
		else if (argumentsList[i] == "-mfs") {
			istringstream iss_d(argumentsList[++i]);
			iss_d >> minFoldToBeShortBranch;
		}
		// -rre\t Report Removed Edges to contigs (default: false)
//		else if (argumentsList[i] == "-rre")
//		{
//			reportRemovedEdgesToContigs = true;
//		}
		// -rur\t Report Unoverlapped Reads to contigs (default: false)
		else if (argumentsList[i] == "-rur") {
			reportUnoverlappedReadsToContigs = true;
		}
		else
		{
			Config::printHelp();
			return false;
		}
	}
	return true;
}


////=============================================================================
//// Get read filenames
////=============================================================================
//vector<string> Config::getReadFilenames()
//{
//	return Config::readFilenamesList;
//}
//
//
////=============================================================================
//// Get edge filenames
////=============================================================================
//vector<string> Config::getEdgeFilenames()
//{
//	return Config::edgeFilenamesList;
//}
//
//
////=============================================================================
//// Get output prefix name
////=============================================================================
//string Config::getOutputFilenamePrefix()
//{
//	return Config::outFilenamePrefix;
//}
//
//UINT64 Config::getMinOvl()
//{
//	return Config::minOvl;
//}
//
