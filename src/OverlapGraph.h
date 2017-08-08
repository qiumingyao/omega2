#ifndef OVERLAPGRAPH_H
#define OVERLAPGRAPH_H

/*
 * ===== CLASS HEADER ========================================================
 * Name        : OverlapGraph.cpp
 * Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
 * Version     : v1.2
 * Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
 * Description : OverlapGraph header file
 *============================================================================
 */

#include "Config.h"
#include "DataSet.h"
#include "Edge.h"

typedef vector<Edge *> t_edge_vec;	// vector of pointers to Edge
class OverlapGraph
{
	private:
		/* ====================  DATA MEMBERS  ======================================= */
		DataSet 		*m_dataset;
		vector<t_edge_vec*> 	*m_graph;
		t_edge_vec 		*m_deletedEdges;	// For saving deleted edges
		UINT64 			m_numberOfNodes;
		UINT64 			m_numberOfEdges;
		UINT64 			m_minOvl;
		UINT64			m_maxSubs;
		UINT64			m_maxEdits;
		bool 			m_flowComputed;

		/* ====================  METHODS       ======================================= */

		// Get mismatch information (type t_vpair) from mismatch string
		t_vpair* getMismatchFromString(const std::string &mismath_string);

		// Insert an edge in the overlap graph
		// Does not insert its twin edge
		void insertEdge( Edge * edge);   
		void insertFwdEdge( Edge * edge);   

		// Remove an edge from the edge list of the edge's source read, 
		// but do not delete this edge from memory yet.
		void removeEdgeFromSourceRead(Edge *edge);

		// Remove an edge from the overlap graph.
		void removeEdge(Edge *edge, bool save_for_contig=false);

		void removeFwdEdge(Edge *edge);
		// Contract composite paths in the overlap graph.
		UINT64 contractCompositeEdges(void);

		// Remove dead-ends from the overlap graph.
		UINT64 removeDeadEndNodes(void);

		// Remove dead-end like short branches, if they are much shorter 
		// than other longer edges even though they are longer than
		// the dead-end length threshold
		UINT64 removeShortBranches(void);

		// remove multi-edges with similar strings
		UINT64 removeSimilarEdges(void);

		// Clip branches with much shorter overlap length
		// comparing to other branches.
		UINT64 clipBranches(void);

		// Loops that can be traversed only one way
		UINT64 reduceLoops(void);

		// Calculate bounds and costs of flow for minimum cost flow in the overlap graph.
		void calculateBoundAndCost(const Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST);

		// Find an edge from source to destination in the overlap graph.
		Edge *findEdge(const UINT64 & source, const UINT64 & destination);

		// Find edges between two reads
		vector<Edge *> findEdges(const UINT64 & source, const UINT64 & destination);

		// Given lists of edge IDs (positive and negative), merge the corresponding
		// edges in each list, and insert them in the graph.
		// Also delete these edges at the end
		UINT64 merge_edges_with_IDlists(const vector<vector<INT64>> *lists_edgeIDs_to_merge, 
				const map<INT64, Edge*> & id_to_edge);

		// Thread Illumina reads using information from one PacBio read
		void readThreading(const vector<UINT64> & readNumbers, 
				vector< t_edge_vec > & list_edges_to_merge);

		// locate a read given the p_read (composite numebr - last bit for orientation)
		Edge * locateRead(UINT64 p_read) const;

		// Merge a list of edges (instead of only two edges)
		UINT64 mergeListOfEdges(const vector< t_edge_vec > & list_edges_to_merge);

		// Load edge file
		void loadEdgesFromEdgeFile(const std::string &readFilename);

		// Sort edges of each read based on ID of the destination read. 
		// This is only for ordering edges for convenience in the output file
		void sortEdgesByLength();

		void sortEdgesByDestID();

	public:
		/* ====================  DATA MEMBERS  ======================================= */
		// Keep track of number of reads contained in the edges with 1 unit of flow assigned 
		static int s_nReads_in_goodEdges; 

		// keep track of number of edges with 1 unit of flow assigned 
		static int s_nGoodEdges;        

		/* ====================  LIFECYCLE     ======================================= */
		OverlapGraph(void);

		OverlapGraph(const vector<std::string> &edge_files, 
				const vector<std::string> &read_files, 
				const UINT64 &minOvl, const UINT64 &maxSubs, const UINT64 &maxEdits);

		~OverlapGraph();

		/* ====================  ACCESSORS     ======================================= */
		UINT64 getNumberOfEdges(void) const {return m_numberOfEdges;}

		UINT64 getNumberOfNodes(void) const {return m_numberOfNodes;}

		/* ====================  MUTATORS      ======================================= */
		void setMinOvl(const UINT64 & minOvl = 50){m_minOvl = minOvl;}

		// Remove all simple edges without flow
		UINT64 removeAllEdgesWithoutFlow();

		void loadStringFromReadsFile(const std::string &read_file, UINT64 &readID);

		void populate_read(const UINT64 &readID, const std::string & read_str);

		// Pacbio read threading from sam file
		UINT64 readsThreading(ifstream & inputSamStream);

		/* ====================  OPERATORS     ======================================= */
		friend ostream& operator<< (ostream &out, const OverlapGraph & graph);


//		// Load read file again to print contigs
//		UINT64 streamReadFileSequences(const std::string &readFilename, UINT64 &readID, 
//				vector<vector<char> > & contigStrings, 
//				vector<unordered_map<UINT32, vector<char> > > & mismatchMap);


		// Some simple simplification.
		void simplifyGraph(void);

		// Calculate the minimum cost flow of the overlap graph using file
		void calculateFlowStream(void);

		// Find all the edges in the graph, to be used in print graph and contigs
		void getEdges(t_edge_vec & contigEdges) const;
		
		// Print contigs to file, only the ones longer than the specified printing threshold
		void printContigs(ostream & out, const vector<std::string> &readFilenameList);
};


bool compareEdgesByDestID (const Edge *edge1, const Edge* edge2);
bool compareEdgesByLength (const Edge *edge1, const Edge* edge2);
#endif /* -----  end of class OverlapGraph  ----- */
