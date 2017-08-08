/*
 * ===== CLASS IMPLEMENTATION ================================================
 * Name        : OverlapGraph.cpp
 * Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
 * Version     : v1.2
 * Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
 * Description : OverlapGraph cpp file
 *============================================================================
 */

#include <set>
#include "Edge.h"
#include "OverlapGraph.h"
#include "CS2_stream/cs2.h"
#include "tRead.h"
#include "namespace_jj.h"
#include "tDataSet.h"
#include "tEdge.h"
#include "tOverlapGraph.h"
#include "Utils.h"
extern TLogLevel loglevel;                      /* verbosity level of logging */
extern std::string outputFilenamePrefix;


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getMismatchFromString
 *  Description:  From the mismatch string in the edge file, get the mismatch into t_vpair
 * =====================================================================================
 */
t_vpair* OverlapGraph::getMismatchFromString(const std::string &mismatch_string)
{
	if (mismatch_string == "NA")
		return nullptr;
	else{
	// mismatches to list
		// loop mismatches
		t_vpair *mismatches = new t_vpair;
		std::stringstream mismatches_ss(mismatch_string);
		// JJ: The source and destination characters are not used?
		// Can I remove them, and only get the position, like below?
		std::string element;
		getline(mismatches_ss, element, ':');
		std::istringstream element_ss(element);
		UINT32 mismatchPosition;
		element_ss >> mismatchPosition;
		mismatches->push_back(make_pair(0, mismatchPosition));
		return mismatches;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  insertEdge
 *  Description:  Insert an edge in the overlap graph.
 *                Does not automatically insert its twin edge.
 * =====================================================================================
 */
void OverlapGraph::insertEdge( Edge *edge)
{
	insertFwdEdge(edge);
	insertFwdEdge(edge->getReverseEdge());
}

void OverlapGraph::insertFwdEdge( Edge *edge)
{
	UINT64 ID = edge->getSourceRead()->getReadID();
	if(m_graph->at(ID)->empty())            /* A new node in the graph */
		m_numberOfNodes++;
	m_graph->at(ID)->push_back(edge);
	++m_numberOfEdges;
	edge->updateReadsLocations(INSERTION);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeEdgeFromSourceRead
 *  Description:  Remove an edge from the edge list of the edge's source read, 
 *  		  but do not delete this edge from memory yet.
 * =====================================================================================
 */
void OverlapGraph::removeEdgeFromSourceRead(Edge *edge)
{
	if (!edge)
		return;
	// all edges incident to source read
	t_edge_vec *fwd_edges = m_graph->at(edge->getSourceRead()->getReadID());	

	fwd_edges->erase(std::remove(fwd_edges->begin(), fwd_edges->end(), edge), fwd_edges->end());
	if (fwd_edges->empty())
		--m_numberOfNodes;
	--m_numberOfEdges;
}

void OverlapGraph::removeEdge(Edge *edge, bool save_for_contig){
	if(edge == nullptr)
		return;
	removeFwdEdge(edge->getReverseEdge());
	removeFwdEdge(edge);
	//TODO: need to redo
	if (save_for_contig){
		// Between the edge and its reverse edge, only save the one with smaller source read number
		if(edge->isSmallerEdge()) {
			m_deletedEdges->push_back(edge);
		}
		else{
			m_deletedEdges->push_back(edge->getReverseEdge());
		}
		FILE_LOG(logDEBUG) << "Insert deleted edge " 
			<< *(m_deletedEdges->back()) << " to vector m_deletedEdges" << "\n";
	}
	else{
		delete edge->getReverseEdge();
		delete edge;
	}
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeEdge
 *  Description:  Remove an edge from the overlap graph, release its memory if we don't 
 *  		  need to save this edge for the contigs.
 *  		  This removes and deletes from memory both the edge and its twin edge
 * =====================================================================================
 */
void OverlapGraph::removeFwdEdge(Edge *edge)
{
	// If edge points to NULL, there is nothing to to
	if(edge == nullptr) {
		return;
	}
	// If the current edge contains some reads. We have to update their location formation.
	edge->updateReadsLocations(DELETION);

	removeEdgeFromSourceRead(edge);

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  contractCompositeEdges
 *  Description:  Contract composite paths in the overlap graph.
 *                u*-------->v>---------*w  => u*--------------------*w
 *                u*--------<v<---------*w  => u*--------------------*w
 * =====================================================================================
 */
UINT64 OverlapGraph::contractCompositeEdges(void)
{
	CLOCKSTART;
	UINT64 counter(0);
	for(UINT64 index = 1 ; index < m_graph->size(); ++index)
	{
		if(m_graph->at(index)->size() == 2) // Check if the node has only two edges.
		{
			// First edge, going into Read 
			Edge *edge1 = m_graph->at(index)->at(0)->getReverseEdge();
			// Second edge, going out from Read
			Edge *edge2 = m_graph->at(index)->at(1);
			// One incoming edge and one outgoing edge.
			// And do not merge if either of the edges is a loop
			if( is_mergeable(edge1, edge2) && !(edge1->isLoop()) && !(edge2->isLoop()) ) 
			{
				Edge *new_edge = Add(edge1, edge2);
				insertEdge(new_edge);
				removeEdge(edge1, false);
				if(edge2 != edge1->getReverseEdge()){
					removeEdge(edge2, false);
				}
				++counter;	// Counter how many edges merged.
			}
		}

	}
	if(counter > 0){
		FILE_LOG(logINFO) << setw(10) << counter << " composite Edges merged." << "\n";
	}
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return counter;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeDeadEndNodes
 *  Description:  remove dead end nodes and their incident edges
 * =====================================================================================
 */
UINT64 OverlapGraph::removeDeadEndNodes(void)
{
	CLOCKSTART;
	vector<UINT64> nodes_to_remove;
	for(UINT64 i = 1; i < m_graph->size(); i++) // For each read.
	{
		if(!m_graph->at(i)->empty())	// If the read has some edges.
		{
			bool isDeadEnd = true;	// flag for dead end edge
			UINT64 inEdge = 0; 	// number of incoming edges to this node
			UINT64 outEdge = 0; 	// number of outgoing edges from this node

			// Find number of in- ane out- edges
			for(UINT64 j=0; j < m_graph->at(i)->size(); j++)
			{
				Edge * edge = m_graph->at(i)->at(j);
				/* Break case:
				 * 0. edge already marked as not dead end
				 * 1. composite edge with more than minReadsCountInEdgeToBeNotDeadEnd (deafult: 10)
				 * 2. composite edge longer than minEdgeLengthToBeNotDeadEnd (default: 500)
				 * 3. the edge is loop for the current node
				 * Then flag=1 and exit the loop
				 */
				if (edge->isNotDeadEnd()){
					isDeadEnd = false;
					break;
				}
				if(edge->getListOfReads() && edge->getListOfReads()->size() >= minReadsCountInEdgeToBeNotDeadEnd) {
					edge->markNotDeadEnd();
					isDeadEnd = false;
					break;
				}
				if(edge->getEdgeLength() >= minEdgeLengthToBeNotDeadEnd) {
					edge->markNotDeadEnd();
					isDeadEnd = false;
					break;
				}
				if(edge->isLoop())
				{
					edge->markNotDeadEnd();
					isDeadEnd = false;
					break;
				}

				if((edge->getOrientation() >> 1) & 1)
					++outEdge;
				else
					++inEdge;
			}
			// no good edges incident to the node and only in-edges or out-edges
			if( isDeadEnd && inEdge*outEdge == 0 && inEdge + outEdge > 0){
				nodes_to_remove.push_back(i);
			}
		}
	}
	FILE_LOG(logDEBUG) << "number of dead end nodes found: " << nodes_to_remove.size() << "\n";

	UINT64 deleted_edges(0);
	// Now delete the edges incident to these nodes
	for(auto it = nodes_to_remove.cbegin(); it != nodes_to_remove.cend(); ++it){
		UINT64 nodeID = *it;
		while(!(m_graph->at(nodeID)->empty()))
		{
			removeEdge(m_graph->at(nodeID)->front(), false);
			++deleted_edges;
		}
	}
	FILE_LOG(logDEBUG) << "number of edges deleted: " << deleted_edges << "\n";

	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return deleted_edges;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeShortBranches
 *  Description:  After flow analysis, at late stage of graph simplification.
 *  Now only look at nodes with 1 incident edge, if this edge is much shorter comparing
 *  to other edges incident to its neighbor, then remove this relatively short branch.
 * =====================================================================================
 */
UINT64 OverlapGraph::removeShortBranches(void)
{
	if(!(this->m_flowComputed)) {
		return 0;
	}
	CLOCKSTART;
	UINT64 num_nodes_rm(0);

	typedef vector<UINT64> LongBrLens;
	map<UINT64, LongBrLens> long_brlens_map; /* map from read number to length vector of length 2 */
	for(UINT64 i = 1; i < m_graph->size(); i++) // For each read.
	{
		/* If this read has exactly 1 incident edge, and this edge is not a loop */
		if(m_graph->at(i)->size()==1 && !(m_graph->at(i)->front()->isLoop())){

			/* The incident edge, going out from its only neighbor */
			Edge* one_edge = m_graph->at(i)->at(0)->getReverseEdge(); 

			UINT64 neighbor = one_edge->getSourceRead()->getReadID();
			UINT64 neighbor_degree = m_graph->at(neighbor)->size();

			/* If neighbor has more edges besides the "short branch" */
			if(neighbor_degree > 1){ 
				UINT64 one_length = one_edge->getOverlapOffset();
				// edge is going in (0) or out (1) from the read
				UINT8 in_out = ((one_edge->getOrientation() >> 1) & 1 );	
/* 				FILE_LOG(logDEBUG1) << "Short branch from " << neighbor << " to " << i 
 * 					<< " with orientation " << in_out << " and length " << one_length << "\n";
 */
				// If the longest branch lengths at this neighor hasn't been found yet,
				// look for it
				if(long_brlens_map.count(neighbor)==0){
/* 					FILE_LOG(logDEBUG1) << "Now look for the longest branches at node " 
 * 						<< neighbor << "\n";
 */
					LongBrLens long_brlens;
					// Initialize both longest length in and out to be 0
					long_brlens.push_back(0);
					long_brlens.push_back(0);
					long_brlens.at(in_out) = one_length;
					/* Find the longest branches at neighbor in both directions */
					for(UINT64 j = 0; j < neighbor_degree; ++j){
						Edge* e = m_graph->at(neighbor)->at(j);
						UINT8 direction = (e->getOrientation() >> 1) & 1; 
						if(e->getOverlapOffset() > long_brlens.at(direction))
							long_brlens.at(direction) = e->getOverlapOffset();

					}
					long_brlens_map[neighbor] = long_brlens;
				}
				// overlap is small 
				// edge is not already in the deleting list
				// edge is smaller than its reverse edge
				if(one_length * minFoldToBeShortBranch < long_brlens_map.at(neighbor).at(in_out)){ 
					removeEdge(one_edge,false);
					++num_nodes_rm;
					FILE_LOG(logDEBUG1) << "Delete this edge, length: " << one_length << " and " << long_brlens_map.at(neighbor).at(in_out) << "\n";
				}
			}
		}
	}

	if(num_nodes_rm > 0){
		FILE_LOG(logINFO) << "short-branch nodes removed: " << num_nodes_rm << "\n";
	}
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return num_nodes_rm;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeSimilarEdges
 *  Description:  Remove similar length edges between same pair of nodes
 * =====================================================================================
 */
UINT64 OverlapGraph::removeSimilarEdges(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	// pairs of edges that are found, the two lists have the same length, a pair are in the same index

	for(UINT64 i = 1; i < m_graph->size(); i++){
		if(m_graph->at(i)->size() > 1){
			UINT64 num_edges(m_graph->at(i)->size());
			for(UINT64 j = 0; j < num_edges; j++)	// For all the edges.
			{
				Edge * e1 = m_graph->at(i)->at(j);
				UINT64 destination1(e1->getDestinationRead()->getReadID());

				// Only check if edge's source is smaller than destination, 
				// and not already in the to remove list
				if( !(e1->isLoop()) ) {
					for(UINT64 k = j + 1; k < num_edges; k++) {

						Edge * e2 = m_graph->at(i)->at(k);
						// Source read is the same, check if 
						// 1. destination read is the same
						// 2. orientation is the same
						if(destination1 == e2->getDestinationRead()->getReadID() ){
							if(e1->getOrientation() == e2->getOrientation()) {
								// The lengths are more than 95% similar
								if(abs((int)(e1->getOverlapOffset() - e2->getOverlapOffset())) < (int)(e2->getOverlapOffset()/20)) 
								{
									//FILE_LOG(logDEBUG1) << *e1 << " and\n " << *e2 << " are similar\n\n";
									e1->getBaseByBaseCoverage();
									e2->getBaseByBaseCoverage();
									UINT64 e1_reads = (e1->getListOfReads() ? e1->getListOfReads()->size() : 0);
									UINT64 e2_reads = (e2->getListOfReads() ? e2->getListOfReads()->size() : 0);
									// Check coverage depth and number of reads 
									// to decide which one to keep
									if(e1->getCovDepth() < e2->getCovDepth() ||	
											(e1->getCovDepth() == e2->getCovDepth() && 
											 e1_reads < e2_reads))
									{
										// remove e1
										removeEdge(e1,false);
										e1 = e2;
										// if e1 is removed from the list of edges, 
										// the index of next edge decreases by 1
										--j;
									}
									else{
										// remove e2
										removeEdge(e2,false);
									}
									--k;
									++counter;
									--num_edges;
								}
							}
						}
						// destination is not the same any more
						else{
							break;
						}
					}
				}
			}
		}
	}
	if(counter > 0){
		FILE_LOG(logINFO) << counter << " edges removed during bubble popping." << "\n";
	}
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return counter;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  clipBranches
 *  Description:  When a node has multiple in- or out- edges, clip the edges with small
 *  		  overlap lengths
 * =====================================================================================
 */
UINT64 OverlapGraph::clipBranches(void)
{
	CLOCKSTART;
	UINT64 num_clip_branches = 0;

	UINT64 max_in_ovl, max_out_ovl, ovl;

	for(UINT64 i = 1; i < m_graph->size(); i++){
		if(m_graph->at(i)->size() > 1){
			max_in_ovl = 0; max_out_ovl = 0;
			t_edge_vec inEdges, outEdges;
			vector<UINT64> inOvls, outOvls;
			for(UINT64 j = 0; j < m_graph->at(i)->size(); j++){
				Edge* e = m_graph->at(i)->at(j);
				// Find the first overlap length
				ovl = e->getOverlapLen();
				// Do not consider loop for now TODO might later
				if(!e->isLoop()){
					// In edges
					if(!((e->getOrientation() >> 1) & 1)){
						inEdges.push_back(e);
						inOvls.push_back(ovl);
						if(ovl > max_in_ovl){
							max_in_ovl = ovl;
						}
					}
					// Out edges
					else{
						outEdges.push_back(e);
						outOvls.push_back(ovl);
						if(ovl > max_out_ovl){
							max_out_ovl = ovl;
						}
					}
				}
			}
			if(inEdges.size() > 1){
				for(UINT64 k = 0; k < inEdges.size(); k++){
					if((inOvls.at(k) + minOvlDiffToClip) < max_in_ovl){
//						FILE_LOG(logDEBUG1) << "Break edge " << *(inEdges.at(k)) << " with overlap length " << inOvls.at(k) << " comparing to " << max_in_ovl << "\n";
						t_edge_vec sub_edges = inEdges.at(k)->breakEdge(0);
						removeEdge(inEdges.at(k));
						for(auto it = sub_edges.begin(); it != sub_edges.end(); ++it)
							insertEdge(*it);
						++num_clip_branches;
					}
				}
			}
			if(outEdges.size() > 1){
				for(UINT64 k = 0; k < outEdges.size(); k++){
					if((outOvls.at(k) + minOvlDiffToClip) < max_out_ovl){
//						FILE_LOG(logDEBUG1) << "Break edge " << *(outEdges.at(k)) << " with overlap length " << outOvls.at(k) << " comparing to " << max_out_ovl << "\n";
						t_edge_vec sub_edges = outEdges.at(k)->breakEdge(0);
						removeEdge(outEdges.at(k));
						for(auto it = sub_edges.begin(); it != sub_edges.end(); ++it)
							insertEdge(*it);
						++num_clip_branches;
					}
				}
			}
//			FILE_LOG(logDEBUG1) << "Node " << i << " done with clipBranches" << "\n" << "\n";
		}
	}
	FILE_LOG(logINFO) << "Short overlap branches clipped: " << num_clip_branches << "\n";
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return num_clip_branches;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  reduceLoops
 *  Description:  This function remove loops
 *                a>--->b>--->b>--->c
 *                a<---<b<--->b>--->c
 * =====================================================================================
 */
UINT64 OverlapGraph::reduceLoops(void)
{
	if(this->m_flowComputed == false) {
		return 0;
	}
	CLOCKSTART;
	UINT64 counter = 0, remove_counter = 0;
	Edge *ab,*bb,*bc;
	for(UINT64 i = 1; i < m_graph->size(); i++)
	{
		if(m_graph->at(i)->size() == 4) // only four edges. The loop is counted twice.
		{
			UINT64 loopCount = 0, incomingEdgeCount = 0, outgoingEdgeCount = 0;
			for(UINT64 j = 0; j< m_graph->at(i)->size(); j++)
			{
				if(m_graph->at(i)->at(j)->isLoop()) // This is a loop
				{
					loopCount++;
					bb = m_graph->at(i)->at(j);
				}
				else if(((m_graph->at(i)->at(j)->getOrientation() >> 1) & 1 ) == 0) // incoming edge
				{
					incomingEdgeCount++;
					ab = m_graph->at(i)->at(j)->getReverseEdge();
				}
				else if(((m_graph->at(i)->at(j)->getOrientation() >> 1) & 1) == 1) // outgoing edge
				{
					outgoingEdgeCount++;
					bc = m_graph->at(i)->at(j);
				}
			}
			/* Be aware that the loops without flow are not removed in function removeAllEdgesWithoutFlow.
			 * Therefore if merge an edge with flow and a loop with flow 0, the resulted edge will have flow 0,
			 * which is wrong. Here we will fix it, reassign the correct number of flow to the loop.
			 * If later reduce tree with flow 0 and two other edges, 
			 * then one edge with flow 0 would have been deleted before 
			 * it's contracted with the other branch, wrong!
			 * *** Now do not consider flow any more since it's not reliable with the removal of edges with flow
			 * *** Whenever an edge with flow is removed, it creates an imbalance of flow at its incident nodes.
			 */
			
			// two in the loop and one incoming and one outgoing
			if(loopCount==2 && incomingEdgeCount == 1 && outgoingEdgeCount == 1)
			{  
				if(bb->getOrientation() == 0)
				{
					++counter;
					Edge *new_edge = Add(ab,  bb->getReverseEdge());
					insertEdge(new_edge);
					removeEdge(ab,false);
					removeEdge(bb,false);
				}
				else if(bb->getOrientation() == 3)
				{
					++counter;
					Edge *new_edge = Add(ab, bb);
					insertEdge(new_edge);
					removeEdge(ab);
					removeEdge(bb);
				}
				else{
					FILE_LOG(logDEBUG1) << "Loop has wrong orientation, cannot be used to reduce the graph" << "\n";
					++remove_counter;
					removeEdge(bb);
				}
			}
			/* two in the loop and two incoming
			 * reduce in the case: *--->b>---<b<---* */
			else if (loopCount==2 && incomingEdgeCount == 2 && outgoingEdgeCount == 0 && bb->getOrientation() == 2){ // 
				counter++;
				Edge *new_edge = Add(ab, bb);
				insertEdge(new_edge);
				removeEdge(ab);
				removeEdge(bb);
			}
			/* two in the loop and two incoming
			 * reduce in the case: *---<b<--->b>---* */
			else if (loopCount==2 && incomingEdgeCount == 0 && outgoingEdgeCount == 2 && bb->getOrientation() == 1){ // 
				counter++; 
				Edge *new_edge = Add(bb, bc);
				insertEdge(new_edge);
				removeEdge(bc);
				removeEdge(bb);
			}
			/* Cannot reduce graph */
			else if (loopCount == 2){
				remove_counter++;
				removeEdge(bb);
			}
		}
	}
	if(counter > 0 || remove_counter > 0){
		FILE_LOG(logINFO) <<" Loops reduced: " << counter << "\n"; // Number of loop we were able to reduce
		FILE_LOG(logINFO) <<" Loops removed: " << remove_counter << "\n"; // Number of loop we were able to reduce
	}
	CLOCKSTOP;
	FILE_LOG(logINFO) << "numberOfEdges = " << m_numberOfEdges << "\n";
	return (counter + remove_counter);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculateBoundAndCost
 *  Description:  This function calculates the cost and bounds for an edge in the overlap graph.
 *                This function is very sensitive to the assembled contigs and to the cost function parameters
 * =====================================================================================
 */
void OverlapGraph::calculateBoundAndCost(const Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST)
{
	for(UINT64 i = 0; i < 3; i++)		// For the simple edges we put very high cost
	{
		FLOWLB[i] = 0; FLOWUB[i] = 10; COST[i] = 500000;
	}

	if(edge->getListOfReads() && !edge->getListOfReads()->empty()) // Composite Edge
	{
		// Case1: Composite Edge of at least minFlowReadCountThreshold (default: 20) reads. Must have at least one unit of flow.
		// Case2: Composite Edge length is longer than minFlowEdgeLengthThreshold (default: 1000)
		if(edge->getListOfReads()->size() >= minReadsCountInEdgeToBe1MinFlow || edge->getEdgeLength() > minEdgeLengthToBe1MinFlow)
		{
			s_nGoodEdges++;
			s_nReads_in_goodEdges += edge->getListOfReads()->size();
			// the first 1 flow must be used and has a cost of 1
			// each additional flow up to 8 flows has a cost of 100000

			// this edge carry the first 1 flow
			FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
			// this edge carries the second 1 flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			// this edge provides additional flow after the second flow
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
		else // Short composite edge containing less than 20 reads. May have zero flow.
		{
			// the first 1 flow may not be required, but has a low cost of 1
			// each additional flow up to 8 flows has a cost of 100000

			// this edge carries the first 1 flow
			FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1;
			// this edge carries the second unit of flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			// this edge provides additional flow after the two units of flow.
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  findEdge
 *  Description:  Return edge with most reads on it, between source and destination
 * =====================================================================================
 */
Edge* OverlapGraph::findEdge(const UINT64 & source, const UINT64 & destination)
{
	t_edge_vec edges = findEdges(source, destination);
	if(edges.empty())
		return nullptr;
	else
		return edges.front();
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  findEdges
 *  Description:  Return all edges between source and destination, sorted by number of reads
 *  		  in decreasing order.
 * =====================================================================================
 */
t_edge_vec OverlapGraph::findEdges(const UINT64 & source, const UINT64 & destination)
{
	t_edge_vec edges;
	for(UINT64 i = 0; i < m_graph->at(source)->size(); i++) // For the list of edges of the source node.
	{
		// check if there is an edge to destination
		if(m_graph->at(source)->at(i)->getDestinationRead()->getReadID() == destination){
			edges.push_back(m_graph->at(source)->at(i));	// return the edge.
		}
	}
	if(edges.empty()){
		FILE_LOG(logDEBUG1) << "Cannot find edge from " << source << " to " << destination << "\n";
	}
	// Sort the edges by number of reads contained in decreasing order
	else
		sort(edges.begin(),edges.end(),compareEdgesByReads);
	return edges;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  merge_edges_with_IDlists
 *  Description:  Given lists of edge IDs (positive and negative), merge the corresponding
 *                edges in each list, and insert them in the graph.
 *                Also delete these edges at the end.
 * =====================================================================================
 */
UINT64 OverlapGraph::merge_edges_with_IDlists(const vector<vector<INT64>> *lists_edgeIDs_to_merge,
		const map<INT64, Edge*> & id_to_edge)
{
	CLOCKSTART;
	UINT64 num_merges = 0;
	// Edges that will be merged, and then deleted in the end
	set<Edge*> edges_to_remove;
	// list of edges to merge into one edge
	t_edge_vec edges_to_merge;

	Edge *e = nullptr;
	// For each vecotr of edges, merge them into one edge
	for(UINT64 i = 0; i < lists_edgeIDs_to_merge->size(); ++i){
		FILE_LOG(logDEBUG1) << i << "th list in the lists" << "\n";
		edges_to_merge.clear();
		const vector<INT64> &current_id_list = lists_edgeIDs_to_merge->at(i);
		for(UINT64 j = 0; j < current_id_list.size(); ++j){
			INT64 edge_id = current_id_list.at(j);
			FILE_LOG(logDEBUG1) << "\t";
			// id to edge map only has postive edge ids as keys
			// So negative edge ids need to be accessed through their
			// reverse edges
			if(edge_id > 0)
				e = id_to_edge.at(current_id_list.at(j));
			else
				e = id_to_edge.at(-current_id_list.at(j))->getReverseEdge();
			FILE_LOG(logDEBUG1) << *e << "\n";
			edges_to_merge.push_back(e);
			// the negative edges' twin edges will be inserted
			if(edge_id > 0)
				edges_to_remove.insert(e);
			else
				edges_to_remove.insert(e->getReverseEdge());
		}
		Edge *merged_edge = mergeEdges(edges_to_merge);
		++num_merges;
		insertEdge(merged_edge);
	}
	FILE_LOG(logDEBUG) << "Total edges to remove: " << edges_to_remove.size() << "\n";
	for(auto it = edges_to_remove.cbegin(); it != edges_to_remove.cend(); ++it){
		FILE_LOG(logDEBUG1) << "Now remove edge " << *(*it) << "\n";
		removeEdge(*it);
	}
	CLOCKSTOP;
	return num_merges;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  readThreading
 *  Description:  Compare the Illumina read tiling in a Pacbio read to the edges in the
 *  		  assembly graph
 * =====================================================================================
 */
UINT64 OverlapGraph::readsThreading(ifstream & inputSamStream)
{
	CLOCKSTART;
	UINT64 edgesThreaded(0);

	// Initialize the pacbio read name to NA,
	// it'll be an indicator of a new pacbio read 
	// in the sam file streaming
	std::string pacbio_name("na"), pacbio_read_name;
	UINT64 readName, p_read;
	
	UINT64 orientation;
	int flag;
	
	// each vector element is a list of edges to merge
	vector< t_edge_vec > list_edges_to_merge;

	// list of read (numbers) tiling the pacbio read
	vector<UINT64> readNumbers;

	// sample SAM record, separated by \t
	// 385350	0	2219562	924	255	360M	
	// *	0	0	AACGTATCTGATTTCCCCAGGAAATTTC
	//                      IIIIIIIIIIIIIIIIIIIIIIIIIIII
	// AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	
	// MD:Z:360	YT:Z:UU
	while(!inputSamStream.eof()){
		inputSamStream >> readName >> flag >> pacbio_read_name;
		inputSamStream.ignore(numeric_limits<streamsize>::max(),'\n');

		// Condense the read number (name) and orientation in one integer
		p_read = (readName << 1);

		// flag 0x10 indicates reverse complement of read is aligned
		orientation = ( ((flag & 0x10) == 0x10) ? 0 : 1);

		p_read += orientation;

		// A new pacbio read started
		if(pacbio_read_name != pacbio_name){
			// Print the tiling for the PacBio read
			if(readNumbers.size() > 0){
				FILE_LOG(logDEBUG1) << "\nPacbio read: " << pacbio_name  
					<< ", " << readNumbers.size() << " reads ### "
					<< "\n";
				if(loglevel > 2){
					for(size_t i = 0; i < readNumbers.size(); i++)
						//TODO If Illumina reads are reprinted for bowtie alignment,
						//we do not need readIDMap to get the read IDs in the graph TODO
						FILE_LOG(logDEBUG1) << m_dataset->getReadIDMap()->at(readNumbers.at(i) >> 1) << " ";
					FILE_LOG(logDEBUG1) << "\n";
				}
				readThreading(readNumbers, list_edges_to_merge);
			}
			pacbio_name = pacbio_read_name;
			readNumbers.clear();
		}
		// If not a new pacbio read, grow the vector of Illumina reads
		readNumbers.push_back(p_read);
	}
	// Print the tiling for the last PacBio read, and do threading with this PacBio read
	if(loglevel > 2){
		FILE_LOG(logDEBUG1) << "\nPacbio read: " << pacbio_name  
			<< " " << readNumbers.size() << " reads" << "\n";
		for(size_t i = 0; i < readNumbers.size(); i++)
			FILE_LOG(logDEBUG1) << m_dataset->getReadIDMap()->at(readNumbers.at(i) >> 1) << "\t";
		FILE_LOG(logDEBUG1) << "\n";
	}
	readThreading(readNumbers, list_edges_to_merge);

	// Set edge ID, used for constructing overlap graph with these edges as nodes
	FILE_LOG(logDEBUG) << "Set edge IDs for the mergeable edges." << "\n";
	INT64 edge_id = 0;
	for(size_t i = 0; i < list_edges_to_merge.size(); i++){
		if (list_edges_to_merge.at(i).size() > 1){
			for(size_t j = 0; j < list_edges_to_merge.at(i).size(); j++){
				Edge *e = list_edges_to_merge.at(i).at(j);
				if(e->getEdgeID() == 0){
					e->setEdgeID(++edge_id);
				}
				FILE_LOG(logDEBUG1) << *e << "\n";
			}
		}
		FILE_LOG(logDEBUG1) << "\n";
	}

	edgesThreaded = mergeListOfEdges(list_edges_to_merge);
	
	CLOCKSTOP;
	return edgesThreaded;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  locateRead
 *  Description:  Find edge that has read on it.
 *  		  p_read: (readID in the file << 1) + orientation (in last bit)
 * =====================================================================================
 */
Edge * OverlapGraph::locateRead(UINT64 p_read) const
{
	vector<t_edge_loc_pair> *r_list_edges;
	UINT64 read_number = m_dataset->getReadIDMap()->at(p_read >> 1);
	Read * r = m_dataset->at(read_number);
	UINT8 read_orient = (p_read & 0x1);
	if(read_orient == 1)
		r_list_edges = r->getFwdEdges();
	else
		r_list_edges = r->getBwdEdges();
	if(r_list_edges->size() > 0)            /* If this read is contained in some edge */
		return r_list_edges->at(0).first;
	else
		return nullptr;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  mergeListOfEdges
 *  Description:  Merge the edges that are supported by the long read threading.
 *                This will involve the new module, using overlap graph.
 * =====================================================================================
 */
UINT64 OverlapGraph::mergeListOfEdges(const vector< t_edge_vec > & list_edges_to_merge)
{

	UINT64 num_merges = 0;
	tDataSet<tRead<INT64> > *e_data = new tDataSet< tRead<INT64> >;
	map<INT64, Edge*> id_to_edge;

	for(size_t i = 0; i < list_edges_to_merge.size(); i++){
		vector<INT64> v;
		for(size_t j = 0; j < list_edges_to_merge.at(i).size(); j++){
			Edge *e = list_edges_to_merge.at(i).at(j);
			INT64 id = e->getEdgeID();
			v.push_back(id);
			if (id > 0){
				id_to_edge[id] = e;
			}
		}
		// If PacBio read can only successfully thread one edge, then no need to use it
		if(v.size() > 1){
			tRead<INT64> e_read(v);
			e_data->addRead(e_read);
		}
	}
	// Check the status of the map from edge IDs to the edges
	FILE_LOG(logDEBUG) << "Number of pairs in the map from edge ID to edge address: " << id_to_edge.size() << "\n";
	if(loglevel > 3){
		for(auto it = id_to_edge.cbegin(); it != id_to_edge.cend(); ++it){
			FILE_LOG(logDEBUG) << setw(5) << it->first << ": " << it->second << " --> " 
				<< *(it->second) << "\n";
		}
	}
	// Remove duplicated list of edges
	e_data->rmDupReads();
	if (loglevel > 3){
		FILE_LOG(logDEBUG1) << "reads in the data set are sorted? " 
			<< std::boolalpha << e_data->isSorted() << "\n";
		FILE_LOG(logDEBUG1) << "data does not have duplicate reads? " 
			<< std::boolalpha << e_data->isDupRemoved() << "\n";
		FILE_LOG(logDEBUG1) << *e_data << "\n";
	}

	// Construct overlap between the lists of edges with 1 edge overlap
	tOverlapGraph<INT64> * e_graph = new tOverlapGraph<INT64>(e_data, 1);
	vector<tEdge<tRead<INT64>>*> contigEdges;
	e_graph->getEdges(contigEdges);
	if (loglevel > 2){
		e_graph->printGraph(outputFilenamePrefix + "_edge_overlap.gdl", contigEdges);
	}

	vector<vector<INT64>> *lists_edgeIDs_to_merge = e_graph->getContigs(contigEdges, true);
	num_merges = merge_edges_with_IDlists(lists_edgeIDs_to_merge, id_to_edge); 

	delete e_graph;
	delete e_data;
	delete lists_edgeIDs_to_merge;

	return num_merges;
}

//=============================================================================
// Default constructor
//=============================================================================
OverlapGraph::OverlapGraph(void)
{
	m_dataset	= nullptr;
	m_graph		= nullptr;
	m_deletedEdges 	= nullptr;

	m_numberOfNodes	= 0;
	m_numberOfEdges	= 0;
	m_minOvl 	= 0;
	m_maxSubs	= 0;
	m_maxEdits	= 0;
	m_flowComputed 	= false;
}


OverlapGraph::OverlapGraph(const vector<std::string> &edge_files, 
		const vector<std::string> &read_files, 
		const UINT64 &minOvl = 0, const UINT64 &maxSubs = 0, const UINT64 &maxEdits = 0)
	: m_minOvl(minOvl), m_maxSubs(maxSubs), m_maxEdits(maxEdits)
{
	CLOCKSTART;
	m_numberOfNodes	= 0;
	m_numberOfEdges	= 0;
	m_flowComputed 	= false;

	/*  JJ: we only need Edge file to get all the reads and edges,
	 *  at least the edges that have edges connected to them.
	 *  However, then we need another data structure to keep track of
	 *  which reads are already loaded in the m_dataset, which can be
	 *  expensive.
	 *  When  with a better way is come up, we should have a function that only reads
	 *  edges files.
	 */

	m_dataset = new DataSet(read_files, true, false); // construct dataset from reads file(s)

	FILE_LOG(logINFO) << "Total number of unique reads loaded from read file(s): " 
		<< m_dataset->size() << "\n";

	UINT64 numOfUniqueRead(m_dataset->size());

	m_graph = new vector< t_edge_vec* >;

	m_graph->reserve(numOfUniqueRead + 1);

	for(UINT64 i = 0; i <= numOfUniqueRead; i++) {
		t_edge_vec *newList = new t_edge_vec;
		m_graph->push_back(newList);
	}

	FILE_LOG(logDEBUG1) << "graph has " << m_graph->size() << " unique number of reads.\n";
	// loop edgeFilenameList
	for (vector<std::string>::const_iterator it = edge_files.begin(); 
			it != edge_files.end(); ++it) {
		loadEdgesFromEdgeFile(*it);
	}
	FILE_LOG(logERROR) << "numberOfEdges loaded from edge file(s): " << m_numberOfEdges << "\n";

	if(loglevel > 4){
		string idmap_file = "idmap.txt";
		ofstream out(idmap_file.c_str());
		printUnorderedMap(*(m_dataset->getReadIDMap()), out);
		out.close();
	}
	// Sort all its incident edges according to the overlap offset 
	// Useful in bubble popping
	FILE_LOG(logINFO) << "Sort edges incident to each node by destination read ID.\n";
	sortEdgesByDestID();

	FILE_LOG(logINFO) << "Initial simplification: contract composite edges, remove dead end nodes,"
		<< " and clip branches with very short overlap length.\n";
		
	// Composite edge contraction with remove dead end nodes
	UINT64 counter(0);
	do {
		counter = contractCompositeEdges();
		counter += removeDeadEndNodes();
	} while (counter > 0);
	FILE_LOG(logERROR) << "numberOfEdges = " << m_numberOfEdges << "\n";
	/* disconnect the edges incident to nodes and have small overlap lengths */
	clipBranches();
	CLOCKSTOP;
}
//=============================================================================
// Default destructor
//=============================================================================
OverlapGraph::~OverlapGraph()
{
	// Free the memory used by the overlap graph.
//	if (m_deletedEdges != nullptr){
//		for(UINT64 i = 0 ; i < m_deletedEdges->size(); ++i){
//			if (m_deletedEdges->at(i) != nullptr){
//				delete m_deletedEdges->at(i);
//				m_deletedEdges->at(i) = nullptr;
//			}
//		}
//		delete m_deletedEdges;
//		m_deletedEdges = nullptr;
//	}
	if (m_graph!= nullptr){
		for(UINT64 i = 0; i < m_graph->size(); i++) {
			if (m_graph->at(i) != nullptr){
				for(UINT64 j = 0; j< m_graph->at(i)->size(); j++) {
					if (m_graph->at(i)->at(j) != nullptr){
						delete m_graph->at(i)->at(j);
						m_graph->at(i)->at(j) = nullptr;
					}
				}
				delete m_graph->at(i);
				m_graph->at(i) = nullptr;
			}
		}
		delete m_graph;
		m_graph = nullptr;
	}
	if (!m_dataset)
		delete m_dataset;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sortEdgesByLength
 *  Description:  For each node, sort its incident edges by their lengths in increasing order
 * =====================================================================================
 */
void OverlapGraph::sortEdgesByDestID()
{
	for(UINT64 i = 1; i < m_graph->size(); i++) {
		if(!m_graph->at(i)->empty()){
			sort(m_graph->at(i)->begin(), m_graph->at(i)->end(), compareEdgesByDestID);
		}
	}
}
void OverlapGraph::sortEdgesByLength()
{
	for(UINT64 i = 1; i < m_graph->size(); i++) {
		if(!m_graph->at(i)->empty()){
			sort(m_graph->at(i)->begin(), m_graph->at(i)->end(), compareEdgesByLength);
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  simplifyGraph
 *  Description:  simplify graph iteratively
 * =====================================================================================
 */
void OverlapGraph::simplifyGraph(void)
{
	UINT64 counter = 0;
	do
	{
		counter = contractCompositeEdges();
		counter += removeSimilarEdges();
		counter += removeDeadEndNodes();
		// the three steps below requires flow to be computed
		// if simplifyGraph is called in the unitig graph, these two functions will just return.
		counter += removeShortBranches();	// Remove short branches
		counter += reduceLoops();	// Reduce loops

	} while (counter > 0);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  calculateFlowStream
 *  Description:  Calculate min cost flow
 *  An sample input to the CS2 algorithm
 * 
 *  p min       3840      13449	// p min numberOfNodes numberOfEdges
 *  n          1         0	// initial flow of 0 in node 1, node 1 is the supersource
 *  n       3840         0	// initial flow of 0 in node 3840, which is the supersink (equal to the number of nodes)
 *  // edge from supersink to supersource, LowBound(1), UpperBound(1000000), Cost per unit of flow(1000000)
 *  a       3840          1          1    1000000    1000000	
 *  // edge from supersource to node 2, with the defined cost function
 *  a          1          2          0    1000000          0	
 *
 *  connect each node to supersource and supersink
 *  connect every edge in the original graph
 *
 * =====================================================================================
 */
void OverlapGraph::calculateFlowStream(void) 
{
	CLOCKSTART;
	// Add super source and super sink nodes, add edge from super sink to super source with very big cost
	// Add edge from super source to every node in the graph, also from every node in the graph to the super sink
	// Every edge will be assigned a lower bound and an upper bound of the flow (capacity), and the cost associated with the edge

	UINT64 V = m_numberOfNodes*2 + 2, E = m_numberOfEdges * 3 + m_numberOfNodes * 4 + 1 , SUPERSOURCE = 1, SUPERSINK = V;

	// Flow bounds and cost of the edges, 
	// cost function originally is a piecewise function with 3 segments
	INT64 FLOWLB[3], FLOWUB[3], COST[3];            
	stringstream ss;
	// Problem description: Number of nodes and edges in the graph.
	ss << "p min " << setw(10) << V << " " << setw(10) << E << "\n";    

	// Flow in the super source
	ss << "n " << setw(10) << SUPERSOURCE << setw(10) << " 0" << "\n";  

	// Flow in the super sink.
	ss << "n " << setw(10) << SUPERSINK << setw(10) << " 0" << "\n";    

	// Add an edge from super sink to super source with very high cost (almost infinity), also at most can be used once
	FLOWLB[0] = 1; 
	FLOWUB[0] = std::numeric_limits<UINT64>::max(); 
	COST[0]   = 1000000;
	ss << "a " << setw(10) << SUPERSINK << " " << setw(10) << SUPERSOURCE 
		<< " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] 
		<< " " << setw(10) << COST[0] << "\n"; 


	// If the ID of a node in the original graph is 100 and directed graph is 5
	// Then listOfNodes->at(100) is equal to 5
	// and ListOfNodesReverse->at(5) is equal to 100.
	vector<UINT64> *listOfNodes = new vector<UINT64>;
	vector<UINT64> *listOfNodesReverse = new vector<UINT64>;

	// For n nodes in the graph, CS2 requires that the nodes are numbered from 1 to n. 
	// In the overlap graph, the nodes does not have sequencinal ID. We need to convert them to 1 - n
	for(UINT64 i = 0; i <= m_graph->size(); i++)      {
		listOfNodes->push_back(0);
		listOfNodesReverse->push_back(0);
	}

	// This loop set lower bound and upper bound from super source to each node, and from each node to super sink. All costs are 0.
	UINT64 currentIndex = 1;
	for(UINT64 i = 1; i < m_graph->size(); i++)
	{
		// edges to and from the super source and super sink
		if(!m_graph->at(i)->empty())
		{
			// Mapping between original node ID and cs2 node ID
			listOfNodes->at(i) = currentIndex;                                      

			// Mapping between original node ID and cs2 node ID
			listOfNodesReverse->at(currentIndex) = i;                       
			FLOWLB[0] = 0; FLOWUB[0] = 1000000; COST[0] = 0;
			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex << " " 
				<< setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex + 1 << " " 
				<< setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

			ss << "a " << setw(10) << 2 * currentIndex << " " << setw(10) << SUPERSINK << " " 
				<< setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";
			
			ss << "a " << setw(10) << 2 * currentIndex + 1 << " " << setw(10) << SUPERSINK << " " 
				<< setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

			currentIndex++;
		}
	}

	// This loop converts the original bi-directed edges to directed edges (1 becomes 6).
	for(UINT64 i = 1; i < m_graph->size(); i++) {
		// edges to and from the super source and super sink
		if(!m_graph->at(i)->empty()) 
		{
			for(UINT64 j = 0; j < m_graph->at(i)->size(); j++)
			{
				Edge *edge = m_graph->at(i)->at(j);
				UINT64 u = listOfNodes->at(edge->getSourceRead()->getReadID());
				UINT64 v = listOfNodes->at(edge->getDestinationRead()->getReadID());

				if(u < v || (u == v && edge < edge->getReverseEdge())) {
					calculateBoundAndCost(edge, FLOWLB, FLOWUB, COST);
					// Here for each edge we add three edges with different values of cost and bounds.
					// Total 6 edges considering the reverse edge too.
					// For details on how to convert the edges off different types please see my thesis.
					UINT64 u1 = 2 * u, u2 = 2 * u + 1, v1 =  2 * v, v2 = 2 * v + 1;
					if(edge->getOrientation() == 0)
					{
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << 
							FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << 
							FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << 
							FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << 
							FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";

						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << 
							FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << 
							FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
					}
					else if(edge->getOrientation() == 1)
					{
						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << 
							FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << 
							FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << 
							FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << 
							FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";

						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << 
							FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << 
							FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";

					}
					else if(edge->getOrientation() == 2)
					{
						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << 
							FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << 
							FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << 
							FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << 
							FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";

						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << 
							FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << 
							FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";

					}
					else if(edge->getOrientation() == 3)
					{
						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << 
							FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << 
							FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << "\n";

						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << 
							FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << 
							FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << "\n";

						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << 
							FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << 
							FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << "\n";
					}
				}
			}
		}
	}
	FILE_LOG(logINFO) <<"Finished initializing flow to edges.\n";

	FILE_LOG(logDEBUG) << "Number of edges with flow 1 set is " << OverlapGraph::s_nGoodEdges << "\n";
	FILE_LOG(logDEBUG) << "Number of reads contained in these edges is " << OverlapGraph::s_nReads_in_goodEdges << "\n";
	FILE_LOG(logINFO) << "Calling CS2 for flow analysis\n";
	stringstream oss;
	main_cs2(&ss, oss);
	FILE_LOG(logINFO) << "Flow analysis finished\n";

	if(loglevel > 2){
		std::string flowfile = outputFilenamePrefix + "_init.flow";
		FILE_LOG(logDEBUG) << "Print result flow in graph to " << flowfile << "\n";
		ofstream flowin(flowfile.c_str());
		flowin << ss.str();
		flowin.close();

		flowfile = outputFilenamePrefix + "_calc.flow";
		FILE_LOG(logDEBUG) << "Print result flow in graph to " << flowfile << "\n";
		ofstream flowout(flowfile.c_str());
		flowout << oss.str();
		flowout.close();
	}

	FILE_LOG(logINFO) <<"Start assigning calculated flow to edges.\n";
	std::string s, d, f;
	UINT64 lineNum = 0;

	while(!oss.eof()) {
		lineNum++;
		UINT64 source, destination, flow;
		// get the flow from CS2
		oss >> source >> destination >> flow;	

		// Map the source to the original graph
		UINT64 mySource = listOfNodesReverse->at(source/2);	

		// Map the destination in the original graph
		UINT64 myDestination = listOfNodesReverse->at(destination/2);	

		if(source != SUPERSINK && source != SUPERSOURCE && destination != SUPERSOURCE && destination != SUPERSINK && flow!=0)
		{
			Edge *edge = findEdge(mySource, myDestination);	// Find the edge in the original graph.
			if(edge){
				edge->m_flow += flow;	// Add the flow in the original graph.
				edge->getReverseEdge()->m_flow += flow;	// Also add the same flow to the twin edge to make sure that the flow is the same for the twin edges (flow is doubled if both directions have same flow)
			}
			else{
				FILE_LOG(logDEBUG1) << "edge does not exist in graph!\n";
			}
		}
	}
	delete listOfNodes;
	delete listOfNodesReverse;
	this->m_flowComputed = true;

	CLOCKSTOP;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  removeAllEdgesWithoutFlow
 *  Description:  After flow analysis, remove all edges with flow 0
 * =====================================================================================
 */
UINT64 OverlapGraph::removeAllEdgesWithoutFlow()
{
	if(!m_flowComputed)
		return 0;
	CLOCKSTART;
	UINT64 num_edge_rm(0);
	for(UINT64 i = 1; i < m_graph->size(); i++) // For each read.
	{
		if(!m_graph->at(i)->empty())	// If the read has some edges.
		{
			for(UINT64 j=0; j < m_graph->at(i)->size(); j++) // For each edge
			{
				Edge * edge = m_graph->at(i)->at(j);

				//Also remove loops without flow. This means, by default, 
				//the edges formed by loops that do not contain many reads will not be used.
				if(edge->m_flow == 0 && !edge->isLoop()) {
					removeEdge(edge,false);
					++num_edge_rm;
				}
			}
		}
	}
	FILE_LOG(logINFO) <<"No flow edges removed: " << num_edge_rm << "\n";
	CLOCKSTOP;
	return num_edge_rm;
}

//=============================================================================
// Function to compare two edges by the destination read number. 
// The read number was determined by the read std::string lexicographically.
// Used for sorting.
//=============================================================================
bool compareEdgesByDestID (const Edge *edge1, const Edge* edge2)
{
	if (edge1->getDestinationRead()->getReadID() < edge2->getDestinationRead()->getReadID())
	{
		return true;
	}
	else if (edge1->getDestinationRead()->getReadID() == edge2->getDestinationRead()->getReadID()){
		return (compareEdgesByLength(edge1, edge2));
	}
	else
		return false;
}


//=============================================================================
// Function to compare two edges by overlap offset. 
// Used for transitive edge removal (which is not included any more)
//=============================================================================
bool compareEdgesByLength (const Edge *edge1, const Edge* edge2)
{
	return (edge1->getEdgeLength() < edge2->getEdgeLength());
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  compareEdgesByReads
 *  Description:  Compare edges by number of reads contained in the them
 * =====================================================================================
 */
bool compareEdgesByReads (const Edge *edge1, const Edge* edge2)
{
	UINT64 num_reads1 = (edge1->getListOfReads() ? edge1->getListOfReads()->size() : 0);
	UINT64 num_reads2 = (edge2->getListOfReads() ? edge2->getListOfReads()->size() : 0);
	return (num_reads1 > num_reads2);
}


//=============================================================================
// Function to compare two std::string length
//=============================================================================
bool compareStringByLength (const std::string & seq1, const std::string & seq2)
{
	return (seq1.length() < seq2.length());
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  loadEdgesFromEdgeFile
 *  Description:  Load edge file
 *                15	1653045	3,90,0,0,370,280,369,399,0,89,NA
 *                15	3610522	1,169,0,0,370,201,369,414,0,168,NA
 *                9	95931	0,127,0,0,353,226,352,365,-12,114,NA
 *                8	950354	0,54,0,0,399,345,398,400,-1,52,NA
 *                source destination orientation, overlap length, substitution, edits,
 *			length1, start1, stop1, length2, start2, stop2, err info
 * =====================================================================================
 */
void OverlapGraph::loadEdgesFromEdgeFile(const std::string &edgeFilename)
{
	CLOCKSTART;
	FILE_LOG(logINFO) << "Load edge file: " << edgeFilename << "\n";

	// Open file
	ifstream filePointer;
	filePointer.open(edgeFilename.c_str());
	if(!filePointer.is_open() )
		MYEXIT("Unable to open file: "+edgeFilename);

	UINT64 edgeNumber(0);
	// Read file
	std::string line_text;
	while(getline(filePointer,line_text)) {
		// FILE_LOG(logDEBUG1) << "Edge " << line_text << "\n";
		// 15	3610522	1,169,0,0,370,201,369,414,0,168,NA
		vector<std::string> line_elements;
		std::stringstream text_ss(line_text);
		std::string element;

		// get the three parts separated by tabs in each line
		while (getline(text_ss, element, '\t')) {
			line_elements.push_back(element);
		}

		// Get source readID in the edge file (corresponding to the ID in reads file),
		// then use the readIDMap to find the ID in the graph
		UINT64 source_org, source;
		istringstream t1(line_elements.at(0)); 
		t1 >> source_org;
		// Find source readID from readIDMap
		t_idmap::const_iterator got_s = m_dataset->getReadIDMap()->find (source_org);
		if ( got_s == m_dataset->getReadIDMap()->end() ) {
			MYEXIT("not found");
		}
		else {
			source = got_s->second;
		}

		// Do the same for the destination ID
		UINT64 destination_org, destination;	// Tab delimited Column2: destination read ID
		istringstream t2(line_elements.at(1)); t2 >> destination_org;
		// Find destination readID from readIDMap
		t_idmap::const_iterator got_d = m_dataset->getReadIDMap()->find (destination_org);
		if ( got_d == m_dataset->getReadIDMap()->end() ) {
			MYEXIT("not found");
		}
		else {
			destination = got_d->second;
		}

		// Properties
		std::string properties = line_elements.at(2);	// Tab delimited Column3: edge attributes

		// For properties list
		vector<std::string> propertiesList;
		std::stringstream properties_ss(properties);
		while (getline(properties_ss, element, ',')) {
			propertiesList.push_back(element);
		}

		// 15	3610522	1,169,0,0,370,201,369,414,0,168,NA
		// properties
		UINT8 orientation;	// Property Col1: orientation
		int dec_orient;
		std::istringstream p1(propertiesList.at(0)); p1 >> dec_orient;
		orientation = static_cast<UINT8>(dec_orient);

		UINT32 overlapLength;	// Property Col2: overlap length
		std::istringstream p2(propertiesList.at(1)); p2 >> overlapLength;

		UINT64 substitutions;	// Property Col3: substitutions
		std::istringstream p3(propertiesList.at(2)); p3 >> substitutions;

		UINT64 edits;	// Property Col4: edit distance
		std::istringstream p4(propertiesList.at(3)); p4 >> edits;

		// If the edge has overlap that satisfies the requirements
		if (overlapLength >= m_minOvl && substitutions <= m_maxSubs && 
				edits <= m_maxEdits){
			UINT32 length1;	// Property Col5: read1 length
			std::istringstream p5(propertiesList.at(4)); p5 >> length1;

			UINT32 start1;	// Property Col6: read1 overlap start
			std::istringstream p6(propertiesList.at(5)); p6 >> start1;

//			UINT32 stop1;	// Property Col7: read1 overlap end
//			std::istringstream p7(propertiesList.at(6)); p7 >> stop1;

			UINT32 length2;	// Property Col8: read2 length
			std::istringstream p8(propertiesList.at(7)); p8 >> length2;

//			UINT32 start2;	// Property Col9: read2 overlap start
//			std::istringstream p9(propertiesList.at(8)); p9 >> start2;
//
//			UINT32 stop2;	// Property Col10: read2 overlap length
//			std::istringstream p10(propertiesList.at(9)); p10 >> stop2;

			// mismatch pair<index of edge from the contained read, position>
			t_vpair *mismatchesInput = getMismatchFromString(propertiesList.at(10));

			/*  get overlap offset
			 *  Example edge list:
			 *  0 = u<-----------<v		reverse of u to reverse of v
			 *  1 = u<----------->v		reverse of u to forward of v
			 *  2 = u>-----------<v		forward of u to reverse of v
			 *  3 = u>----------->v		forward of u to forward of v
			 * 1	2   2,34,0,0,35,1,34,35,0,33,NA
			 * 2	3   0,33,0,0,35,2,34,35,0,32,NA
			 * 3	4   3,34,1,1,35,1,34,35,0,33,1CA
			 */

			// overlap offset
			UINT32 overlapOffset = start1;	// correct, but not ready yet JJ: why not ready yet?
			// Insert an edge
			Read *read1 = m_dataset->at(source);
			if(read1->getReadLength() != length1){
				FILE_LOG(logERROR) << "ERROR: Read number " << source 
					<< " has length difference "
					<< read1->getReadLength() << " from reads file and "
					<< length1 << " from edge file\n";
			}
			Read *read2 = m_dataset->at(destination);
			if(read2->getReadLength() != length2){
				FILE_LOG(logERROR) << "ERROR: Read number " << destination 
					<< " has length difference "
					<< read2->getReadLength() << " from reads file and "
					<< length2 << " from edge file\n";
			}
			Edge *an_edge = new Edge(read1,read2,orientation,overlapOffset,mismatchesInput);
			an_edge->make_reverseEdge();
			insertEdge(an_edge);	// insert edge
			++edgeNumber;
			if(edgeNumber % 1000000 == 0){
				FILE_LOG(logINFO) << setw(10) << (edgeNumber / 1000000) << ",000,000"  
					<< " edges loaded to memory, "
					<< setw(7) << checkMemoryUsage() << " MB" << "\n";
			}
		}
	}
	filePointer.close();

	CLOCKSTOP;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getEdges
 *  Description:  save all the edges in the graph in a vector of pointers to these edges
 * =====================================================================================
 */
void OverlapGraph::getEdges(t_edge_vec & contigEdges) const
{
	CLOCKSTART;
	contigEdges.clear();
	for(UINT64 i = 1; i<= m_dataset->size(); i++)
	{
		if(!m_graph->at(i)->empty()) // if this read has some edge(s) going out of it (since now the graph is directed)
		{
			for(UINT64 j = 0; j < m_graph->at(i)->size(); j++)
			{
				Edge * e = m_graph->at(i)->at(j);
				if(e->isSmallerEdge()){  /* Only consider the edges between non-contained reads */
					contigEdges.push_back(e); // List of contigs.
				}
			}
		}
	}
	FILE_LOG(logINFO) << "Number of edges(contigs) in the graph to print: " << contigEdges.size() << "\n";
//	if (m_deletedEdges && m_deletedEdges->size() > 0){
//		for (UINT64 i = 0; i < m_deletedEdges->size(); ++i){
//			Edge * e = m_deletedEdges->at(i);
//			if(e->isSmallerEdge() ){  /* Only consider the edges between non-contained reads */
//				contigEdges.push_back(e); // List of contigs.
//				e->setContigEdgeID(contigEdgeID);
//				++contigEdgeID;
//			}
//		}
//	}
	// Sort the contigs by lengths in decreasing order
	std::sort(contigEdges.begin(), contigEdges.end(), compareEdgesByLength);
	std::reverse(contigEdges.begin(), contigEdges.end());
	FILE_LOG(logINFO) << "Number of edges(contigs) including deleted to print: " << contigEdges.size() << "\n";
	CLOCKSTOP;
}

/*===== FUNCTION: operator<< ==================================================
 *  This function prints the overlap graph in overlap_graph->gdl file. The graph can be viewed by
 *  aisee (free software available at http://www.aisee.com/)
 *  It also stores the contigs in a file.
 *
 *graph: {
 *layoutalgorithm :forcedir
 *fdmax:704
 *tempmax:254
 *tempmin:0
 *temptreshold:3
 *tempscheme:3
 *tempfactor:1.08
 *randomfactor:100
 *gravity:0.0
 *repulsion:161
 *attraction:43
 *ignore_singles:yes
 *node.fontname:"helvB10"
 *edge.fontname:"helvB10"
 *node.shape:box
 *node.width:80
 *node.height:20
 *node.borderwidth:1
 *node.bordercolor:31
 *node: { title:"43" label: "43" }	// node, Title and label are both node ID 43 (and read number)
 *node: { title:"65" label: "65" }
 *............................................
 *............................................
 * edges from source node 43 to destination node 32217, thickness of 3 means composite edge, thickness of 1 for simple edge
 * edge type of backarrowstyle:solid arrowstyle:solid color: green is >----------------<
 * edge type of arrowstyle:solid color: red is <----------------<
 * edge type of arrowstyle: none color: blue  is <------------------->
 * (1,0x,206,30) means (Flow, coverageDepth, OverlapOffset, numberOfReads)
 *
 *edge: { source:"43" target:"32217" thickness: 3 backarrowstyle:solid arrowstyle:solid color: green label: "(1,0x,206,30)" }
 *edge: { source:"65" target:"38076" thickness: 3 arrowstyle:solid color: red label: "(0,0x,75,11)" }
 *edge: { source:"280" target:"47580" thickness: 3 arrowstyle: none color: blue label: "(0,0x,123,11)" }
 *}
 *=============================================================================
 */
ostream& operator<< (ostream &out, const OverlapGraph & graph)
{
	CLOCKSTART;
	t_edge_vec contig_edges;
	graph.getEdges(contig_edges);

	UINT64 thickness;
	// Graph specification before the nodes and edges
	out << "graph: {\n\
		layoutalgorithm :forcedir\n\
		fdmax:704\n\
		tempmax:254\n\
		tempmin:0\n\
		temptreshold:3\n\
		tempscheme:3\n\
		tempfactor:1.08\n\
		randomfactor:100\n\
		gravity:0.0\n\
		repulsion:161\n\
		attraction:43\n\
		ignore_singles:yes\n\
		node.fontname:\"helvB10\"\n\
		edge.fontname:\"helvB10\"\n\
		node.shape:box\n\
		node.width:80\n\
		node.height:20\n\
		node.borderwidth:1\n\
		node.bordercolor:31\n";

	// All the nodes, title and label
	for(UINT64 i = 1; i< graph.m_dataset->size(); i++)
	{
		if(!graph.m_graph->at(i)->empty())
			// Print nodes even if there is some edge incident to it
			out << "node: { title:\""<< i <<"\" label: \"" << i << "\" }\n";	
	}

	// All the edges
	for (UINT64 i = 0; i < contig_edges.size(); i++)
	{
		Edge *e = contig_edges.at(i);
		int num_reads = 0;
		if(e->getListOfReads())
			num_reads = e->getListOfReads()->size();
		thickness = (num_reads > 0) ? 1: 3;
		UINT64 source = e->getSourceRead()->getReadID(), destination = e->getDestinationRead()->getReadID();
		// Edge label: (first overlap length, edge length, number of reads, overlap offset, last overlap length)
		if(e->isSmallerEdge()) {
			if(e->getOrientation() == 0)
				out << "edge: { source:\"" << source << "\" target:\"" << destination
					<< "\" thickness: " << thickness << " arrowstyle: none backarrowstyle: solid color: red label: \"(" 
					<< e->getOverlapLen() << "," << e->getEdgeLength() << "," << num_reads << "," << e->getOverlapOffset()
					<< "," << e->getReverseEdge()->getOverlapLen()
					<< ")\" }" << "\n";
			else if(e->getOrientation() == 1)
				out << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
					<< thickness << " backarrowstyle:solid arrowstyle:solid color: green label: \"(" 
					<< e->getOverlapLen() << "," << e->getEdgeLength() << "," << num_reads << "," << e->getOverlapOffset()
					<< "," << e->getReverseEdge()->getOverlapLen()
					<< ")\" }" << "\n";
			else if(e->getOrientation() == 2)
				out << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
					<< thickness << " arrowstyle: none color: blue label: \"(" 
					<< e->getOverlapLen() << "," << e->getEdgeLength() << "," << num_reads << "," << e->getOverlapOffset()
					<< "," << e->getReverseEdge()->getOverlapLen()
					<< ")\" }" << "\n";
			else if(e->getOrientation() == 3)
				out << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
					<< thickness << " arrowstyle:solid color: red label: \"(" 
					<< e->getOverlapLen() << "," << e->getEdgeLength() << "," << num_reads << "," << e->getOverlapOffset()
					<< "," << e->getReverseEdge()->getOverlapLen()
					<< ")\" }" << "\n";
		}
	}
	out << "}";
	CLOCKSTOP;
	return out;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  readThreading
 *  Description:  Given the vector of reads that form the tiling for a long pacbio read,
 *  		  find mergings of edges in the graph that are supported by the tiling.
 *  		  These merges will then be done altogether next, and the edges to merge
 *  		  will be deleted from the graph after all the merges are done.
 * =====================================================================================
 */
void OverlapGraph::readThreading(const vector<UINT64> & readNumbers, vector< t_edge_vec > & list_edges_to_merge)
{
	if (readNumbers.size() <= 1){
		FILE_LOG(logDEBUG1) << "number of reads tiling PacBio read is " << readNumbers.size() 
			<< ", no threading" << "\n";
		return;
	}

	/* If there are more than two reads, do the following */
	
	// list of edges where the reads in the pacbio read reside
	t_edge_vec list_edges;
	// list of nodes included in the pacbio read
	vector<UINT64> list_nodes;
	// Keep track of indices of reads that are nodes
	vector<UINT64> index_of_nodes;

	vector<UINT64>::const_iterator it_p = readNumbers.cbegin();
	Edge *e1;
	UINT64 read_number, next_read_number;
	
	FILE_LOG(logDEBUG1) << "Loop through reads on pacbio read" << "\n";
	// Loop through all the reads in the pacbio read, look for edges to merge
	while(it_p != readNumbers.cend()){
		read_number = m_dataset->getReadIDMap()->at((*it_p) >> 1);
		/* If a read is a node in graph, instead of being contained in some edge
		 * Check the next read:
		 * 1. if next read is in an edge or not in the graph, 
		 *    no special treatment for the current node, continue
		 * 2. if next read is also a node, check and see if there is an edge between them
		 * 3. if this node is the last node in the tiling, done
		 * Then move on to check next read
		 */
		if (m_graph->at(read_number)->size() > 0){ 
			if((++it_p) != readNumbers.cend()){
				next_read_number = m_dataset->getReadIDMap()->at((*(it_p)) >> 1);
				if ( m_graph->at(next_read_number)->size()>0){
					e1 = findEdge(read_number, next_read_number);
					if(e1 != nullptr){
						FILE_LOG(logDEBUG1) << "two node reads: " << read_number 
							<< " " << next_read_number 
							<< " has an edge." << "\n";
						list_edges.push_back(e1);
					}
				}
			}
		}
		/* try to see if a read is located on an edge of the graph */
		else{
			e1 = locateRead(*it_p);
			/* read is not in the graph at all */
			if (e1 == nullptr){
				FILE_LOG(logDEBUG1) << "pacbio's read " << read_number << " is not in the graph" << "\n";
			}
			/* read is not a node but contained in some edge */
			else{                   
				FILE_LOG(logDEBUG1) << "pacbio's read " << read_number 
					<< " is on the edge" << *e1 << "\n";
				auto it_g = find(e1->getListOfReads()->cbegin(), e1->getListOfReads()->cend(), m_dataset->at(read_number));
				while(it_g != e1->getListOfReads()->cend() && it_p != readNumbers.cend()){
					read_number = m_dataset->getReadIDMap()->at(*it_p >> 1);
					if((*it_g)->getReadID() == read_number){
						++it_g;
						++it_p;
					}
					else{
						break;
					}
				}
				list_edges.push_back(e1);
				/* pacbio reaches end of an edge */
				if(it_g == e1->getListOfReads()->cend()){ 
					FILE_LOG(logDEBUG1) << "Pacbio read reached edge's end: " << *e1 << "\n";
				}
				/* the end of the pacbio read tiling is reached */
				if(it_p == readNumbers.cend()){
					FILE_LOG(logDEBUG1) << "Pacbio reached end, all agree with graph" << "\n";
					break;
				}
				else if(it_g != e1->getListOfReads()->cend()){
					FILE_LOG(logDEBUG1) << "Pacbio read does not agree with overlap graph at position: " << "\n"
						<< distance(readNumbers.cbegin(), it_p)  << " "
						<< read_number << " and " << "\n"
						<< *it_g << " edge: " << *e1
						<< "\n";
				}
			}
			++it_p;
		}
	}
	FILE_LOG(logDEBUG1) << "All Illumina reads traced!" << "\n";

	// split the list of edges into sets of edges that can be merged together
	if(list_edges.size() > 1){
		/* for debug */
		if(loglevel > 2){
			for(UINT64 i = 0; i < list_edges.size(); i++)
				FILE_LOG(logDEBUG1) << *(list_edges.at(i)) << "\n";
		}
		t_edge_vec mergeable_edges;
		mergeable_edges.push_back(list_edges.at(0));
		list_edges.erase(list_edges.begin());
		while(list_edges.size() > 0){
			Edge * edge1 = mergeable_edges.back();
			Edge * edge2 = list_edges.front();
			if( is_mergeable(edge1, edge2)){
				mergeable_edges.push_back(edge2);
			}
			else{
				// Make new empty vector of edges, swap content with mergeable_edges
				t_edge_vec m_edges;
				m_edges.swap(mergeable_edges);
				mergeable_edges.push_back(edge2);
				if(m_edges.size() > 1){
					FILE_LOG(logDEBUG1) << "edges can be merged:\n";
					for(UINT64 i = 0; i < m_edges.size(); i++)
						FILE_LOG(logDEBUG1) << *(m_edges.at(i)) << "\n";
					FILE_LOG(logDEBUG1) << "\n";
					list_edges_to_merge.push_back(m_edges);
				}
			}
			list_edges.erase(list_edges.begin());
			if(list_edges.size() == 0){
				t_edge_vec m_edges;
				m_edges.swap(mergeable_edges);
				mergeable_edges.push_back(edge2);
				if(m_edges.size() > 1){
					FILE_LOG(logDEBUG1) << "edges can be merged:\n";
					for(UINT64 i = 0; i < m_edges.size(); i++)
						FILE_LOG(logDEBUG1) << *(m_edges.at(i)) << "\n";
					FILE_LOG(logDEBUG1) << "\n";
					list_edges_to_merge.push_back(m_edges);
				}

			}
		}
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  printContigs
 *  Description:  Print contigs for all the edges in the graph, by streaming all the reads files.
 * =====================================================================================
 */
void OverlapGraph::printContigs(ostream & out, const vector<std::string> &readFilenameList)
{
	CLOCKSTART;
	UINT64 readID = 1;
	for(auto it = readFilenameList.cbegin(); it != readFilenameList.cend(); ++it){
		loadStringFromReadsFile(*it, readID);
	}

	t_edge_vec contigEdges; 
	getEdges(contigEdges);

	UINT64 printed_contigs(0);
	for ( auto it = contigEdges.cbegin(); it != contigEdges.cend(); ++it )
	{
//		FILE_LOG(logDEBUG1) << *(*it) << "\n";
		string contigString = (*it)->getEdgeString();

		if(contigString.length() >= minContigLengthTobeReported){
			++printed_contigs;
			out << ">contig_" << setfill('0') << setw(10) << printed_contigs
				<< " Edge ("  << (*it)->getSourceRead()->getReadID() << ", " 
				<< (*it)->getDestinationRead()->getReadID() 
				<< ") String Length: " << contigString.length() << "\n";

			UINT32 start=0;
			do
			{
				out << contigString.substr(start, 100) << "\n";  // save 100 BP in each line.
				start+=100;
			} while (start < contigString.length());
		}
		// If an edge is already shorter than the threshold, the following ones are
		// even shorter, since they are sorted by length
		else
			break;
	}
	FILE_LOG(logINFO) << "Total number of contigs printed: " << printed_contigs << "\n";
	CLOCKSTOP;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  loadStringFromReadsFile
 *  Description:  Fill the bases for all the edges in the graph, by streaming the reads in the
 *  		  reads file. Whenever a read is read from the reads file, its nucleotide sequence
 *  		  is streamed in, and the edges where this read resides will have the corresponding
 *  		  bases populated(determined). In the meantime, the positions with mismatches will
 *  		  have a count of each base option.
 * =====================================================================================
 */
void OverlapGraph::loadStringFromReadsFile(const std::string &read_file, UINT64 & readID)
{
	CLOCKSTART;
	FILE_LOG(logINFO) << "load read strings and fill the bases for the edges: " << read_file << "\n";

	// To count of reads in this file
	UINT64 readCount = 0;

	// Open file
	ifstream filePointer;
	filePointer.open(read_file.c_str());
	if(!filePointer.is_open()){
		FILE_LOG(logERROR) << "Unable to open file: " << read_file << "\n";
		return;
	}

	// Variables
	vector<std::string> line;
	std::string line0,line1, text;
	enum FileType {FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;

	while(!filePointer.eof())
	{
		// Check FASTA and FASTQ
		if(fileType == UNDEFINED) {
			getline (filePointer,text);
			if (text.length() > 0){
				if(text[0] == '>')
					fileType = FASTA;
				else if(text[0] == '@')
					fileType = FASTQ;
				else{
					FILE_LOG(logERROR) << "Unknown input file format."<<"\n";
					break;
				}
				filePointer.seekg(0, ios::beg);
			}
		}

		line.clear();

		// FASTA file read
		if(fileType == FASTA) {
			getline (filePointer,line0);	// get ID line
			getline (filePointer,line1,'>');	// get string line

			line1.erase(std::remove(line1.begin(), line1.end(), '\n'), 
					line1.end());
		}
		// FASTQ file read
		else if(fileType == FASTQ) {
			getline(filePointer, line0);	// ID
			getline(filePointer, line1);	// String
			// Ignore the next two lines
			filePointer.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			filePointer.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		// Get ReadID after removing the > or @ identifier and convert string to UINT64
		std::string readName="";
		if(line0[0] == '>' || line0[0] == '@')
			readName = line0.substr(1);
		else 
			readName = line0;

		populate_read(readID, line1);
		++readID;
		++readCount;
		if(readID % 1000000 == 0 && readID > 0){
			FILE_LOG(logINFO) << setw(10) << (readID / 1000000) << ",000,000"  
				<< " reads streamed, "
				<< setw(7) << checkMemoryUsage() << " MB\n";
		}
	}

	filePointer.close();
	FILE_LOG(logINFO) << setw(10) << readCount << " reads streamed from this read file\n";
	CLOCKSTOP;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  populate_read
 *  Description:  Given a read in the dataset, and its sequence, populate the strings for
 *   		  the edges.
 * =====================================================================================
 */
void OverlapGraph::populate_read(const UINT64 &readID, const std::string & read_str)
{
	Read *read = m_dataset->at(readID);
	std::string read_str_reverse = reverseComplement(read_str);

	// Edges with read as source or destination
	if(!(m_graph->at(readID)->empty())){
		for (auto it = m_graph->at(readID)->begin(); it != m_graph->at(readID)->end(); ++it){
			if((*it)->isSmallerEdge()){
				Edge *edge = *it;
				if (((edge->getOrientation() >> 1) & 1))
					edge->loadReadString(read_str, -1);
				else
					edge->loadReadString(read_str_reverse, -1);
			}
			else{
				Edge *edge = (*it)->getReverseEdge();
				if ((edge->getOrientation() & 1))
					edge->loadReadString(read_str, -2);
				else
					edge->loadReadString(read_str_reverse, -2);
			}
		}
	}
	// Edges with read on it
	vector< t_edge_loc_pair > *fwd_edges = read->getFwdEdges();
	vector< t_edge_loc_pair > *bwd_edges = read->getBwdEdges();

	for(auto it = fwd_edges->cbegin(); it != fwd_edges->cend(); ++it){
		it->first->loadReadString(read_str, it->second);
	}

	for(auto it = bwd_edges->cbegin(); it != bwd_edges->cend(); ++it){
		it->first->loadReadString(read_str_reverse, it->second);
	}
}
