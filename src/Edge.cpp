//============================================================================
// Name        : Edge.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Edge cpp file
//============================================================================

#include "Config.h"
#include "Edge.h"
#include <numeric>
extern TLogLevel loglevel;                      /* verbosity level of logging */


Edge::Edge(void)
{
	m_orient		= 0;
	m_overlapOffset		= 0;
	m_flag			= 0;
	m_string		= std::string();
	m_mismatch_basecount	= vector< vector<UINT64> >();
	m_edgeID		= 0;
	m_contigEdgeID		= 0;
	m_coverageDepth		= 0;
	m_SD			= 0;
	m_source 		= new Read;
	m_destination		= new Read;

	m_mismatches		= new t_vpair;
	m_listOfReads		= new vector<Read*>;
	m_listOfOverlapOffsets	= new vector<UINT32>;
	m_listOfOrientations	= new vector<UINT8>;
}

// Constructor for a (composite) edge, given
Edge::Edge(Read *source, Read *destination, UINT8 orient, UINT32 overlapOffset, t_vpair *mismatches, 
		vector<Read *> *listOfReads , vector<UINT32> *listOfOverlapOffsets, 
		vector<UINT8> * listOfOrientations )
:m_source(source), m_destination(destination), m_orient(orient), m_overlapOffset(overlapOffset),
	m_flag(0), m_string(std::string()), m_mismatches(mismatches), m_listOfReads(listOfReads), 
	m_listOfOverlapOffsets(listOfOverlapOffsets), m_listOfOrientations(listOfOrientations), 
	m_reverseEdge(nullptr),m_edgeID(0), m_contigEdgeID(0),m_coverageDepth(0.0), m_SD(0), m_flow(0)
{
	if(m_listOfOverlapOffsets && !m_listOfOverlapOffsets->empty() &&
			m_listOfOverlapOffsets->at(0) > m_source->getReadLength()){
		FILE_LOG(logERROR) << "ERROR: overlap offset " << m_listOfOverlapOffsets->at(0)
			<< " is greater than source read length " << m_source->getReadLength() << "\n";
	}
	// If this edge is a simple edge, make its reverse edge
	// For composite edge, leave the reverse edge to the sum of reverse edges
	if(source == destination)
		m_flag |= (1 << 1);
}

Edge::Edge(Read *source, Read *destination, UINT8 orient, UINT32 overlapOffset, t_vpair *mismatches)
:m_source(source), m_destination(destination), m_orient(orient), m_overlapOffset(overlapOffset),
	m_flag(0), m_string(std::string()),m_mismatches(mismatches), m_listOfReads(nullptr), 
	m_listOfOverlapOffsets(nullptr), m_listOfOrientations(nullptr),
	m_reverseEdge(nullptr),m_edgeID(0), m_contigEdgeID(0),m_coverageDepth(0.0), m_SD(0), m_flow(0)
{
	if(m_overlapOffset >= m_source->getReadLength()){
		FILE_LOG(logERROR) << "ERROR: overlap offset " << m_overlapOffset 
			<< " is greater than source read length " << m_source->getReadLength() << "\n";
	}
	// If this edge is a simple edge, make its reverse edge
	// For composite edge, leave the reverse edge to the sum of reverse edges
	if(source == destination)
		m_flag |= (1 << 1);
}

void Edge::copyEdge(const Edge &edge)
{
	//CLOCKSTART;
//	FILE_LOG(logDEBUG1) << edge << "\n";
	//FILE_LOG(logDEBUG1) << "Copy number variables\n";
	m_source = edge.m_source;
	m_destination = edge.m_destination;
	m_orient = edge.m_orient;
	m_overlapOffset = edge.m_overlapOffset;
	m_flag = edge.m_flag;
	m_string = edge.m_string;
	m_mismatch_basecount = edge.m_mismatch_basecount;
	m_edgeID = edge.m_edgeID;
	m_contigEdgeID = edge.m_contigEdgeID;
	m_coverageDepth = edge.m_coverageDepth;
	m_SD = edge.m_SD;
	m_flow = edge.m_flow;
	//FILE_LOG(logDEBUG1) << "Copy mismatches\n";
	if(edge.m_mismatches){
		m_mismatches = new t_vpair;
		for(auto it = edge.m_mismatches->cbegin(); it != edge.m_mismatches->cend(); ++it)
			m_mismatches->push_back(*it);
	}
	else
		m_mismatches = nullptr;

	//FILE_LOG(logDEBUG1) << "Copy list of reads\n";
	if(edge.m_listOfReads){
		m_listOfReads = new vector< Read* >;
		for(auto it = edge.m_listOfReads->cbegin(); it != edge.m_listOfReads->cend(); ++it)
			m_listOfReads->push_back(*it);
	}
	else 
		m_listOfReads = nullptr;

	//FILE_LOG(logDEBUG1) << "Copy list of overlap offsets\n";
	if(edge.m_listOfOverlapOffsets){
		m_listOfOverlapOffsets = new vector< UINT32 >;
		for(auto it = edge.m_listOfOverlapOffsets->cbegin(); it != edge.m_listOfOverlapOffsets->cend(); ++it)
			m_listOfOverlapOffsets->push_back(*it);
	}
	else 
		m_listOfOverlapOffsets = nullptr;

	//FILE_LOG(logDEBUG1) << "Copy list of orientations\n";
	if(edge.m_listOfOrientations){
		m_listOfOrientations = new vector< UINT8 >;
		for(auto it = edge.m_listOfOrientations->cbegin(); it != edge.m_listOfOrientations->cend(); ++it)
			m_listOfOrientations->push_back(*it);
	}
	else 
		m_listOfOrientations = nullptr;
	//CLOCKSTOP;
}

Edge::Edge(const Edge &edge)
{
	copyEdge(edge);
	m_reverseEdge = new Edge;
	m_reverseEdge->copyEdge(*(edge.m_reverseEdge));
	m_reverseEdge->setReverseEdge(this);
}

Edge& Edge::operator=(const Edge &edge)
{
	if(this == &edge)
		return *this;
	clearEdge();
	//CLOCKSTART;
	copyEdge(edge);
	m_reverseEdge = new Edge;
	m_reverseEdge->copyEdge(*(edge.m_reverseEdge));
	m_reverseEdge->setReverseEdge(this);
	//CLOCKSTOP;
	return *this;
}

Edge::~Edge()
{
	if(m_mismatches != nullptr){
		delete m_mismatches;
		m_mismatches = nullptr;
	}
	if(m_listOfReads != nullptr){
		delete m_listOfReads;
		m_listOfReads = nullptr;
	}
	if(m_listOfOverlapOffsets != nullptr){
		delete m_listOfOverlapOffsets;
		m_listOfOverlapOffsets = nullptr;
	}
	if(m_listOfOrientations != nullptr){
		delete m_listOfOrientations;
		m_listOfOrientations = nullptr;
	}
}

void Edge::clearEdge()
{
	//CLOCKSTART;
//	FILE_LOG(logDEBUG1) << "clearing edge " << *this << "\n"\n"
	if(m_mismatches != nullptr){
		delete m_mismatches;
		m_mismatches = nullptr;
	}
	if(m_listOfReads != nullptr){
		delete m_listOfReads;
		m_listOfReads = nullptr;
	}
	if(m_listOfOverlapOffsets != nullptr){
		delete m_listOfOverlapOffsets;
		m_listOfOverlapOffsets = nullptr;
	}
	if(m_listOfOrientations != nullptr){
		delete m_listOfOrientations;
		m_listOfOrientations = nullptr;
	}
	// This function does not delete the reverse edge
	//CLOCKSTOP;
}


void Edge::make_reverseEdge()
{
	// If reverse edge is already made, do not make it again
	if(m_reverseEdge)
		return;
	// Get the lists of info contained in the composite edges.
	// For simple edges, these will simply return nullptr.
	t_vpair* reverse_mismatches = getReverseMismatches();
	vector<Read* > *reverseListOfReads = getReverseListOfReads();
	vector<UINT32> *reverseListOfOvlapOffsets = getReverseListOfOverlapOffsets();
	vector<UINT8> *reverseListOfOrient = getReverseListOfOrient();
	// Make the new reverse edge
	// And set its reverse edge to the current edge
	m_reverseEdge = new Edge(m_destination, m_source, get_twin_orient(m_orient), 
			m_destination->getReadLength() + m_overlapOffset - m_source->getReadLength(),
			reverse_mismatches, reverseListOfReads, reverseListOfOvlapOffsets, reverseListOfOrient);

//	FILE_LOG(logDEBUG1) << "Edge: " << *this << "\n";
//	FILE_LOG(logDEBUG1) << "Reverse edge: " << *m_reverseEdge << "\n";
	m_reverseEdge->setReverseEdge(this);
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getReverseMismatches
 *  Description:  Ex) calculate reverseMismatchPosition (5AC -> 8CA)
 *                0123456789
 *                ACGTTAGTTA
 *                    |x|||||
 *                    TCGTTACAAG
 *                    9876543210
 *                reverseMismatchPosition = m_overlapOffset(4) + Read2Length(10) - mismatchPosition(5) - 1
 * =====================================================================================
 */
t_vpair * Edge::getReverseMismatches() const
{
	if(!m_mismatches){
		return nullptr;
	}
	t_vpair *reverseMismatches = new t_vpair;
	UINT32 mismatchIndex, reverseMismatchIndex, mismatchPosition, reverseMismatchPosition;
	for(UINT64 i = 0; i < m_mismatches->size(); i++) {
		mismatchIndex 		= m_mismatches->at(i).first;
		mismatchPosition 	= m_mismatches->at(i).second;

		reverseMismatchIndex 	= m_listOfReads->size() - mismatchIndex;
		reverseMismatchPosition = m_overlapOffset + m_destination->getReadLength() - mismatchPosition - 1;
		reverseMismatches->push_back(make_pair(reverseMismatchIndex, reverseMismatchPosition));
	}
	return reverseMismatches;
}

vector<UINT32> * Edge::getReverseListOfOverlapOffsets() const
{
	if (!m_listOfOverlapOffsets){
		return nullptr;
	}
	vector<UINT32> *reverseList = new vector<UINT32>;
	// First overlap offset in the reverse edge
	reverseList->push_back(getLastOverlapOffset()+m_destination->getReadLength() - m_listOfReads->back()->getReadLength());
	for(int i = m_listOfReads->size() - 2 ; i >= 0; --i){
		reverseList->push_back(m_listOfOverlapOffsets->at(i+1) + m_listOfReads->at(i+1)->getReadLength() - m_listOfReads->at(i)->getReadLength());
	}
	return reverseList;
}

/* Get the overlap length for the first link in the edge */
UINT32 Edge::getOverlapLen() const
{
	if(!m_listOfReads || m_listOfReads->empty())
		return (m_source->getReadLength() - m_overlapOffset);
	else
		return (m_source->getReadLength() - m_listOfOverlapOffsets->at(0));
}


/* Get the overlap length for the last link in the edge */
UINT32 Edge::getLastOverlapOffset() const
{
	if(!m_listOfReads || m_listOfReads->empty())
		return m_overlapOffset;
	UINT32 overlap_sum = std::accumulate(m_listOfOverlapOffsets->begin(), m_listOfOverlapOffsets->end(), 0);
	return (m_overlapOffset - overlap_sum);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  breakForwardEdge
 *  Description:  Break an edge at a specified link, only the forward edge
 * =====================================================================================
 */
vector<Edge*> Edge::breakForwardEdge(const UINT32 &link) const
{
	// sub edges generated by breaking a link from an edge
	vector<Edge*> sub_edges;

	/* If the edge is simple, no new edge is generated, return empty vector.
	 * If the edge is composite, insert two new edges, and then delete the current edge */
	if ( !m_listOfReads && !m_listOfReads->empty()){

		UINT64 num_reads = m_listOfReads->size();	// number of contained reads

		assert(link <= num_reads);// Link index should not exceed number of contained reads

		/* the position is not the first link, 
		 * there will be an edge with source the same as the original source */
		if(link != 0){                   
			UINT32 ovloffset = m_listOfOverlapOffsets->at(0);
			// Make the edge before the link
			Read* e1_dest_read = m_listOfReads->at(link - 1);
			UINT8 e1_orient = (m_orient & 2) + m_listOfOrientations->at(link-1);
			t_vpair * e1_mismatches = new t_vpair;
			vector<Read*> * e1_list_reads = nullptr;
			vector<UINT32> * e1_list_ovloffsets = nullptr;
			vector<UINT8> * e1_list_orients = nullptr;
			/* breaking at link 1, e1 is a simple edge */
			if(link == 1){
				// e1 mismatches should have mismatchIndex=0, 
				// do not need to update the mismatch position. 
				// New mismatch position should be same.
				for(auto it = m_mismatches->cbegin(); it != m_mismatches->cend(); ++it){
					if(it->first == 0){
						e1_mismatches->push_back(*it);
					}
				}
			}
			/* breaking after link 1, e1 is a composite edge */
			else{
				e1_list_reads = new vector<Read*>;
				e1_list_ovloffsets = new vector<UINT32>;
				e1_list_orients = new vector<UINT8>;
				UINT64 i = 0;
				for(i = 0; i < (link-1); i++){
					e1_list_reads->push_back(m_listOfReads->at(i));
					e1_list_ovloffsets->push_back(m_listOfOverlapOffsets->at(i));
					e1_list_orients->push_back(m_listOfOrientations->at(i));
					ovloffset += m_listOfOverlapOffsets->at(i+1);
				}
				// e1 mismatches should have mismatchIndex<=(link-1), 
				for(auto it = m_mismatches->cbegin(); it != m_mismatches->cend(); ++it){
					if(it->first <= link - 1){
						e1_mismatches->push_back(*it);
					}
				}
			}
			Edge * e1 = new Edge(m_source, e1_dest_read, e1_orient, 
					ovloffset, e1_mismatches, e1_list_reads, e1_list_ovloffsets, 
					e1_list_orients);
			sub_edges.push_back(e1);
			FILE_LOG(logDEBUG1) << "First sub edge is " << *e1 << "\n";
		}
		/* the position is not the last link, 
		 * there will be an edge with destination the same as the original destination */
		if(link != num_reads){          
			UINT32 mismatchIndex, mismatchPosition;
			UINT32 ovloffset = getLastOverlapOffset(); /* last overlap offset */

			Read* e2_source_read = m_listOfReads->at(link);
			UINT8 e2_orient = (m_orient & 1) + (m_listOfOrientations->at(link)  << 1); 
			t_vpair * e2_mismatches = new t_vpair;
			vector<Read*> * e2_list_reads = nullptr;
			vector<UINT32> * e2_list_ovloffsets = nullptr;
			vector<UINT8> * e2_list_orients = nullptr;

			/* broken at the second last link,
			 * e2 is a simple edge */
			if(link == (num_reads-1)){
				// e2 mismatches should have mismatchIndex=link+1, 
				// update the index and mismatch position.
				for(auto it = m_mismatches->cbegin(); it != m_mismatches->cend(); ++it){
					if(it->first == link+1){
						mismatchIndex = 0;
						mismatchPosition = it->second - std::accumulate(m_listOfOverlapOffsets->begin(), m_listOfOverlapOffsets->end(), 0);
						e2_mismatches->push_back(make_pair(mismatchIndex, mismatchPosition));
					}
				}
			}
			/* e2 is a composite edge */
			else{
				e2_list_reads = new vector<Read*>;
				e2_list_ovloffsets = new vector<UINT32>;
				e2_list_orients = new vector<UINT8>;
				UINT64 i = 0;
				for(i = link+1; i < num_reads; i++){
					e2_list_reads->push_back(m_listOfReads->at(i));
					e2_list_ovloffsets->push_back(m_listOfOverlapOffsets->at(i));
					e2_list_orients->push_back(m_listOfOrientations->at(i));
					ovloffset += m_listOfOverlapOffsets->at(i);
				}
				// e2 mismatches should have mismatchIndex > link, 
				// update the index and mismatch position.
				for(auto it = m_mismatches->cbegin(); it != m_mismatches->cend(); ++it){
					if(it->first > link){
						mismatchIndex = it->first - (link+1);
						mismatchPosition = it->second - std::accumulate(m_listOfOverlapOffsets->begin(), m_listOfOverlapOffsets->begin() + link + 1, 0);
						e2_mismatches->push_back(make_pair(mismatchIndex, mismatchPosition));
					}
				}
			}
			Edge *e2 = new Edge(e2_source_read, m_destination, e2_orient, ovloffset, e2_mismatches, e2_list_reads, e2_list_ovloffsets, e2_list_orients);
			sub_edges.push_back(e2);
			FILE_LOG(logDEBUG1) << "Second sub edge is " << *e2 << "\n";
		}
	}
	return sub_edges;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  breakEdge
 *  Description:  break a given edge into (at most) 2 bi-directed edges.
 *  First create sub-edges after the specified link is taken out, both forward and reverse;
 *  then insert these sub-edges in the graph, then remove the original edge.
 * =====================================================================================
 */
vector<Edge*> Edge::breakEdge(const UINT32 &link) const
{
	vector<Edge*> sub_edges;
	if(m_listOfReads && !m_listOfReads->empty())
	{
		/* get new forward subedges */
		vector<Edge*> forward_sub_edges = this->breakForwardEdge(link);           

		/* get new backward subedges */
		vector<Edge*> backward_sub_edges = m_reverseEdge->breakForwardEdge(m_listOfReads->size() - link); 
		UINT64 num_edges = forward_sub_edges.size();
		for(UINT64 i = 0; i < num_edges; i++) {
			// set twin edges
			forward_sub_edges.at(i)->setReverseEdge(backward_sub_edges.at(num_edges-i-1));
			backward_sub_edges.at(num_edges-i-1)->setReverseEdge(forward_sub_edges.at(i));
			sub_edges.push_back(forward_sub_edges.at(i));
		}
	}
	return sub_edges;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getBaseByBaseCoverage
 *  Description:  Calculate the coverage depth of an edge for every basepair and then update
 *                the Mean and SD of coverage depth in the edge. Only consider reads that are
 *                unique to the edge.
 * =====================================================================================
 */
void Edge::getBaseByBaseCoverage()
{
	// Array length is the same as the std::string length in the edge.
	UINT64 length = getEdgeLength();
	vector<UINT64> coverageBaseByBase;
	coverageBaseByBase.reserve(length);

	for(UINT64 i = 0; i < length; i++) {
		coverageBaseByBase.push_back(0);	// At first all the bases are covered 0 times.
	}
	UINT64 overlapOffset(0);

	// Increment the coverage of the section that each read covers,
	// NOT counting the source read and destination read because they are shared among multiple edges
	// JJ: so? if they are shared, their coverage should just be discarded?
	//     If indeed discarded, then the coverage vector can be shorter
	if(m_listOfReads){
		for(UINT64 i = 0; i < m_listOfReads->size(); i++)	// For each read in the edge.
		{
			Read *read = m_listOfReads->at(i);
			// Where the current read starts in the std::string.
			overlapOffset += m_listOfOverlapOffsets->at(i);	
			if(read->getFwdEdges()->size() == 1){
				for(UINT64 j = overlapOffset; j < overlapOffset + read->getReadLength(); ++j) {
					coverageBaseByBase.at(j) += 1;
				}
			}
		}
	}
	m_coverageDepth = get_mean(coverageBaseByBase);
	m_SD = get_sd(coverageBaseByBase);
}

//=============================================================================
// Merge two edges in the overlap graph.
//=============================================================================
Edge* Add( const Edge *edge1, const Edge *edge2)
{
	assert (is_mergeable(edge1, edge2));
	Edge *merge_forward = merge_forward_edges(*edge1, *edge2);
	Edge *merge_reverse = merge_forward_edges(*(edge2->m_reverseEdge), *(edge1->m_reverseEdge));
	merge_forward->setReverseEdge(merge_reverse);
	merge_reverse->setReverseEdge(merge_forward);
	return merge_forward;
}

/* Edge operator+( const Edge & edge1, const Edge & edge2)
 * {
 * 	assert (is_mergeable(&edge1, &edge2));
 * 	CLOCKSTART;
 * 	Edge merge_forward = merge_forward_edges(edge1, edge2);
 * 	Edge merge_reverse = merge_forward_edges(*(edge2.m_reverseEdge), *(edge1.m_reverseEdge));
 * 	merge_forward.setReverseEdge(&merge_reverse);
 * 	merge_reverse.setReverseEdge(&merge_forward);
 * 	CLOCKSTOP;
 * 	return merge_forward;
 * }
 */

void Edge::setReverseEdge(Edge * edge)
{ 
	if(m_reverseEdge == edge)
		return;
	else if(!m_reverseEdge) {
		m_reverseEdge = edge;
	}
	else{
		m_reverseEdge = edge;
	}
}

// Assisting function for just merging forward edges
Edge*  merge_forward_edges(const Edge & edge1, const Edge & edge2)
{

	// Orientation
	UINT8 orientationForward = mergedEdgeOrientation(edge1.m_orient, edge2.m_orient);

	// overlap offset
	UINT32 overlapOffsetForward = edge1.m_overlapOffset + edge2.m_overlapOffset;

	t_vpair * mismatchesForward = new t_vpair;
	vector<Read* > * listReadsForward = new vector<Read* >;
	vector<UINT32> * listOverlapOffsetsForward= new vector<UINT32>;
	vector<UINT8> * listOrientationsForward = new vector<UINT8>;

	// Merge the lists from the two edges.
	mergeList(&edge1, &edge2, listReadsForward, listOverlapOffsetsForward, listOrientationsForward, mismatchesForward); 

	// Make the forward edge
	Edge *edgeForward = new Edge(edge1.m_source,edge2.m_destination,orientationForward, 
			overlapOffsetForward, mismatchesForward, listReadsForward, 
			listOverlapOffsetsForward, listOrientationsForward);

	return edgeForward;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  mergeEdges
 *  Description:  merge a list of edges, not restricted to only two edges
 * =====================================================================================
 */
Edge * mergeEdges(const vector<Edge *> & list_edges)
{
	if (list_edges.size() < 2){
		FILE_LOG(logERROR) << "There are fewer than 2 edges in the list of edges to merge." << "\n";
		return nullptr;
	}

	if(loglevel > 2){
		FILE_LOG(logDEBUG1) << "Merging following list of edges: " << "\n";
		for(auto it = list_edges.cbegin(); it != list_edges.cend(); ++it){
			FILE_LOG(logDEBUG1) << *(*it) << "\n";
		}
		FILE_LOG(logDEBUG1) << "\n";
	}


	FILE_LOG(logDEBUG1) << "Add first two edges \n";
	Edge *sum_edge = Add(list_edges.at(0), list_edges.at(1));
	FILE_LOG(logDEBUG1) << *sum_edge << "\n";
	for(size_t i = 2; i < list_edges.size(); ++i){
		Edge *temp_edge = Add(sum_edge, list_edges.at(i));
		FILE_LOG(logDEBUG1) << *temp_edge << "\n";
		delete sum_edge;
		sum_edge = temp_edge;
		FILE_LOG(logDEBUG1) << *sum_edge << "\n\n";
	}
	FILE_LOG(logDEBUG1) << "\nFinal edge resulted is: \n" << *sum_edge << "\n";
	return sum_edge;
}

//=============================================================================
// Merge the list of reads, list of overlap offsets and list of orientations of two edges.
//=============================================================================
void mergeList(const Edge *edge1, const Edge *edge2, 
		vector<Read* > *listReads, vector<UINT32> *listOverlapOffsets, vector<UINT8> *listOrientations, 
		t_vpair *mismatchesForward)
{
//	CLOCKSTART;
	// Copy the list from edge1.
	if (edge1->getListOfReads() && !edge1->getListOfReads()->empty() ){
		for(auto it = edge1->getListOfReads()->cbegin(); it != edge1->getListOfReads()->cend(); ++it){
			listReads->push_back(*it);
		}
		for(auto it = edge1->getListOfOverlapOffsets()->cbegin(); it != edge1->getListOfOverlapOffsets()->cend(); ++it){
			listOverlapOffsets->push_back(*it);
		}
		for(auto it = edge1->getListOfOrientations()->cbegin(); it != edge1->getListOfOrientations()->cend(); ++it){
			listOrientations->push_back(*it);
		}
	}

	if (edge1->getMismatches() && !edge1->getMismatches()->empty())
		*mismatchesForward = *(edge1->getMismatches());

	// Insert the common node of the two edges
	listReads->push_back(edge1->getDestinationRead()); 	
	// last overlap offset from the first edge
	listOverlapOffsets->push_back(edge1->getLastOverlapOffset());	
	// Orientation of the common node. 
	listOrientations->push_back((edge1->getOrientation() & 1));

	// Concatenate the list from edge2.
	if (edge2->getListOfReads() && !edge2->getListOfReads()->empty()){
		for(auto it = edge2->getListOfReads()->cbegin(); it != edge2->getListOfReads()->cend(); ++it){
			listReads->push_back(*it);
		}
		for(auto it = edge2->getListOfOverlapOffsets()->cbegin(); it != edge2->getListOfOverlapOffsets()->cend(); ++it){
			listOverlapOffsets->push_back(*it);
		}
		for(auto it = edge2->getListOfOrientations()->cbegin(); it != edge2->getListOfOrientations()->cend(); ++it){
			listOrientations->push_back(*it);
		}
	}

	// index of reads with mismatches need to be updated (increased by number of reads in edge1 + 1)
	if (edge2->getMismatches() && !edge2->getMismatches()->empty()){
		for (UINT32 i = 0; i < edge2->getMismatches()->size(); i++) {
			UINT32 mismatchIndex = edge1->getListOfReads()->size() + 1 + edge2->getMismatches()->at(i).first;
			UINT32 mismatchPosition = edge1->getOverlapOffset() + edge2->getMismatches()->at(i).second;
			mismatchesForward->push_back(make_pair(mismatchIndex, mismatchPosition));
		}
	}

//	CLOCKSTOP;
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  operator<<
 *  Description:  overloading of operator <<, for printing edge information
 * =====================================================================================
 */
ostream& operator<< (ostream &out, const Edge &edge)
{
	out << "Edge ID: " << setfill(' ') << setw(6) << edge.m_edgeID
		<< " Orient: "<< setfill(' ') << setw(3) << (int)edge.m_orient 
//		<< " Flow: " << setfill(' ') << setw(3) << edge.m_flow 
		<< " Reads: " << setfill(' ') << setw(7) << (edge.m_listOfReads ? (edge.m_listOfReads)->size() : 0)
		<< " Length: " << setfill(' ') << setw(10) << edge.getEdgeLength() << "\t";
	out << setw(10) << setfill(' ') << (edge.m_source)->getReadID() << " (--" << setw(7) << setfill(' ') << edge.m_overlapOffset << ", " 
		<< setfill(' ') << setw(3)<< (int)edge.m_orient << "--) " << setw(10) << setfill(' ') << (edge.m_destination)->getReadID();
//	if(edge.m_listOfReads && edge.m_listOfReads->size() > 0){
//		out << "\t";
//		for(auto it = edge.m_listOfReads->begin(); it != edge.m_listOfReads->end(); ++it)
//			out << setfill(' ') << setw(10) << (*it)->getReadID();
//	}
	return out;
}


//=============================================================================
// Check if two edges can be merged into one edge
// For two edges e1(u,v) and e2(v,w), at node v, 
// one of the edges should be an incoming edge and the other should be an outgoing
// edge to match.
//=============================================================================
bool is_mergeable(const Edge *edge1, const Edge *edge2)
{
	// First, the destination of edge1 has to be the same as 
	// the source read of edge2
	if (edge1->getDestinationRead()->getReadID() != edge2->getSourceRead()->getReadID()){
		return false;
	}

	// *-----> and >------* *1 and 1*
	// *-----< and <------* *0 and 0*
	else if ((edge1->getOrientation() & 1) == ((edge2->getOrientation() >>1) & 1))
		return true;
	else{
		return false;
	}
}

// Orientation of the edge when two edges are merged.
UINT8 mergedEdgeOrientation(const UINT8 &orient1, const UINT8 &orient2)
{
//	assert(orient1 < 4 && orient2 < 4);
	return ((orient1 & 2) | (orient2 & 1));
}

UINT8 get_twin_orient(const UINT8 &orient)
{
	assert(static_cast<int>(orient) < 4); // Quit if orient is not 0, 1, 2, or 3
	// Exchange the last two bits of the orient, and flip them
	UINT8 twin_orient = ((orient >> 1) ^ 1) | (((orient & 1) ^ 1) << 1) ;
	assert(static_cast<int>(twin_orient) < 4);
	return twin_orient;
}

vector<UINT8>* Edge::getReverseListOfOrient()
{
	if(!m_listOfOrientations){
		return nullptr;
	}
	// Same number of elements, but orientations are flipped
	vector<UINT8> *reverseList = new vector<UINT8>;
	reverseList->reserve(m_listOfOrientations->size());
	for( auto rit = m_listOfOrientations->crend(); rit != m_listOfOrientations->crbegin(); ++rit)
		reverseList->push_back((*rit) ^ 1);
	return reverseList;
}


vector<Read *> * Edge::getReverseListOfReads()
{
	if(!m_listOfReads){
		return nullptr;
	}
	vector<Read* > * reverseList = new vector<Read* >;
	reverseList->reserve(m_listOfReads->size());
	for(auto it = m_listOfReads->rbegin(); it != m_listOfReads->rend(); ++it){
		reverseList->push_back(*it);
	}
	return reverseList;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  updateEdgeInfo
 *  Description:  Update a read's residing edge information according to the 
 *  		  insertion of new edge of deletion of existing edge
 *  		  TODO: in the case of loop reduction, it's possible that one read 
 *  		  appears more than once in an edge. Should this function be called more
 *  		  than once in such cases?
 * =====================================================================================
 */
void Read::updateEdgeInfo(Edge *edge, UINT32 read_index, EdgeOP operation)
{
	//CLOCKSTART;
	// Insert an edge with this read included on it
	if (operation == INSERTION){
		if((edge->getListOfOrientations()->at(read_index) & 1) == 1){
			m_fwdEdges->push_back(std::make_pair(edge, read_index));
		}
		else{
			m_bwdEdges->push_back(std::make_pair(edge, read_index));
		}
	}
	// Delete an edge with this read included on it
	else{
		t_edge_loc_pair pair_to_rm = std::make_pair(edge, read_index);
		if((edge->getListOfOrientations()->at(read_index) & 1) == 1){
			m_fwdEdges->erase(
			std::remove(m_fwdEdges->begin(),m_fwdEdges->end(), pair_to_rm),
					m_fwdEdges->end());
		}
		else{
			m_bwdEdges->erase(
			std::remove(m_bwdEdges->begin(),m_bwdEdges->end(), pair_to_rm),
					m_bwdEdges->end());
		}
	}
	//CLOCKSTOP;
}


void Edge::updateReadsLocations(EdgeOP operation) 
{
	if(m_listOfReads){
		for(UINT64 i = 0; i < m_listOfReads->size(); ++i){
			m_listOfReads->at(i)->updateEdgeInfo(this, i, operation);
		}
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  loadReadString
 *  Description:  Fill the string for the current edge from a read and its string
 * =====================================================================================
 */
void Edge::loadReadString(const std::string read_str, int index)
{
	// If this edge is not smaller edge ( source smaller than the destination ID)
	// or its length is smaller than the reporting threshold,
	// do not load string in this edge.
	if (!isSmallerEdge() || getEdgeLength() < minContigLengthTobeReported){
		return;
	}
	// TODO: need to update the mismatch base count
	assert(index >= -2);
	if(m_listOfReads)
		assert(index < static_cast<int>(m_listOfReads->size()) );
	if(m_string.length()==0)
		m_string = std::string(getEdgeLength(), 'N');
	// Read is the source of the edge, get the first overlap offset length
	if(index == -1){
		UINT32 len = getFirstOverlapOffset();
		m_string.replace(0, len, read_str, 0, len);
	}
	// Read is the destination of the edge
	else if (index == -2){
		m_string.replace(m_overlapOffset, read_str.length(), read_str,
				0, read_str.length());
	}
	// Read is neither source nor destination, but lives on the edge
	else{
		UINT32 start;
		UINT32 len;
		if(index < static_cast<int>(m_listOfOverlapOffsets->size() - 1))  // Not the last read on edge
		{
			len = m_listOfOverlapOffsets->at(index+1);
			start = std::accumulate(m_listOfOverlapOffsets->begin(), m_listOfOverlapOffsets->begin()+index+1, 0);
		}
		else {
			len = getLastOverlapOffset();
			start = std::accumulate(m_listOfOverlapOffsets->begin(), m_listOfOverlapOffsets->end(), 0);
		}

		m_string.replace(start, len, read_str, 0, len);
	}
	if(m_string.length() != getEdgeLength()){
		FILE_LOG(logERROR) << "ERROR: edge length changed!!! " << m_string.length() << " and " << getEdgeLength() << "\n";
	}
}


template <typename T>
float get_mean(const vector<T> &numbers)
{
	if (numbers.empty())
		return 0;
	float sum(0);
	sum = std::accumulate(numbers.begin(), numbers.end(), sum);
	return (sum/static_cast<float>(numbers.size()));
}

template <typename T>
float get_sd(const vector<T> &numbers)
{
	if(numbers.empty())
		return 0;
	float mean = get_mean(numbers);
	float variance(0.0f);
	for(auto it = numbers.begin(); it != numbers.end(); ++it){
		variance += (*it - mean) * (*it - mean);
	}
	return (sqrt(variance/numbers.size()));
}
