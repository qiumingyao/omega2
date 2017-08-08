/*
 * =====================================================================================
 *
 *       Filename:  tOverlapGraph.h
 *
 *    Description:  template class tOverlapGraph
 *
 *        Version:  1.0
 *        Created:  07/24/2015 11:55:29
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  JJ Crosskey (cjg), jjchai01@gmail.com
 *   Organization:  ORNL
 *
 * =====================================================================================
 */

#ifndef TOVERLAPGRAPH_H_
#define TOVERLAPGRAPH_H_

#include "Config.h"
#include "tEdge.h"
#include "tDataSet.h"
#include "tRead.h"
extern TLogLevel loglevel;/*  verbosity level of logging */
extern string outputFilenamePrefix;

/*
 * =====================================================================================
 *        Class:  tOverlapGraph
 *  Description:  overlap graph class with template, where read is vector of element of
 *                type T
 * =====================================================================================
 */

enum nodeType {
	NODE_UNEXPLORED = 0, // Current node u is not explored yet. Meaning that there is no edge (u,v) in the graph.
	NODE_EXPLORED = 1, //  Current node u is explored. Meaning that all edges (u,v) are inserted in the dataset.
	NODE_EXPLORED_AND_TRANSITIVE_EDGES_MARKED = 2 // Meaning that all transitive edges (u,v) of current node u is marked and its neighbors transitive edges are also marked. Now it is safe to remove any transitive edge (u,v) from node u.
};

enum markType{
	MARK_VACANT = 0,
	MARK_INPLAY = 1,
	MARK_ELIMINATED = 2
};

template <typename T>
class tOverlapGraph
{
	public:
		/* ====================  LIFECYCLE     ======================================= */

		tOverlapGraph ();                             /* constructor      */

		~tOverlapGraph ();                            /* destructor       */

		tOverlapGraph( tDataSet<tRead<T>> * data_set, UINT64 min_ovl=1);

		/* ====================  ACCESSORS     ======================================= */

		UINT64 getNumberOfNodes(void) const {return numberOfNodes;}

		UINT64 getNumberOfEdges(void) const {return numberOfEdges;}

		void print_neighbors(UINT64 readNumber) const;

		bool getEdges(vector<tEdge<tRead<T>> *> & contigEdges);

		void printGraph(string graphFileName, const vector<tEdge<tRead<T>> *> & contigEdges);

		vector<vector<T>>*  getContigs(const vector<tEdge<tRead<T>> *> & contigEdges, 
				bool isolates = true);
		/* ====================  MUTATORS      ======================================= */

//		void setDataSet(tDataSet<T> * dataset){dataSet = dataset;}

		/* ====================  OPERATORS     ======================================= */
		// TODO: depth first search to find all paths
		// void DFS(void) const;
		//

	private:
		/* ====================  METHODS       ======================================= */
		void buildHashTable(void);
		
		void buildOverlapGraph(void);

		void insertEdge( tEdge<tRead<T>> * edge);

		void insertEdge(const tRead<T> &read1, const tRead<T> &read2, char orient, UINT64 overlapOffset);

		void removeEdge( tEdge<tRead<T>> * edge);

		UINT64 contractCompositeEdges(void);

		vector<tEdge<tRead<T>> *> findEdges(UINT64 source, UINT64 destination);
		
		void markContainedReads(void);

		bool checkOverlapForContainedRead(const tRead<T> & r1, const tRead<T> & r2, char orient, UINT64 start);

		void insertAllEdgesOfRead ( UINT64 readNumber, vector<nodeType> * exploredReads );

		bool checkOverlap ( const tRead<T> &read1, const tRead<T> &read2, UINT64 orient, UINT64 start );

		void markTransitiveEdges(UINT64 readNumber, vector<markType> * markedNodes);

		void removeTransitiveEdges ( UINT64 readNumber );

		vector<T> getStringInEdge (const tEdge<tRead<T>> *edge);
		/* ====================  DATA MEMBERS  ======================================= */
		UINT64 minOvl;

		tDataSet<tRead<T>> *dataSet;

		map<vector<T>, vector<UINT64> > *hashTable;

		vector< vector<tEdge<tRead<T>>* > *> * graph;

		vector< vector<T> > * contigStrings;

		UINT64 numberOfNodes;

		UINT64 numberOfEdges;


}; /* -----  end of class tOverlapGraph  ----- */

	template < typename T >
bool compareEdges(const T * const & edge1, const T * const & edge2)
{
	return (*edge1 < *edge2);
}

	template < typename T >
tOverlapGraph< T >::tOverlapGraph ()
{
	minOvl = 0;
	dataSet = nullptr;
	hashTable = nullptr;
	graph = nullptr;
	contigStrings = nullptr;
	numberOfNodes = 0;
	numberOfEdges = 0;
	cout << "finished initialization of overlapgraph without arguments" << endl;
}

	template < typename T >
tOverlapGraph< T >::~tOverlapGraph ()
{
	cout << "deconstructing overlap graph object" << endl;
	if (graph != nullptr){
		for(UINT64 i = 0; i < graph->size(); i++)
		{
			if (graph->at(i) != nullptr)
				delete graph->at(i);
		}
	}
	cout << "graph successfully deleted" << endl;
}


	template < typename T >
tOverlapGraph < T >::tOverlapGraph (tDataSet<tRead<T>> * data_set, UINT64 min_ovl)
{
	minOvl = min_ovl;

	dataSet = data_set;

	if(dataSet->isDupRemoved())
		dataSet->rmDupReads();

	numberOfNodes = 0;

	numberOfEdges = 0;

	buildHashTable();

	buildOverlapGraph();

	contractCompositeEdges();

}


	template < typename T >
void tOverlapGraph<T>::buildHashTable ( void )
{
	CLOCKSTART;
	// TODO: for now the hash table is just a map, instead of a unordered_map.
	// Since the hash function has not been designed for the templated class.
	hashTable = new map<vector<T>, vector<UINT64> >;
	for(UINT64 i = 1; i <= dataSet->numberOfUniqueReads; i++){
		tRead<T> read = dataSet->getReadFromID(i);
		UINT64 read_length = read.size();

		// Forward prefix: 00 = 0	0000 --> 0x0
		// Forward suffix: 01 = 1	0100 --> 0x4
		// Reverse prefix: 10 = 2	1000 --> 0x8
		// Reverse suffix: 11 = 3	1100 --> 0xC
		vector<T> forward_prefix = read.substr(0, minOvl, 0);
		if(hashTable->count(forward_prefix) == 0 ){
			vector<UINT64> readIDs;
			(*hashTable)[forward_prefix] = readIDs;
		}
		(hashTable->at(forward_prefix)).push_back(i + 0x0000000000000000);

		vector<T> forward_suffix = read.substr(read_length - minOvl, minOvl, 0);
		if(hashTable->count(forward_suffix) == 0 ){
			vector<UINT64> readIDs;
			(*hashTable)[forward_suffix] = readIDs;
		}
		(hashTable->at(forward_suffix)).push_back(i + 0x4000000000000000);

		vector<T> backward_prefix = read.substr(0, minOvl, 1);
		if(hashTable->count(backward_prefix) == 0 ){
			vector<UINT64> readIDs;
			(*hashTable)[backward_prefix] = readIDs;
		}
		(hashTable->at(backward_prefix)).push_back(i + 0x8000000000000000);

		vector<T> backward_suffix = read.substr(read_length - minOvl, minOvl, 1);
		if(hashTable->count(backward_suffix) == 0 ){
			vector<UINT64> readIDs;
			(*hashTable)[backward_suffix] = readIDs;
		}
		(hashTable->at(backward_suffix)).push_back(i + 0xC000000000000000);
	}
	cout << "Size of hashTable is " << hashTable->size() << endl;
	CLOCKSTOP;
}


	template < typename T >
void tOverlapGraph<T>::markContainedReads ( void )
{
	CLOCKSTART;
	UINT64 counter = 0;
	for(UINT64 i = 1; i <= dataSet->numberOfUniqueReads; ++i){
		tRead<T> &read_i = dataSet->getReadFromID(i);
		if(read_i.superID == 0){
			for(UINT64 index=0; index < (read_i.size() - minOvl); ++index){
				vector<T> sub_read = read_i.substr(index, minOvl);
				if(hashTable->count(sub_read) > 0){
					for(UINT64 k = 0; k < hashTable->at(sub_read).size(); k++){
						UINT64 hash_val = hashTable->at(sub_read).at(k);
						UINT64 readID = (hash_val & 0x0FFFFFFFFFFFFFFF);
						tRead<T> &read_k = dataSet->getReadFromID(readID);
						if(read_k.superID == 0 && i != readID){
							char orient = (hash_val >> 62);
							if(checkOverlapForContainedRead(read_i, read_k, orient, index)){
								read_k.superID = i;
								FILE_LOG(logDEBUG1) << "Set read " << readID << "'s superID to " << i << endl;
								++counter;
							}
						}
					}
				}
			}
		}
	}
	FILE_LOG(logINFO) << "Number of contained reads: " << counter << endl;
	CLOCKSTOP;
}		/* -----  end of method TOverlapGraph<T>::markContainedReads  ----- */

	template < typename T >
bool tOverlapGraph<T>::checkOverlapForContainedRead ( const tRead<T> & read1, const tRead<T> & read2, 
		char orient, UINT64 start )
{
	FILE_LOG(logDEBUG1) << "Check if " << read1 << " contains " << read2 << endl;
	if(orient == 0 || orient == 2)
		// orient 0 (00): forward prefix
		//   >--------MMMMMMMMMMMMMMM*******------> read1      M means match found by hash table
		//            MMMMMMMMMMMMMMM*******>       read2      * means we need to check these characters for match
		//				OR
		// orient 2 (10): reverse prefix
		//	 >---*****MMMMMMMMMMMMMMM*******------> read1
		//		  MMMMMMMMMMMMMMM*******<	Reverese complement of read2
	{
		UINT64 lengthRemaining1 = read1.size() - start - minOvl; 	// This is the remaining of read1
		UINT64 lengthRemaining2 = read2.size() - minOvl; // This is the remaining of read2
		if (lengthRemaining2 == 0)
			return true;
		else if(lengthRemaining1 >= lengthRemaining2)
		{
			return read1.substr(start + minOvl, lengthRemaining2) == read2.substr(minOvl, lengthRemaining2, (orient >> 1)); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	else							
		// orient 1 (01): forward suffix
		//   >---*****MMMMMMMMMMMMMMM-------------> 	read1      M means match found by hash table
		//      >*****MMMMMMMMMMMMMMM       		read2      * means we need to check these characters for match
		// OR
		// orient 3 (11): reverse suffix
		//	 >---*****MMMMMMMMMMMMMMM-------------> read1
		//	    <*****MMMMMMMMMMMMMMM		Reverse Complement of Read2
	{
		UINT64 lengthRemaining1 = start;
		UINT64 lengthRemaining2 = read2.size() - minOvl;
		if (lengthRemaining2 == 0)
			return true;
		else if(lengthRemaining1 >= lengthRemaining2)
		{
			return read1.substr(start - lengthRemaining2, lengthRemaining2) == read2.substr(0, lengthRemaining2, (orient >> 1)); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	FILE_LOG(logDEBUG1) << read1 << " does not contain " << read2 << endl;
	return false;
}

	template < typename T >
void tOverlapGraph<T>::buildOverlapGraph ( void )
{
	CLOCKSTART;

	UINT64 counter 				= 0;
	vector<nodeType> *exploredReads 	= new vector<nodeType>;
	vector<UINT64> * queue 			= new vector<UINT64>;
	vector<markType> *markedNodes 		= new vector<markType>;
	graph 					= new vector< vector<tEdge<tRead<T>>* > * >;
	exploredReads->reserve(dataSet->numberOfUniqueReads+1);
	queue->reserve(dataSet->numberOfUniqueReads+1);
	markedNodes->reserve(dataSet->numberOfUniqueReads+1);
	graph->reserve(dataSet->numberOfUniqueReads+1);

	for(UINT64 i = 0; i <= dataSet->numberOfUniqueReads; i++) // Initialization
	{
		vector<tEdge<tRead<T>> * > *newList = new vector<tEdge<tRead<T>>* >;
		graph->push_back(newList);
		exploredReads->push_back(NODE_UNEXPLORED);
		queue->push_back(0);
		markedNodes->push_back(MARK_VACANT);
	}

	markContainedReads();

	FILE_LOG(logDEBUG1) << *dataSet << endl;
	for(UINT64 i = 1; i <= dataSet->numberOfUniqueReads; i++)
	{
		if(exploredReads->at(i) == NODE_UNEXPLORED)
		{
			UINT64 start(0), end(0); 	// Initialize queue start and end.
			queue->at(end++) = i;
			while(start < end) 	// This loop will explore all connected component starting from read i.
			{
				counter++;
				UINT64 read1 = queue->at(start++);
				FILE_LOG(logDEBUG1) << "read " << read1 << endl;
				if(exploredReads->at(read1) == NODE_UNEXPLORED)
				{

//					FILE_LOG(logDEBUG) << "explore read (insert all edges of read) " << read1 << endl;
					insertAllEdgesOfRead(read1, exploredReads);	// Explore current node.
					if(loglevel > 3){
						print_neighbors(read1);
					}
					exploredReads->at(read1) = NODE_EXPLORED;
				}
				if(graph->at(read1)->size() != 0) 	// Read has some edges (required only for the first read when a new queue starts.
				{
					if(exploredReads->at(read1) == NODE_EXPLORED) 	// Explore unexplored neighbors first.
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++ )
						{

							UINT64 read2_number = (graph->at(read1)->at(index1)->getDestinationRead()).getReadNumber();
							if(exploredReads->at(read2_number) == NODE_UNEXPLORED) 	// Not explored.
							{
//								FILE_LOG(logDEBUG) << "explore read (insert all edges of read) " << read2_number << endl;
								queue->at(end++) = read2_number; 	// Put in the queue.
								insertAllEdgesOfRead(read2_number, exploredReads);
								if(loglevel > 3){
									print_neighbors(read2_number);
								}
								exploredReads->at(read2_number) = NODE_EXPLORED;
							}
						}
						markTransitiveEdges(read1, markedNodes); // Mark transitive edges
						exploredReads->at(read1) = NODE_EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
					}
					if(exploredReads->at(read1) == NODE_EXPLORED_AND_TRANSITIVE_EDGES_MARKED)
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++) 	// Then explore all neighbour's neighbors
						{
							UINT64 read2 = graph->at(read1)->at(index1)->getDestinationRead().getReadNumber();
							if(exploredReads->at(read2) == NODE_EXPLORED)
							{
								for(UINT64 index2 = 0; index2 < graph->at(read2)->size(); index2++) 	// Explore all neighbors neighbors
								{
									UINT64 read3 = graph->at(read2)->at(index2)->getDestinationRead().getReadNumber();
									if(exploredReads->at(read3) == NODE_UNEXPLORED) 	// Not explored
									{
//										FILE_LOG(logDEBUG) << "explore read (insert all edges of read) " << read3 << endl;
										queue->at(end++) = read3; 	// Put in the queue
										insertAllEdgesOfRead(read3, exploredReads);
										if(loglevel > 3){
											print_neighbors(read3);
										}
										exploredReads->at(read3) = NODE_EXPLORED;
									}
								}
								markTransitiveEdges(read2, markedNodes); // Mark transitive edge
								exploredReads->at(read2) = NODE_EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
							}
						}
						removeTransitiveEdges(read1); // Remove the transitive edges
					}
				}
				if(counter%100000==0)	// Show the progress.
				{
					cout <<"counter: " << setw(10) << counter << " Nodes: " << setw(10) << numberOfNodes << " Edges: " << setw(10) << numberOfEdges/2;
				}
			}
		}
	}

	delete exploredReads;
	delete queue;
	delete markedNodes;
	delete hashTable;	// Do not need the hash table any more.
//	for(UINT64 i = 1; i <= dataSet->numberOfUniqueReads; ++i){
//		print_neighbors(i);
//	}
	if (numberOfEdges > 0)
	{
		do
		{
			counter = contractCompositeEdges();
		} while (counter > 0);
	}

	CLOCKSTOP;
}		/* -----  end of method TOverlapGraph<T>::buildOverlapGraph  ----- */


	template < typename T >
void tOverlapGraph<T>::insertAllEdgesOfRead ( UINT64 readNumber, vector<nodeType> * exploredReads )
{
//	CLOCKSTART;
	tRead<T> read1 = dataSet->getReadFromID(readNumber);
//	FILE_LOG(logDEBUG) << "Getting edges of read with number " 
//		<< readNumber << ": " << read1 << endl;

	// Look up all the substring with length minOvl, starting from the second element
	// until last one
	for(UINT64 j = 1; j < read1.size()-minOvl; j++)
	{
		vector<T> subString;

		subString = read1.substr(j,minOvl);

		FILE_LOG(logDEBUG1) << "Look up substring " << subString << " in hashTable." << endl;

		// TODO: change this to "find" using iterator
		if (hashTable->count(subString) > 0) {
			vector<UINT64> listOfReads=hashTable->at(subString);

			if(loglevel > 3){
				for(size_t i = 0; i < listOfReads.size(); ++i){
					FILE_LOG(logDEBUG1) << ((listOfReads.at(i)) & 0X3FFFFFFFFFFFFFFF)  << " ";
				}
				FILE_LOG(logDEBUG1) << endl;
			}

			for(UINT64 k = 0; k < listOfReads.size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads.at(k);
				UINT64 overlapOffset;
				char orientation;
				tRead<T> &read2 = dataSet->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
				FILE_LOG(logDEBUG1) << "read 2 has read number " << (data & 0X3FFFFFFFFFFFFFFF) << " and orientation " << (data >> 62) << endl;
				if(exploredReads->at(read2.getReadNumber())!= NODE_UNEXPLORED) 	// No need to discover the same edge again. All edges of read2 is already inserted in the graph.
					continue;
				if(read1.superID == 0 && read2.superID == 0 && checkOverlap(read1,read2,(data >> 62),j)) // Both read need to be non contained.
				{
					switch (data >> 62) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
					{
						case 0: orientation = 3; overlapOffset = read1.size() - j; break; 	// 3 = r1>------->r2
						case 1: orientation = 0; overlapOffset = minOvl + j; break; 	// 0 = r1<-------<r2
						case 2: orientation = 2; overlapOffset = read1.size() - j; break; 	// 2 = r1>-------<r2
						case 3: orientation = 1; overlapOffset = minOvl + j; break; 	// 1 = r2<------->r2
					}
					insertEdge(read1,read2,orientation,read1.size()-overlapOffset); 	// Insert the edge in the graph.
				}
			}
		}
	}
	if(graph->at(readNumber)->size() != 0){
		sort(graph->at(readNumber)->begin(),graph->at(readNumber)->end(),compareEdges<tEdge<tRead<T>>>); // Sort the list of edges of the current node according to the overlap offset (ascending).
		if(loglevel > 3){
			print_neighbors(readNumber);
		}
	}
//	CLOCKSTOP;
}		/* -----  end of method tOverlapGraph<T>::insertAllEdgesOfRead  ----- */


	template < typename T >
bool tOverlapGraph<T>::checkOverlap ( const tRead<T> &read1, const tRead<T> &read2, UINT64 orient, UINT64 start )
{
//	FILE_LOG(logDEBUG) << "check overlap with orientation " << orient << " starting from position " << start << endl;
//	FILE_LOG(logDEBUG) << "between \n" << read1 << " and \n" << read2 << endl;
	if(orient == 0 || orient == 2)		
	// orient 0
	//   >--------MMMMMMMMMMMMMMM*************> 		read1      M means match found by hash table
	//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match
	//				OR
	// orient 2
	//	 >---*****MMMMMMMMMMMMMMM*************>		read1
	// 	          MMMMMMMMMMMMMMM*************-------<	Reverese complement of read2
	{
		if(read1.size() - start - minOvl >= read2.size()  - minOvl) // The overlap must continue till the end.
			return false;
		size_t check_len = read1.size() - (start + minOvl);
//		FILE_LOG(logDEBUG) << (read1.substr(start + minOvl, check_len, 0)) << endl;
//		FILE_LOG(logDEBUG) << (read2.substr(minOvl,  check_len, (orient >> 1))) << endl;
		bool overlap = (read1.substr(start + minOvl, check_len, 0) == read2.substr(minOvl,  check_len, (orient >> 1)));
//		FILE_LOG(logDEBUG) << std::boolalpha << "overlap? " << overlap << endl;
		return overlap; // If the remaining strings match.
	}
	else								
	// orient 1
	//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
	//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match
	//				OR
	// orient 3
	//	 	>********MMMMMMMMMMMMMMM-------------> 		read1
	//	    <----********MMMMMMMMMMMMMMM			Reverse Complement of Read2
	{
		if(read2.size() - minOvl < start)
			return false;
		size_t check_len = start;
//		FILE_LOG(logDEBUG) << (read1.substr(0, check_len, 0)) << endl;
//		FILE_LOG(logDEBUG) << (read2.substr(read2.size()-minOvl-start, check_len, (orient >> 1))) << endl;
		bool overlap = (read1.substr(0, check_len, 0) == read2.substr(read2.size()-minOvl-start, check_len, (orient >> 1)));
//		FILE_LOG(logDEBUG) << std::boolalpha << "overlap? " << overlap << endl;
		return overlap; // If the remaining strings match.
	}
}


	template < typename T >
void tOverlapGraph<T>::removeTransitiveEdges ( UINT64 readNumber )
{
	UINT64 counter(0);
	for(UINT64 index = 0; index < graph->at(readNumber)->size(); index++)  	// Go through the list of edges of the current read.
	{
		if(graph->at(readNumber)->at(index)->transitiveRemovalFlag == true)	// This edge is marked as transitive. We will first remove the reverese edge.
		{
			tEdge<tRead<T>> *twinEdge = graph->at(readNumber)->at(index)->getReverseEdge();
			++counter;
			UINT64 ID = twinEdge->getSourceRead().getReadNumber();
			for(UINT64 index1 = 0; index1 < graph->at(ID)->size(); index1++) 	// Get the reverse edge first
			{
				if((graph->at(ID)->at(index1)) == twinEdge)
				{
					delete twinEdge;
					graph->at(ID)->at(index1) = graph->at(ID)->at(graph->at(ID)->size()-1); // Move the transitive edges at the back of the list and remove.
					graph->at(ID)->pop_back();
					if(graph->at(ID)->empty())
						numberOfNodes--;
					numberOfEdges--;
					break;
				}
			}
		}
	}
	UINT64 j=0;
	for(UINT64 index=0; index < graph->at(readNumber)->size(); index++) // Then we will remove all the transitive edges of the current read.
	{
		if(graph->at(readNumber)->at(index)->transitiveRemovalFlag == false)	// We move all the non-transitive edges at the beginning of the list
			graph->at(readNumber)->at(j++) = graph->at(readNumber)->at(index);
		else
		{
			numberOfEdges--;
			delete graph->at(readNumber)->at(index);
		}
	}
	graph->at(readNumber)->resize(j);
	if(graph->at(readNumber)->empty())
		numberOfNodes--;
//	FILE_LOG(logDEBUG) << "Total number of transitive edges removed: " << counter << endl;
}


	template < typename T >
void tOverlapGraph<T>::insertEdge( tEdge<tRead<T>> * edge )
{
	UINT64 ID = edge->getSourceRead().getReadNumber(); // This is the source read.
	if(graph->at(ID)->empty()) 	// If there is no edge incident to the node
		numberOfNodes++;	// Then a new node is inserted in the graph. Number of nodes increased.
	graph->at(ID)->push_back(edge);	// Insert the edge in the list of edges of ID
	numberOfEdges++;	// Increase the number of edges.
}


	template < typename T >
void tOverlapGraph<T>::insertEdge(const tRead<T> &read1, const tRead<T> &read2, char orient, UINT64 overlapOffset)
{
	tEdge<tRead<T>> * edge1 = new tEdge<tRead<T>>(read1,read2,orient,overlapOffset);	// Create a new edge in the graph to insert.
	UINT64 overlapOffsetReverse = read2.size() + overlapOffset - read1.size();	// Set the overlap offset accordingly for the reverse edge. Note that read lengths are different.
	// If read lengths are the same. Then the reverse edge has the same overlap offset.
	tEdge<tRead<T>> * edge2 = new tEdge<tRead<T>>(read2,read1,edge1->getReverseOrientation(),overlapOffsetReverse);	// Create a new edge for the reverses string.

	edge1->setReverseEdge(edge2);	// Set the reverse edge pointer.
	edge2->setReverseEdge(edge1);	// Set the reverse edge pinter.
	insertEdge(edge1);	// Insert the edge in the overlap graph.
	insertEdge(edge2);	// Insert the edge in the overlap graph.
}


	template < typename T >
void tOverlapGraph<T>::removeEdge ( tEdge<tRead<T>> * edge )
{
	tEdge<tRead<T>> *twinEdge = edge->getReverseEdge();
	UINT64 ID1 = edge->getSourceRead().getReadNumber(), ID2 = edge->getDestinationRead().getReadNumber();  // Get the source and destation read IDs.
	for(UINT64 i = 0; i< graph->at(ID2)->size(); i++) // Delete the twin edge first.
	{
		if((graph->at(ID2)->at(i)) == twinEdge)
		{
			graph->at(ID2)->at(i) = graph->at(ID2)->at(graph->at(ID2)->size()-1);
			graph->at(ID2)->pop_back();
			if(graph->at(ID2)->empty())
				numberOfNodes--;
			numberOfEdges--;
			break;
		}
	}
	for(UINT64 i = 0; i< graph->at(ID1)->size(); i++) // Delete the edge then.
	{
		if((graph->at(ID1)->at(i)) == edge)
		{
			graph->at(ID1)->at(i) = graph->at(ID1)->at(graph->at(ID1)->size()-1);
			graph->at(ID1)->pop_back();
			if(graph->at(ID1)->empty())
				numberOfNodes--;
			numberOfEdges--;
			break;
		}
	}
}


	template < typename T >
void tOverlapGraph<T>::markTransitiveEdges(UINT64 readNumber, vector<markType> * markedNodes)
{
	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++) // Mark all the neighbours of the current read as INPLAY
		markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead().getReadNumber()) = MARK_INPLAY; // Inplay

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++) // Traverse through the list of edges according to their overlap offset.
	{
		UINT64 read2 = graph->at(readNumber)->at(i)->getDestinationRead().getReadNumber(); // For each neighbor
		if(markedNodes->at(read2) == MARK_INPLAY) 										// If the neighbor is marked as MARK_INPLAY
		{
			for(UINT64 j = 0; j < graph->at(read2)->size(); j++)
			{
				UINT64 read3 = graph->at(read2)->at(j)->getDestinationRead().getReadNumber(); // Get the neighbors neighbors
				if(markedNodes->at(read3) == MARK_INPLAY)
				{

					char type1 = graph->at(readNumber)->at(i)->getOrientation();
					char type2 = graph->at(read2)->at(j)->getOrientation();
					if((type1 == 0 ||  type1 == 2) && (type2==0 || type2==1)) 	// Check edge orientation
						markedNodes->at(read3) = MARK_ELIMINATED; 					// Mark as MARK_ELIMINATED
					else if((type1==1||type1==3) && (type2==2 || type2==3)) 	// Check edge orientation
						markedNodes->at(read3) = MARK_ELIMINATED; 					// Mark as MARK_ELIMINATED
				}
			}
		}
	}
	for(UINT64 i = 0;i < graph->at(readNumber)->size(); i++)
	{
		if(markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead().getReadNumber()) == MARK_ELIMINATED) // Current read to a node marked as MARK_ELIMINATED
		{
			graph->at(readNumber)->at(i)->transitiveRemovalFlag = true; 					// Mark this edge as transitive edge. Will remove this edge later.
			graph->at(readNumber)->at(i)->getReverseEdge()->transitiveRemovalFlag = true;	// Mark also the reverse edge. Will remove this edge later.
		}
	}

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++)
		markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead().getReadNumber()) = MARK_VACANT; 	// Change back all the variables modified in this function to MARK_VACANT

	markedNodes->at(readNumber) = MARK_VACANT; 		// Mark as vacant.
}


	template < typename T >
UINT64 tOverlapGraph<T>::contractCompositeEdges(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	for(UINT64 index = 1 ; index < graph->size(); index++)
	{
		if(graph->at(index)->size() == 2) // Check if the node has only two edges.
		{
			tEdge<tRead<T>> * edge1 = graph->at(index)->at(1)->getReverseEdge();  // Second Edge.
			tEdge<tRead<T>> * edge2 = graph->at(index)->at(0);  // First edge.
			{
				if( mergeable(*edge1, *edge2) && edge1->getSourceRead().getReadNumber() != edge1->getDestinationRead().getReadNumber()) // One incoming edge and one outgoing edge.
				{
					tEdge<tRead<T>> * merged_edge = *edge1 + *edge2;
					tEdge<tRead<T>> * merged_edge_reverse = *(edge2->getReverseEdge()) + *(edge1->getReverseEdge());
					merged_edge->setReverseEdge(merged_edge_reverse);
					merged_edge_reverse->setReverseEdge(merged_edge);

					removeEdge(edge1);
					removeEdge(edge2);
					insertEdge(merged_edge);
					insertEdge(merged_edge_reverse);
					++counter;	// Counter how many edges merged.
				}
			}
		}

	}
	cout << setw(10) << counter << " composite Edges merged." << endl;
	CLOCKSTOP;
	return counter;
}

	template < typename T >
void tOverlapGraph<T>::print_neighbors(UINT64 readNumber) const
{
	cout << graph->at(readNumber)->size() << " neighbors of read with number " << readNumber << ":" << endl;
	for(UINT64 i = 0; i< graph->at(readNumber)->size(); ++i){
		cout << *(graph->at(readNumber)->at(i)) << endl;
	}
}

	template < typename T >
vector<T> tOverlapGraph<T>::getStringInEdge (const tEdge<tRead<T>> *edge)
{
	FILE_LOG(logDEBUG) << *edge << endl;
	vector<T> returnString, dest_string;
	char source_orient = ((edge->getOrientation()) >> 1);
	char dest_orient = ((edge->getOrientation()) & 1);
	if(dest_orient == 0){
		dest_string = edge->getDestinationRead().getStringReverse();
	}
	else{
		dest_string = edge->getDestinationRead().getStringForward();
	}

	if (edge->getListOfOverlapOffsets()->size() == 0){
		returnString = edge->getSourceRead().substr(0, 
				edge->getOverlapOffset(),
				(1-source_orient));
	}
	else {
		returnString = edge->getSourceRead().substr(0, 
				edge->getListOfOverlapOffsets()->at(0),
				(1-source_orient));
		UINT64 last_offset = edge->getOverlapOffset() - edge->getListOfOverlapOffsets()->at(0);
		UINT64 i = 0;
		for(i = 0; i < (edge->getListOfReadIDs()->size()-1); ++i){
//			FILE_LOG(logDEBUG) << "Read with ID " << edge->getListOfReadIDs()->at(i)
//				<< " : " << dataSet->getReadFromID(edge->getListOfReadIDs()->at(i)) << endl;
//			FILE_LOG(logDEBUG) << "Orientation: " << edge->getListOfOverlapOffsets()->at(i)
//				<< " and overlap offset: " << edge->getListOfOverlapOffsets()->at(i) << endl;
			vector<T> new_string = (dataSet->getReadFromID(edge->getListOfReadIDs()->at(i))).substr(0, 
					edge->getListOfOverlapOffsets()->at(i+1),
					(1-edge->getListOfOrientations()->at(i)));
			for(auto it = new_string.cbegin(); it != new_string.cend(); ++it){
				returnString.emplace_back(*it);
			}
			last_offset -= edge->getListOfOverlapOffsets()->at(i+1);
		}
		// Last read/overlap offset
		vector<T> new_string = (dataSet->getReadFromID(edge->getListOfReadIDs()->at(i))).substr(0, 
				last_offset, (1-edge->getListOfOrientations()->at(i)));
		for(auto it = new_string.cbegin(); it != new_string.cend(); ++it){
			returnString.emplace_back(*it);
		}
	}
	for(auto it = dest_string.cbegin(); it!=dest_string.cend(); ++it){
		returnString.emplace_back(*it);
	}
	return returnString;
}

	template < typename T >
bool tOverlapGraph<T>::getEdges(vector<tEdge<tRead<T>> *> & contigEdges)
{
	CLOCKSTART;
	contigEdges.clear();
	for(UINT64 i = 1; i<= dataSet->numberOfUniqueReads; i++)
	{
		if(!graph->at(i)->empty()) // if this read has some edge(s) going out of it (since now the graph is directed)
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)
			{
				tEdge<tRead<T>> * e 		= graph->at(i)->at(j);
				UINT64 source 			= e->getSourceRead().getReadNumber(); 
				UINT64 destination 		= e->getDestinationRead().getReadNumber();
				if(source < destination || (source == destination && e < e->getReverseEdge()) )  /* Only consider the edges between non-contained reads */
				{
					contigEdges.push_back(e); // List of contigs.
					//e->setEndCorrdinateLimit(dataSet->getPacBioReadLength());
				}
			}
		}
	}
	FILE_LOG(logINFO) << "Number of edges found in graph is " << contigEdges.size() << endl;
	CLOCKSTOP;
	return true;
}


	template < typename T >
vector<vector<T>>* tOverlapGraph<T>::getContigs(const vector<tEdge<tRead<T>> *> & contigEdges, bool isolates)
{
	CLOCKSTART;
	vector<vector<T>>* contig_strings = new vector<vector<T>>;
	for(UINT64 i = 1; i <= dataSet->numberOfUniqueReads; i++){
		// For isolate reads that are not contained in other reads
		if(graph->at(i)->empty() && dataSet->getReadFromID(i).superID == 0 && isolates){
			contig_strings->push_back(dataSet->getReadFromID(i).getStringForward());
		}
	}
	if(isolates){
		FILE_LOG(logINFO) << "Number of strings reads: " 
			<< contig_strings->size() << endl;
	}
	for(UINT64 i = 0; i < contigEdges.size(); ++i){
		contig_strings->push_back(getStringInEdge(contigEdges.at(i)));
	}
	FILE_LOG(logINFO) << "Number of strings from contigs and reads in total: " 
		<< contig_strings->size() << endl;
	CLOCKSTOP;
	return contig_strings;
}
	template < typename T >
void tOverlapGraph<T>::printGraph(string graphFileName, const vector<tEdge<tRead<T>> *> & contigEdges)
{
	CLOCKSTART;
	UINT64 thickness;
	string edgeColor;
	ofstream graphFilePointer; 

	/************************* Store the graph in a file. ************************/
	graphFilePointer.open(graphFileName.c_str());
	if(!graphFilePointer.is_open())
		MYEXIT("Unable to open file: "+graphFileName);

	// Graph specification before the nodes and edges
	graphFilePointer << "graph: {" << endl 
		<<  "layoutalgorithm :forcedir" << endl 
		<<  "fdmax:704" << endl 
		<<  "tempmax:254" << endl 
		<<  "tempmin:0" << endl 
		<<  "temptreshold:3" << endl 
		<<  "tempscheme:3" << endl 
		<<  "tempfactor:1.08" << endl 
		<<  "randomfactor:100" << endl 
		<<  "gravity:0.0" << endl 
		<<  "repulsion:161" << endl 
		<<  "attraction:43" << endl 
		<<  "ignore_singles:yes" << endl 
		<<  "node.fontname:\"helvB10\"" << endl 
		<<  "edge.fontname:\"helvB10\"" << endl 
		<<  "node.shape:box" << endl 
		<<  "node.width:80" << endl 
		<<  "node.height:20" << endl 
		<<  "node.borderwidth:1" << endl 
		<<  "node.bordercolor:31" << endl;

	// All the nodes, title and label
	for(UINT64 i = 1; i<= dataSet->numberOfUniqueReads; i++)
	{
		// Print nodes even if there are no edge connected to it
		if(!graph->at(i)->empty()){
			graphFilePointer << "node: { title:\""<< i 
				<<"\" label: \"" << i << "\" }" << endl;
		}
	}

	// All the edges
	for (UINT64 i = 0; i < contigEdges.size(); i++)
	{
		tEdge<tRead<T>> * e = contigEdges.at(i);
		thickness = e->getListOfReadIDs()->empty() ? 1: 3;
		UINT64 source = e->getSourceRead().getReadNumber(), destination = e->getDestinationRead().getReadNumber();
		// Edge label: (first overlap length, edge length, number of reads, overlap offset, last overlap length)
		if(source < destination || (source == destination && e < e->getReverseEdge()) )
		{
			if(e->getOrientation() == 0)
				graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination
					<< "\" thickness: " << thickness << " arrowstyle: none backarrowstyle: solid color: red label: \"(" 
					<< e->getListOfReadIDs()->size() << "," << e->getOverlapOffset()
					<< ")\" }" << endl;
			else if(e->getOrientation() == 1)
				graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
					<< thickness << " backarrowstyle:solid arrowstyle:solid color: green label: \"(" 
					<< e->getListOfReadIDs()->size() << "," << e->getOverlapOffset()
					<< ")\" }" << endl;
			else if(e->getOrientation() == 2)
				graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
					<< thickness << " arrowstyle: none color: blue label: \"(" 
					<< e->getListOfReadIDs()->size() << "," << e->getOverlapOffset()
					<< ")\" }" << endl;
			else if(e->getOrientation() == 3)
				graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: "
					<< thickness << " arrowstyle:solid color: red label: \"(" 
					<< e->getListOfReadIDs()->size() << "," << e->getOverlapOffset()
					<< ")\" }" << endl;
		}
	}
	graphFilePointer << "}";
	graphFilePointer.close();
	CLOCKSTOP;
	/************************* Store the graph in a file done. ************************/
}

#endif
