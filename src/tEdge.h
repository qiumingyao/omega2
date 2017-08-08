/*
 * =====================================================================================
 *
 *       Filename:  tEdge.h
 *
 *    Description:  template Edge class
 *
 *        Version:  1.0
 *        Created:  07/23/2015 17:58:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  JJ Crosskey (cjg), jjchai01@gmail.com
 *   Organization:  ORNL
 *
 * =====================================================================================
 */

#ifndef TEDGE_H_
#define TEDGE_H_

#include "Config.h"
#include "tRead.h"

template <typename T>
class tEdge
{
	private:
		T  source;

		T  destination;

		char orientation;

		UINT64 overlapOffset;

		vector<UINT64> * listOfReadIDs;

		vector<UINT64> * listOfOverlapOffsets;

		vector<char> * listOfOrientations;

		tEdge * reverseEdge;

		INT64 edgeID;
	public:
		// variables
		bool transitiveRemovalFlag;
		// Constructor and Destructor
		tEdge(void);

		tEdge(const T & from, const T & to, char orient, UINT64 offset);

		tEdge(const T & from, const T & to, char orient, UINT64 offset,
			vector<UINT64> * listOfReadIDs, vector<UINT64> * list_offsets,
			vector<char> * list_orients);

		~tEdge();

		// Accessors
		T getSourceRead() const {return source;}

		T getDestinationRead() const {return destination;}

		char getOrientation() const {return orientation;}

		UINT64 getOverlapOffset() const {return overlapOffset;}

		vector<UINT64> * getListOfReadIDs() const {return listOfReadIDs;}

		vector<UINT64> * getListOfOverlapOffsets() const {return listOfOverlapOffsets;}

		vector<char> * getListOfOrientations() const {return listOfOrientations;}

		tEdge * getReverseEdge() const {return reverseEdge;}

		char getReverseOrientation() const;

		bool operator<(const tEdge<T> &rhs) const {return (overlapOffset < rhs.overlapOffset);}

		// Mutators
		void setEdgeID(INT64 id) {edgeID = id;}

		void setReverseEdge(tEdge * edge){reverseEdge = edge;}

		tEdge * makeReverseEdge();

};
template <typename T>
tEdge<T>::tEdge(void)
{
	/* initialize orientation to -1 (no sense) */
	orientation 		= -1;                       
	overlapOffset		= 0;
	listOfReadIDs 		= new vector<UINT64>;
	listOfOverlapOffsets 	= new vector<UINT64>;
	listOfOrientations 	= new vector<char>;
	reverseEdge 		= nullptr;
	edgeID 			= 0;
	transitiveRemovalFlag 	= false;
}

template <typename T>
tEdge<T>::~tEdge(void)
{
//	cout << "deconstructor of edge" << endl;
//	cout << *this << endl;
	if (listOfReadIDs != nullptr)
		delete listOfReadIDs;
	if (listOfOverlapOffsets != nullptr)
		delete listOfOverlapOffsets;
	if (listOfOrientations != nullptr)
		delete listOfOrientations;
//	cout << "deconstructor of edge finished" << endl;
}

template <typename T>
tEdge<T>::tEdge(const T & from, const T & to, char orient, UINT64 offset)
{
	source 			= from;
	destination		= to;
	orientation		= orient;
	overlapOffset		= offset;
	listOfReadIDs 		= new vector<UINT64>;
	listOfOverlapOffsets 	= new vector<UINT64>;
	listOfOrientations 	= new vector<char>;
	reverseEdge 		= nullptr;
	edgeID 			= 0;
	transitiveRemovalFlag	= false;
}

template <typename T>
tEdge<T>::tEdge(const T &from, const T &to, char orient, UINT64 offset,
			vector<UINT64> * list_readIDs, vector<UINT64> * list_offsets,
			vector<char> * list_orients)
{
	source 			= from;
	destination		= to;
	orientation		= orient;
	overlapOffset		= offset;
	listOfReadIDs 		= list_readIDs;
	listOfOverlapOffsets 	= list_offsets;
	listOfOrientations 	= list_orients;
	reverseEdge 		= nullptr;
	edgeID 			= 0;
	transitiveRemovalFlag	= false;
}

template <typename T>
char tEdge<T>::getReverseOrientation() const
{
	char r_orient = ~orientation;
	r_orient = ((r_orient & 1) << 1) + ((r_orient & 2) >> 1);
	return r_orient;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  operator+
 *  Description:  merge two edges together, if possible, otherwise return null.
 *  		  This function does not take care of the reverse edge merging,
 *  		  need to call this function again to merge the reverse edges
 * =====================================================================================
 */
template <typename T>
tEdge<T>* operator+(const tEdge<T> & lhs, const tEdge<T> & rhs)
{
	if (mergeable( lhs, rhs)){
		T m_source			= lhs.getSourceRead();
		T m_destination			= rhs.getDestinationRead();
		char m_orient			= ((lhs.getOrientation() & 2) | (rhs.getOrientation() & 1));
		UINT64 m_offset			= lhs.getOverlapOffset() + rhs.getOverlapOffset();

		vector<UINT64> * list_readIDs	= new vector<UINT64>;
		vector<UINT64> * list_offsets	= new vector<UINT64>;
		vector<char> *list_orients	= new vector<char>;

		UINT64 offset 			= 0;
		for(UINT64 i = 0; i < lhs.getListOfOrientations()->size(); i++){
			list_readIDs->push_back(lhs.getListOfReadIDs()->at(i));
			list_offsets->push_back(lhs.getListOfOverlapOffsets()->at(i));
			list_orients->push_back(lhs.getListOfOrientations()->at(i));
			offset += lhs.getListOfOverlapOffsets()->at(i);
		}
		list_readIDs->push_back(rhs.getSourceRead().getReadNumber());
		list_offsets->push_back(lhs.getOverlapOffset() - offset);
		list_orients->push_back((lhs.getOrientation() & 1));
		for(UINT64 i = 0; i < rhs.getListOfOrientations()->size(); i++){
			list_readIDs->push_back(rhs.getListOfReadIDs()->at(i));
			list_offsets->push_back(rhs.getListOfOverlapOffsets()->at(i));
			list_orients->push_back(rhs.getListOfOrientations()->at(i));
		}
		tEdge<T> * merged_edge = new tEdge<T>(m_source, m_destination, m_orient, m_offset,
				list_readIDs, list_offsets, list_orients);
		return merged_edge;
	}
	else
		return nullptr;
}

template <typename T>
bool mergeable(const tEdge<T> & lhs, const tEdge<T> & rhs)
{
	if(lhs.getDestinationRead().getReadNumber() != rhs.getSourceRead().getReadNumber())
		return false;
	else if((lhs.getOrientation() & 1) != ((rhs.getOrientation() & 2) >> 1)) 
		return false;
	else
		return true;
}

template <typename T>
tEdge<T> * tEdge<T>::makeReverseEdge()
{
	tEdge<T> * r_edge 	= new tEdge<T>;
	r_edge->source 		= destination;
	r_edge->destination 	= source;
	r_edge->orientation 	= getReverseOrientation();
	r_edge->overlapOffset 	= overlapOffset + destination->size() - source->size();

	// Reverse the list of reads
	for(UINT64 i = listOfReadIDs->size() - 1; i <= 0; --i){
		r_edge->listOfReadIDs->push_back(listOfReadIDs->at(i));
	}
	// Reverse the list of overlap offset
//	UINT64 r_edge_first_ovl = r_edge->overlapOffset;
//	vector<UINT64> r_edge_overlap_offsets;
//	r_edge_overlap_offsets.push_back(listOfOverlapOffsets->at(0) + source->size() - getReadFromID(listOfReadIDs->at(0))->size());
//	r_edge_first_ovl -= r_edge_overlap_offsets.at(0);
//	for(UINT64 i = 0; i < listOfOverlapOffsets->size(); ++i){
//		UINT64 ovl = listOfOverlapOffsets->at(i+1) + getReadFromID(listOfReadIDs->at(i))->size() - 
//			getReadFromID(listOfReadIDs->at(i+1))->size();
//		r_edge_overlap_offsets.push_back(ovl);
//		r_edge_first_ovl -= ovl;
//	}
//	for(UINT64 i = r_edge_overlap_offsets.size() - 1; i < 1; --i){
//		r_edge->listOfOverlapOffsets->push_back(r_edge_overlap_offsets.at(i));
//	}
	// Reverse the list of orientations
	for(UINT64 i = listOfOrientations->size() - 1; i <= 0; --i){
		r_edge->listOfOrientations->push_back(~(listOfOrientations->at(i)));
	}
	r_edge->reverseEdge = this;
	reverseEdge = r_edge;
	r_edge->edgeID = -edgeID;
	return r_edge;
}

inline string get_arrow(char orientation) { return ((orientation & 1) ? ">" : "<"); }

template <typename T>
ostream& operator<<(ostream & out, const tEdge<T> & edge)
{
	out << std::left << setw(10) << setfill(' ') << (edge.getSourceRead()).getReadNumber() 
		<< " " << get_arrow(edge.getOrientation() >> 1) << "--"
		<< setw(4) << std::right << setfill('0') << edge.getOverlapOffset() << "--" 
		<< get_arrow(edge.getOrientation() & 1) << " "
		<< std::right << setw(10) << setfill(' ') << (edge.getDestinationRead()).getReadNumber();
	if((edge.getListOfReadIDs()->size()) > 0){
		out << " includes: ";
		for(auto it = edge.getListOfReadIDs()->cbegin(); it != edge.getListOfReadIDs()->cend(); ++it){
			out << " " << *it;
		}
	}
	return out;
}

#endif // TEDGE_H_
