#ifndef READ_H
#define READ_H

//============================================================================
// Name        : Read.h
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Read header file
//============================================================================

#include "Config.h"
#include<utility> // for std::pair

class Edge;

enum EdgeOP {INSERTION, DELETION};

typedef std::pair<Edge* , UINT32> t_edge_loc_pair;
class Read
{
	private:
		/* ====================  DATA MEMBERS  ======================================= */
		UINT64 m_readID;         // Unique Identification of the read.

		std::string m_seq;                   // String representation of the read.

		std::string m_rev_seq;

		UINT32 m_length;         // Length of the read.

		/* vector of pairs, with each pair storing 
		 * first: pointer to an edge that includes this read
		 * second: index(location) of this read on the edge
		 * One vector stores the forward-string edges, 
		 * and another vector stores the backward-string edges
		 */
		vector< t_edge_loc_pair > * m_fwdEdges;

		vector< t_edge_loc_pair > * m_bwdEdges;

		/* ====================  METHODS      ======================================= */
		void initEdgeInfo();

	public:
		/* ====================  LIFECYCLE     ======================================= */
		Read();

		Read(const UINT32 length, const std::string seq);

		Read(const UINT32 length);

		Read(const std::string & seq);

		Read(const Read &s_read);

		~Read(void);

		/* ====================  OPERATORS     ======================================= */
		char & operator[](size_t i){assert(i < m_seq.length()); return m_seq[i];}

		Read& operator=(const Read &s_read);

		/* ====================  DATA MEMBERS  ======================================= */
		friend class Edge;

		/* ====================  MUTATORS      ======================================= */
		void setReadID(const UINT64 & ID){m_readID = ID;}

		void setSeq(const std::string & seq);

		void setReadLength(UINT32 length){m_length = length;}

		void updateEdgeInfo(Edge *edge, UINT32 read_index, EdgeOP);

		/* ====================  ACCESSORS     ======================================= */ 
		std::string getStringForward(void) const {return m_seq;}

		std::string getStringReverse(void) const {return m_rev_seq;}

		UINT32 getReadLength(void) const {return m_length;}

		UINT64 getReadID(void) const {return m_readID;}

		vector< t_edge_loc_pair >* getFwdEdges() const {return m_fwdEdges;}

		vector< t_edge_loc_pair >* getBwdEdges() const {return m_bwdEdges;}

};/* End of class Read */

std::string reverseComplement(const std::string & seq);

std::ostream &operator<<(std::ostream & out, const Read & read);

#endif /* READS_H */
