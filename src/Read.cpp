//============================================================================
// Name        : Read.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Read cpp file
//============================================================================


#include "Config.h"
#include "Read.h"


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initEdgeInfo
 *  Description:  Initialize vector of pairs that stores the information of edges
 *  		  in which the current read is included.
 * =====================================================================================
 */
void Read::initEdgeInfo()
{
	m_fwdEdges = new vector< t_edge_loc_pair >;
	m_bwdEdges = new vector< t_edge_loc_pair >;
}

Read::Read()
{
	m_readID	= 0;
	m_seq		= std::string();
	m_rev_seq	= std::string();
	m_length	= 0;
	initEdgeInfo();
}

Read::Read(const UINT32 length, const std::string seq)
{
	assert(length == seq.length());
	m_length = length;
	m_seq = seq;
	m_rev_seq = reverseComplement(seq);
	initEdgeInfo();
}

Read::Read(const UINT32 length)
:  m_length(length)
{
	initEdgeInfo();
}


Read::Read(const std::string & seq)
	: m_seq(seq)
{
	m_length = seq.length();
	m_rev_seq = reverseComplement(seq);
	initEdgeInfo();
}


// Technically, this is not needed, since it's just element-wise
// shallow copy.
// But this sets up the copy and move constructor (assignment) for the future
// if any memory allocation is added for the Read class.
Read::Read(const Read &s_read)
	: m_readID(s_read.m_readID),
	m_seq(s_read.m_seq),
	m_rev_seq(s_read.m_rev_seq),
	m_length(s_read.m_length)
{
	m_fwdEdges = new vector< t_edge_loc_pair >(*(s_read.m_fwdEdges));
	m_bwdEdges = new vector< t_edge_loc_pair >(*(s_read.m_bwdEdges));
}


Read& Read::operator=(const Read &s_read)
{
	if (this == &s_read)
		return *this;
	m_readID 	= s_read.m_readID;
	m_seq		= s_read.m_seq;
	m_rev_seq	= s_read.m_rev_seq;
	m_length	= s_read.m_length;
	if(!m_fwdEdges){
		delete m_fwdEdges;
		m_fwdEdges = nullptr;
	}
	if(!m_bwdEdges){
		delete m_bwdEdges;
		m_bwdEdges = nullptr;
	}
	m_fwdEdges = new vector< t_edge_loc_pair >(*(s_read.m_fwdEdges));
	m_bwdEdges = new vector< t_edge_loc_pair >(*(s_read.m_bwdEdges));
	return *this;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Destructor
 *  Description:  Destructor for Read object
 * =====================================================================================
 */
Read::~Read(void)
{
	if(!m_fwdEdges){
		delete m_fwdEdges;
		m_fwdEdges = nullptr;
	}
	if(!m_bwdEdges){
		delete m_bwdEdges;
		m_bwdEdges = nullptr;
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  setSeq
 *  Description:  If a read's sequence was not stored (only got length), set the nucleotide
 *  		  sequence.
 * =====================================================================================
 */
void Read::setSeq(const std::string & seq)
{
	if(seq.length() != m_length){
		FILE_LOG(logWARNING) << "Nucleotide sequence does not have the right string length "
			<< m_length << ", set to new length " << seq.length() << endl;
	}
	m_seq = seq; 
	m_rev_seq = reverseComplement(seq); 
	m_length=seq.length();
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  reverseComplement
 *  Description:  Returns the reverse complement of a read.
 * =====================================================================================
 */
std::string reverseComplement(const std::string & seq)
{
	UINT64 stringLength = seq.length();
	std::string reverse(stringLength,'0');
	for(UINT64 i = 0;i < stringLength; i++)
	{
		if( seq[i] & 0X02 ) // C or G
			reverse.at(stringLength -  i - 1 ) = seq[i] ^ 0X04;
		else // A <==> T
			reverse.at(stringLength -  i - 1 ) = seq[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
}

std::ostream &operator<<(std::ostream & out, const Read & read)
{
	out << "ID: " << setw(10) << setfill(' ') << read.getReadID()
		<< ", String: \n" << read.getStringForward();
	return out;
}
