#ifndef DATASET_H
#define DATASET_H

//============================================================================
// Name        : DataSet.h
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : DataSet header file
//============================================================================

#include "Config.h"
#include "Read.h"
typedef unordered_map<UINT64, UINT64> t_idmap;
class DataSet
{
	private:
		/* ====================  DATA MEMBERS  ======================================= */
		vector<Read*> *m_vec_reads;    /* vector pointers to Reads */

		t_idmap *m_readIDMap; // map from IDs in file to IDs in graph

		/* ====================  METHODS      ======================================= */
		// Load reads from a read file
		void loadReadsFromReadFile(const std::string &read_file, const bool read_string = false);
		//
		// Load reads from an edge file, in this case, only the length is available
		void loadReadsFromEdgeFile(const std::string &edge_file);

		bool do_sizes_match()
		{ return m_vec_reads->size() == m_readIDMap->size();}
	public:
		/* ====================  LIFECYCLE     ======================================= */
		DataSet();

		// Load a read file or an edge file
		DataSet(const std::string &input_file, const bool is_reads_file = true, 
				const bool read_string = false);

		// Load a list of files (either reads or edges)
		DataSet(const vector<std::string> &input_files, const bool is_reads_file = true, 
				const bool read_string = false);

		// Copy constructor
		DataSet(const DataSet &s_dataset);

		~DataSet();

		/* ====================  OPERATORS     ======================================= */
		friend std::ostream& operator<< (std::ostream &out, DataSet & a_data_set);

		DataSet& operator= (const DataSet &s_dataset);

		/* ====================  MUTATORS      ======================================= */
		void addRead(Read *r){m_vec_reads->push_back(r);}

		void rmRead(Read *r);

		/* ====================  ACCESSORS     ======================================= */ 
		UINT64 size() const{ return m_vec_reads->size();}

		Read* at(UINT64 ID) const;

		t_idmap* getReadIDMap() const {return m_readIDMap;}
};

bool compareEdgesByReads (const Edge *edge1, const Edge* edge2);

void printUnorderedMap(const unordered_map<UINT64, UINT64> & readIDMap, ostream & out);
#endif /* DATASET_H */
