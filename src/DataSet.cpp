//============================================================================
// Name        : DataSet.cpp
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : DataSet cpp file
//============================================================================

#include "DataSet.h"

DataSet::DataSet(void)
{
	m_vec_reads = new vector<Read *>;
	m_readIDMap = new t_idmap;
}

void DataSet::loadReadsFromReadFile(const std::string &read_file, const bool read_string)
{
	CLOCKSTART;
	FILE_LOG(logINFO) << "load reads from read file: " << read_file << "\n";
	UINT64 readID(m_vec_reads->size() + 1);

	// To count of reads in this file
	UINT64 readCount = 0;

	// Open file
	ifstream filePointer;
	filePointer.open(read_file.c_str());
	if(!filePointer.is_open()){
		FILE_LOG(logWARNING) << "Unable to open file: " << read_file << "\n";
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
				if(text[0] == '>'){
					FILE_LOG(logINFO) << "Input reads file format: FASTA\n";
					fileType = FASTA;
				}
				else if(text[0] == '@'){
					FILE_LOG(logINFO) << "Input reads file format: FASTA\n";
					fileType = FASTQ;
				}
				else{
					FILE_LOG(logERROR) << "Unknown input file format."<<endl;
					break;
				}
				filePointer.seekg(0, filePointer.beg);
			}
		}

		line.clear();

		// FASTA file read
		if(fileType == FASTA) {
			getline (filePointer,line0);	// get ID line
			getline (filePointer,line1,'>');	// get string line
			line1.erase(std::remove(line1.begin(), line1.end(), '\n'), 
					line1.end());
			//std::transform(line1.begin(), line1.end(), line1.begin(), ::toupper);
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

		// Map assign, key=readID from alinger, val=readID on graph (sequencial)
		istringstream readNameStream(readName);
		UINT64 readIDInFile;
		// Assuming that the read ID in the read file is an integer
		try{
			readNameStream >> readIDInFile;
			(*m_readIDMap)[readIDInFile] = readID;
		}
		catch(std::ios::failure e){
			FILE_LOG(logERROR) << "Read ID in the read file " << read_file 
				<< " is not integer.\n";
			exit(1);
		}

		Read *r = nullptr;

		if (read_string)	// Read the string to memory
			r = new Read(line1);
		else	// Only store the length
			r = new Read(line1.length());

		// readID, read(sequence), readLength
		r->setReadID(readID);

		// add read to the dataset
		addRead(r);

		++readID;
		++readCount;
		if(readCount % 1000000 == 0){
			FILE_LOG(logDEBUG) << setw(10) << (readCount/1000000)  << ",000,000"  
				<< " reads loaded to memory, "
				<< setw(7) << checkMemoryUsage() << " MB\n";
		}
	}

	filePointer.close();
	FILE_LOG(logDEBUG) << setw(10) << readCount << " reads loaded from this read file\n";
	CLOCKSTOP;
}

void DataSet::loadReadsFromEdgeFile(const std::string &edge_file)
{
	// TODO
	// edge files contain both read and edge information, they can be loaded
	// separately for reads/edges, or together for overlap graph
	FILE_LOG(logINFO) << "load reads from edge file: " << edge_file << "\n"
		<< "currently not implemented, please load from reads file instead.\n";
}

DataSet::DataSet(const std::string &input_file, const bool is_reads_file, const bool read_string)
{
	m_vec_reads = new vector<Read *>;
	m_readIDMap = new t_idmap;
	if (is_reads_file)
		loadReadsFromReadFile(input_file, read_string);
	else
		loadReadsFromEdgeFile(input_file);
	assert(do_sizes_match() == true);
}

DataSet::DataSet(const vector<std::string> &input_files, const bool is_reads_file, const bool read_string)
{
	CLOCKSTART;
	m_vec_reads = new vector<Read *>;
	m_readIDMap = new t_idmap;
	for(auto it = input_files.cbegin(); it != input_files.cend(); ++it){
		if (is_reads_file)
			loadReadsFromReadFile(*it, read_string);
		else
			loadReadsFromEdgeFile(*it);
	}
	assert(do_sizes_match() == true);
	CLOCKSTOP;
}

DataSet::DataSet(const DataSet &s_dataset)
{
	m_vec_reads = new vector<Read*>;
	m_vec_reads->reserve((s_dataset.m_vec_reads)->size());
	for(auto it = s_dataset.m_vec_reads->cbegin(); it != s_dataset.m_vec_reads->cend(); ++it){
		m_vec_reads->push_back(*it);
	}
	m_readIDMap = new t_idmap(*(s_dataset.m_readIDMap));
}

DataSet::~DataSet()
{
	if (m_vec_reads != nullptr){
		// The DataSet class should not delete reads
		for(UINT64 i = 0; i < m_vec_reads->size(); i++)
		{
			if (m_vec_reads->at(i) != nullptr){
				delete m_vec_reads->at(i);
				m_vec_reads->at(i) = nullptr;
			}
		}
		delete m_vec_reads;
		m_vec_reads = nullptr;
	}
	if (!m_readIDMap){
		delete m_readIDMap;
		m_readIDMap = nullptr;
	}
}

std::ostream& operator<< (std::ostream &out, DataSet & a_data_set)
{
	out << "Dataset with size " << a_data_set.size() << endl;
	return out;
}

DataSet& DataSet::operator= (const DataSet &s_dataset)
{
	if(this == &s_dataset)
		return *this;
	delete m_vec_reads;
	delete m_readIDMap;
	m_vec_reads = new vector<Read*>;
	m_vec_reads->reserve((s_dataset.m_vec_reads)->size());
	for(auto it = s_dataset.m_vec_reads->cbegin(); it != s_dataset.m_vec_reads->cend(); ++it){
		m_vec_reads->push_back(*it);
	}
	m_readIDMap = new t_idmap(*(s_dataset.m_readIDMap));
	return *this;
}

void DataSet::rmRead(Read *r)
{
	m_vec_reads->erase(std::remove(m_vec_reads->begin(), m_vec_reads->end(), r),m_vec_reads->end());
}


Read * DataSet::at(UINT64 ID) const
{
	assert(ID >= 1 && ID <= m_vec_reads->size());
	return m_vec_reads->at(ID - 1);
}


//=============================================================================
// print readIDMap (unordered map) key=readID from the aligner, val=readID in graph
//=============================================================================
void printUnorderedMap(const unordered_map<UINT64, UINT64> & readIDMap, ostream & out)
{
	out << "ID_in_file\tID_in_omega\n";
	for(unordered_map<UINT64, UINT64>::const_iterator it=readIDMap.begin(); it!=readIDMap.end(); ++it){
		out << it->first << "\t" << it->second << "\n";
	}
}

