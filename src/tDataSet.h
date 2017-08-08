/*
 * =====================================================================================
 *
 *       Filename:  tDataSet.h
 *
 *    Description:  template DataSet class
 *    		    The data set includes "reads" with typename T,
 *    		    which corresponds to "Read" in the original 
 *    		    omega implementation.
 *
 *        Version:  1.0
 *        Created:  07/23/2015 16:20:59
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  JJ Crosskey (cjg), jjchai01@gmail.com
 *   Organization:  ORNL
 *
 * =====================================================================================
 */

#ifndef TDATASET_H
#define TDATASET_H
#include "Config.h"
#include "tRead.h"

template <typename T>
class tDataSet
{
	private:
		vector<T> *reads;
		
		bool is_sorted;

		bool has_duplicates;
	public:
		tDataSet(void);
		
		~tDataSet();

		UINT64 numberOfUniqueReads;

		bool isSorted(void) const {return is_sorted;}

		bool isDupRemoved(void) const {return (!has_duplicates);}

		T& getReadFromID(UINT64 id);

		void addRead(T r){reads->push_back(r);}

		void sortReads(void);

		void rmDupReads(void);

		UINT64 size() const {return reads->size();}
};

template <typename T>
tDataSet<T>::tDataSet(void)
{
	reads = new vector<T>;
	numberOfUniqueReads = 0;
	is_sorted = false;
	has_duplicates = true;
}

template <typename T>
tDataSet<T>::~tDataSet(void)
{
//	for(UINT64 i = 0; i < reads->size(); i++)
//		delete reads->at(i);
	cout << "deconstructor of tDataSet called " << endl;
	delete reads;
}

template <typename T>
T& tDataSet<T>::getReadFromID(UINT64 ID)
{
	// ID 0 is reserved for Super read of all reads
	assert(ID > 0 && ID <= reads->size());
	return reads->at(ID - 1);
}

template <typename T>
void tDataSet<T>::sortReads()
{
	std::sort(reads->begin(), reads->end());
	is_sorted = true;
}

template <typename T>
void tDataSet<T>::rmDupReads()
{
	if(!is_sorted)
		sortReads();
	reads->erase(unique(reads->begin(), reads->end()), reads->end());
	has_duplicates = false;
	numberOfUniqueReads = reads->size();

	UINT64 i = 1;
	// After removing duplicate reads, assign IDs to them,
	// starting with 1. 0 is reserved for super read
	typename vector<T>::iterator it;
	for(it = reads->begin(); it!= reads->end(); ++it){
		it->setReadNumber(i);
		++i;
	}
}

template <typename T>
ostream& operator<< (ostream &out, tDataSet<T> &data)
{
	for(UINT64 j = 0; j < data.size(); j++){
		out << setw(7) << j+1 << ": ";
		out << setw(7) << data.getReadFromID(j+1).getReadNumber() << " # "
			<< setw(7) << data.getReadFromID(j+1).superID 
			<< " # " << data.getReadFromID(j+1) <<endl;
	}
	return out;
}
#endif /* TDATASET_H */
