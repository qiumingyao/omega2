/*
 * =====================================================================================
 *
 *       Filename:  tRead.h
 *
 *    Description:  template Read class
 *    		    typename T represent the type of elements that form a read
 *
 *        Version:  1.0
 *        Created:  07/23/2015 11:11:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  JJ Crosskey (cjg), jjchai01@gmail.com
 *   Organization:  ORNL
 *
 * =====================================================================================
 */
#ifndef TREAD_H_
#define TREAD_H_

#include "Config.h"
#include "namespace_jj.cpp"
template <typename T>
class tRead
{
	private:
		UINT64 readNumber;

		// Forward "string" of the read
		vector<T> read;

		// Reverse complement "string" of the read
		vector<T> readReverse;

	public:
		// ID of the read that contains this read
		UINT64 superID;
		
		tRead(void)
		{
			readNumber = 0;
		}

		tRead(const vector<T> & r)
		{
			readNumber = 0;
			superID = 0;
			read = r;
			for(INT64 i = read.size() - 1; i >= 0; i--){
				readReverse.push_back(jj::complement(read.at(i)));
			}
		}

		~tRead()
		{
//			cout << "tRead<T> deconstructor called on read " << readNumber << " at address " << this << endl;
		}

		void setReadNumber(UINT64 id)
		{
			assert(id > 0);
			readNumber = id;
		}

		T& operator[](size_t nIndex)
		{
			assert(nIndex < read.size());
			return read[nIndex];
		}

		const T& at (size_t nIndex) const
		{
			assert(nIndex < read.size());
			return read[nIndex];
		}

		// Get substring of the read or reverse read
		vector<T> substr(size_t pos, size_t len, char orient = 0) const;

		vector<T> getStringForward(void) const {return read;} 

		vector<T> getStringReverse(void) const {return readReverse;} 

		UINT64 getReadLength(void) const {return read.size();} 

		UINT64 size(void) const {return read.size();} 

		UINT64 getReadNumber(void) const {return readNumber;} 

		bool operator==(const tRead<T> &rhs) const { return(read == rhs.read);}
		
		bool operator<(const tRead<T> &rhs) const { return ( read < rhs.read); }
};

template <class T>
vector<T> tRead<T>::substr(size_t pos, size_t len, char orient) const
{
//	cerr << "Substring of read " << *this << " at pos " << pos 
//		<< " with length " << len << " and orientation " << static_cast<int>(orient) <<  endl;
	assert(len > 0 && (pos+len) <= read.size());

	vector<T> sub_read;
	// orient = 0: forward string of the read
	// orient = 1: reverse string of the read
	if (orient == 0){
		for(size_t i = 0; i < len; i++){
			sub_read.push_back(read.at(pos + i));
		}
	}
	else{
		for(size_t i = 0; i < len; i++){
			sub_read.push_back(readReverse.at(pos + i));
		}
	}
	return sub_read;
}


template <class T>
bool operator==(vector<T> &lhs, vector<T> &rhs)
{
	if(lhs.size() != rhs.size())
		return false;

	for(UINT64 i = 0; i < lhs.size(); ++i){
		if(lhs.at(i) != rhs.at(i))
			return false;
	}
	return true;
}

template <class T>
ostream & operator<<(ostream & out, const vector<T> & v)
{
	for(UINT64 i = 0; i < v.size(); i++)
		out << v.at(i) << "\t";
	return out;
}

template <class T>
bool operator<(tRead<T> &lhs, tRead<T> &rhs){
	return (lhs.getStringForward() < rhs.getStringForward());
}

template <class T>
bool tRead_comp(tRead<T> &lhs, tRead<T> &rhs){
	return (lhs.getStringForward() < rhs.getStringForward());
}

template <class T>
bool is_equal_tRead(tRead<T> * lhs, tRead<T> * rhs){
	return (lhs->getStringForward() == rhs->getStringForward());
}

template <class T>
ostream& operator<< (ostream &out, const tRead<T> &r)
{
	for(UINT64 i = 0; i < r.size(); i++){
		out << setw(8) << r.at(i) << "\t" ;
	}
	return out;
}

#endif // TREAD_H_
