#ifndef CONFIG_H
#define CONFIG_H

//============================================================================
// Name        : Config.h
// Author      : Tae-Hyuk (Ted) Ahn, JJ Crosskey
// Version     : v1.2
// Copyright   : 2015 Oak Ridge National Lab (ORNL). All rights reserved.
// Description : Config header file
//============================================================================

//multi-thread library OPENMP
//#include <omp.h>

// C headers:
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <assert.h> // for assert()

// C++ headers:
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>
#include <utility>
#include <limits>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <map>
#include <unordered_map>
#include <typeinfo>
#include <tuple>

// Local headers
#include "logcpp/log.h"
using namespace std;


#ifndef BUFFER_SIZE
#define BUFFER_SIZE 100
#endif

#define DEBUG


//============================================================================
// global variables with default value
//============================================================================
// minimum reads count in an edge to be not dead end edge (default: 10)
extern unsigned int minReadsCountInEdgeToBeNotDeadEnd;	

// minimum edge length to be not dead end edge (default: 500)
extern unsigned int minEdgeLengthToBeNotDeadEnd;		

// minimum reads count in an edge to be 1 minimum flow (default: 20)
extern unsigned int minReadsCountInEdgeToBe1MinFlow;	

// minimum edge length to be 1 minimum flow (default: 1000)
extern unsigned int minEdgeLengthToBe1MinFlow;			

// minimum contig length to be reported (default: 500)
extern unsigned int minContigLengthTobeReported;			

// -mdc: Minimum overlap length difference to clip branches (default: 50)
extern unsigned int minOvlDiffToClip;	

// -mfs: Minimum fold difference to consider branches to be short (default: 5)
extern unsigned int minFoldToBeShortBranch;	

// report unoverlapped reads to contigs (default: false)
extern bool reportUnoverlappedReadsToContigs;			
//============================================================================


//============================================================================
// variables typedef
//============================================================================
typedef unsigned char UINT8;
typedef signed char INT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef unsigned long UINT32;
typedef long INT32;
typedef unsigned long long UINT64;
typedef long long INT64;


//============================================================================
//	Exit code that displays the place of exit and message.
//============================================================================
#define MYEXIT(a) { cerr << endl << "Exit from File: " << __FILE__ << " Line: " << __LINE__ << " Function: " << __FUNCTION__ << "()" << endl << "Message: " << a << endl; exit(0);}


//============================================================================
// Clock logging for debugging
//============================================================================
#define CLOCKSTART clock_t begin = clock(); INT64 mem_start = checkMemoryUsage(); \
		FILE_LOG(logDEBUG) << endl << ">>> Function start: "<< __FUNCTION__ << "()" << endl;
#define CLOCKSTOP clock_t end = clock();  INT64 mem_end = checkMemoryUsage(); \
		FILE_LOG(logINFO) << "<<< Function stop: " << __FUNCTION__ << "(), Elapsed time: " << double(end - begin) / CLOCKS_PER_SEC<< " seconds, Memory usage: " \
		<< mem_end << " - " <<  mem_start << " = "<< mem_end - mem_start << " MB."<< endl << "----" << endl;
#define PRINTMEM FILE_LOG(logDEBUG) << checkMemoryUsage() << " MB used" << endl;

//#define CLOCKSTART double begin = omp_get_wtime(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
//#define CLOCKSTOP double end = omp_get_wtime(); cout << "Function " << __FUNCTION__ << "() finished in " << double(end - begin) << " Seconds." << endl << endl;

// Get the memory usage with a Linux kernel.
inline unsigned int checkMemoryUsage()
{
    // get KB memory into count
    unsigned int count=0;

    #if defined(__linux__)
    ifstream f("/proc/self/status"); // read the linux file
    while(!f.eof()){
        string key;
        f>>key;
        if(key=="VmData:"){     // size of data
            f>>count;
        break;
        }
    }
    f.close();
    #endif

    // return MBs memory (size of data)
    return (count/1024);
}


//============================================================================
// Class of Config
//=============================================================================
class Config
{

public:
    // Default constructor
    Config();

    // Default destructor
    ~Config();

    // Variables
    static vector<string> readFilenamesList;
    static vector<string> edgeFilenamesList;
    static string outFilenamePrefix;
    static UINT64 minOvl;
    static UINT64 maxSubs;
    static UINT64 maxEdits;
    static string pacbioSamFile;

    // Get options
    static bool setConfig(int argc, char **argv);

    // Print help
    static void printHelp();

    // Get read filenames
    static vector<string> getReadFilenames() {return readFilenamesList;}

    // Get edge filenames
    static vector<string> getEdgeFilenames() {return edgeFilenamesList;}

    // Get output prefix name
    static string getOutputFilenamePrefix() {return outFilenamePrefix;}

    // Get minimum overlap length
    static UINT64 getMinOvl() {return minOvl;}

    // Get maximum allowed substitutions
    static UINT64 getMaxEdits() {return maxEdits;}

    // Get maximum allowed substitutions
    static UINT64 getMaxSubs() {return maxSubs;}

    // Get verbosity level of log messages
    static std::string getLogLevel(){return FILELog::ToString(FILELog::ReportingLevel());}
};


#endif
