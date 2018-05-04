/*****************************************************************************
 * 
 * File:    tbbhashmap.cpp
 * Author:  Alex Stivala,  Intel TBB code based on example in manual.
 * Created: May 2009
 *
 * C callable interface to Intel Thread Building Blocks conccurrent hashmap.
 * 
 * $Id: tbbhashmap.cpp 2404 2009-05-17 02:34:13Z astivala $
 *
 * Developed on:
 * Linux charikar.csse.unimelb.edu.au 2.6.9-78.0.5.ELsmp #1 SMP Wed Sep 24 05:40:24 EDT 2008 x86_64 x86_64 x86_64 GNU/Linux
 * TBB21_INSTALL_DIR=/home/charikar/pgrad/astivala/tbb21_20080605oss
 * . ~/tbb21_20080605oss/em64t/cc3.4.3_libc2.3.4_kernel2.6.9/bin/tbbvars.sh
 * TBB_ARCH_PLATFORM=em64t/cc3.4.3_libc2.3.4_kernel2.6.9
 *
 *
 * Preprocessor symbols:
 *
 * USE_GOOD_HASH  - use mixing hash function rather than trivial one
 *
 *****************************************************************************/

#include "tbb/concurrent_hash_map.h"
#include "tbb/scalable_allocator.h"
#include <pthread.h>
#include "tbbhashmap.h"

using namespace tbb;
using namespace std;

#define USE_GOOD_HASH


/***************************************************************************
 *
 * concurrent hash map classes
 *
 ***************************************************************************/


#ifdef USE_GOOD_HASH
/*
  hash a 64 bit value into 32 bits. From:
  (Thomas Wang, Jan 1997, Last update Mar 2007, Version 3.1)
  http://www.concentric.net/~Ttwang/tech/inthash.htm
  (found by reference in NIST Dictionary of Algorithms and Data Structures)
*/
static unsigned long hash6432shift(unsigned long long key)
{
  key = (~key) + (key << 18); /* key = (key << 18) - key - 1; */
  key = key ^ (key >> 31);
  key = key * 21; /* key = (key + (key << 2)) + (key << 4); */
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);
  return (unsigned long) key;
}

#endif



//! Structure that defines hashing and comparison operations for user's type.
struct MyHashCompare {
    static size_t hash( const _SET& x ) {
        size_t h = 0;
        unsigned int i;
        unsigned long highhash,lowhash,q;
        
#ifdef USE_GOOD_HASH
        highhash = hash6432shift(x.high);
        lowhash = hash6432shift(x.low);
        q = lowhash ^ highhash; /* FIXME: should have a hash12832shift() instead */
#else
        q = x.low;
#endif
        h = (size_t)q;
        return h;        
    }

    //! True if strings are equal
    static bool equal( const _SET& x, const _SET& y ) {
      return x.low == y.low && x.high == y.high;
    }
};

//! A concurrent hash table that maps _SETs to ints.
typedef concurrent_hash_map<_SET,int,MyHashCompare,
                            scalable_allocator<std::pair<_SET,int> > > SetHashTable;


/***************************************************************************
 *
 * global data
 *
 ***************************************************************************/


SetHashTable table; // The hash table



/***************************************************************************
 *
 * external functions
 *
 ***************************************************************************/

extern "C" void tbbhashmap_insert(_SET key, int value)
{
  SetHashTable::accessor a; // destroyed on exiting this scope
  table.insert(a, key);
  a->second = value;
}

extern "C" int tbbhashmap_lookup(_SET key)
{
  SetHashTable::const_accessor ca;
  if (table.find(ca, key))
    return ca->second;
  else
    return 0;
}



extern "C" int tbbhashmap_haskey(_SET key)
{
  SetHashTable::const_accessor ca; // destroyed on exiting this scope
  return table.find(ca, key) ? 1 : 0;
}

