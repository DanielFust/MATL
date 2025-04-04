#ifndef _MATL_BLOCKING__HPP
#define _MATL_BLOCKING__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains shared information for cache blocking for the matl library
 **************************************************************************************************/
//#include "_matl_forward_declarations_.hpp"

namespace matl
{
    // Anonymous namespace to protect accidental parameter modifications
    namespace 
    {
        // L1 Cache size of Intel 8-core i9 CPU is 32768 Bytes

        const static thread_local index_t DEFAULT_CHUNK_SIZE = 64;

        static thread_local index_t CHUNK_SIZE = 64;
    }    


   /***********************************************************************************************
    * void set_chunk_size(index_t)
    * 
    * Sets the number of bytes to be considered in a chunk of cache memory. 
    * A chunk_size x chunk_size block should fit within cache memory
    * 
    * Parameters:
    *     chunk_size        Number of bytes allowed in a cache block
    ***********************************************************************************************/
    inline void set_chunk_size(index_t N)
    {
        CHUNK_SIZE = N;
    }

   /***********************************************************************************************
    * void reset_chunk_size()
    * 
    * Sets the chunk size to its default value
    ***********************************************************************************************/
    inline void reset_chunk_size()
    {
        CHUNK_SIZE = DEFAULT_CHUNK_SIZE;
    }


   /***********************************************************************************************
    * index_t chunk_size()
    * 
    * Get the thread's chunk size
    ***********************************************************************************************/
    inline index_t chunk_size()
    {
        return CHUNK_SIZE;
    }

}


#endif





