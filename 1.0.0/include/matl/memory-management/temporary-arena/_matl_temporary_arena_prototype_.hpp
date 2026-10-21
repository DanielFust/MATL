/**************************************************************************************************
 * File: _matl_temporary_arena_prototype_.hpp
 *
 * Description: Prototype for the TemporaryArena class.
 **************************************************************************************************/
#ifndef _MATL_TEMPORARY_ARENA_PROTOTYPE__HPP
#define _MATL_TEMPORARY_ARENA_PROTOTYPE__HPP
namespace matl
{
   /***********************************************************************************************
    * class TemporaryArena
    *
    * Stack-ordered memory for the temporary matrices created while an expression is evaluated
    * (e.g. the result of A*B in 'D = A*B + C'). Each thread's MemoryManager owns one arena.
    *
    * Temporaries live for exactly one assignment and are released in the reverse order of their
    * allocation, so the arena allocates by advancing an offset through a chunk of memory and
    * releases by moving it back. Once the arena has grown to fit the largest set of temporaries
    * used at one time, evaluation performs no heap allocation.
    *
    *   - Growth: if a request does not fit, a new chunk is added (heap allocation). When the arena
    *     next becomes empty, its chunks are merged into a single chunk of their total size, so
    *     the same work then fits without growing.
    *
    *   - Out-of-order release: a block released while blocks allocated after it are still in use
    *     is marked free, and its memory is reclaimed once every block above it has been released.
    *     Evaluation always releases in stack order; this keeps out-of-order release safe.
    *
    *   - reserve() preallocates capacity up front, so that even the first evaluation needs no
    *     heap allocation.
    *
    * Blocks are aligned to alignof(std::max_align_t).
    ***********************************************************************************************/
    class TemporaryArena
    {
      private:
        // A contiguous chunk of memory
        struct Chunk
        {
            std::unique_ptr<std::byte[]> data;
            index_t                      size;
            index_t                      used;
        };

        // One allocation, in stack order
        struct Record
        {
            std::byte* ptr;
            index_t    chunk;          // index of the chunk holding the block
            index_t    offset_before;  // chunk's 'used' offset before the block was allocated
            index_t    bytes;          // size of the block (after rounding for alignment)
            bool       released;       // released, but not yet reclaimed (out-of-order release)
        };

        std::vector<Chunk>  chunks_;
        std::vector<Record> records_;
        index_t             current_;      // index of the chunk at the top of the stack
        index_t             bytes_in_use_;
        index_t             high_water_;

        // Allocate a new chunk at the end of the chunk list
        void add_chunk(index_t bytes);

        // Merge all chunks into one (only when the arena is empty)
        void consolidate();

      public:
        // Alignment of every block
        static constexpr index_t ALIGNMENT = alignof(std::max_align_t);

        // Smallest chunk the arena allocates
        static constexpr index_t MIN_CHUNK_SIZE = 64 * 1024;

        // Constructor (allocates no memory until first use)
        TemporaryArena();

        // Not copyable (blocks are referenced by address)
        TemporaryArena(TemporaryArena const&) = delete;
        TemporaryArena& operator=(TemporaryArena const&) = delete;

        // Allocate a block of at least 'bytes' bytes
        std::byte* allocate(index_t bytes);

        // Release a block returned by allocate()
        void release(std::byte* ptr);

        // Ensure at least 'bytes' bytes can be allocated without heap allocation
        void reserve(index_t bytes);

        // Free all memory held by the arena (the arena must be empty)
        void release_memory();

        // Bytes currently allocated to temporaries
        index_t bytes_in_use() const;

        // Total bytes held by the arena
        index_t capacity() const;

        // Largest number of bytes in use at one time since the last reset
        index_t high_water_mark() const;

        // Reset the high water mark to the current usage
        void reset_high_water_mark();

        // Number of blocks currently allocated (including released blocks not yet reclaimed)
        index_t blocks_in_use() const;
    };
}
#endif
