/**************************************************************************************************
 * File: _matl_temporary_arena_implementation_.hpp
 *
 * Description: Implementation of the TemporaryArena class.
 **************************************************************************************************/
#ifndef _MATL_TEMPORARY_ARENA_IMPLEMENTATION__HPP
#define _MATL_TEMPORARY_ARENA_IMPLEMENTATION__HPP
namespace matl
{
   /***********************************************************************************************
    * Constructor
    *
    * No memory is allocated for the arena until it is first used.
    ***********************************************************************************************/
    inline TemporaryArena::TemporaryArena() : current_(0), bytes_in_use_(0), high_water_(0)
    {}

   /***********************************************************************************************
    * void add_chunk(index_t bytes)
    *
    * Appends a new chunk of at least 'bytes' bytes. Chunks grow geometrically (each new chunk is
    * at least as large as the current capacity) so the number of chunks stays small.
    ***********************************************************************************************/
    inline void TemporaryArena::add_chunk(index_t bytes)
    {
        index_t size = bytes;
        if (size < MIN_CHUNK_SIZE)   size = MIN_CHUNK_SIZE;
        if (size < this->capacity()) size = this->capacity();
        this->chunks_.push_back(Chunk{std::make_unique<std::byte[]>(size), size, 0});
    }

   /***********************************************************************************************
    * void consolidate()
    *
    * Replaces several chunks with a single chunk of their total size. Only valid while no blocks
    * are allocated.
    ***********************************************************************************************/
    inline void TemporaryArena::consolidate()
    {
        if (this->chunks_.size() <= 1) return;
        const index_t total = this->capacity();
        this->chunks_.clear();
        this->chunks_.push_back(Chunk{std::make_unique<std::byte[]>(total), total, 0});
        this->current_ = 0;
    }

   /***********************************************************************************************
    * std::byte* allocate(index_t bytes)
    *
    * Allocates a block at the top of the stack.
    *
    * Parameters:
    *     bytes        Requested size in bytes
    * Returns:
    *     Pointer to a block of at least 'bytes' bytes, aligned to ALIGNMENT
    ***********************************************************************************************/
    inline std::byte* TemporaryArena::allocate(index_t bytes)
    {
        // Round up to the alignment (every block has a distinct address, even for 0 bytes)
        const index_t n = ((bytes + ALIGNMENT - 1) / ALIGNMENT) * ALIGNMENT + (bytes == 0 ? ALIGNMENT : 0);

        // Find a chunk, starting at the top of the stack, with room for the block. Chunks after
        // the current one hold no blocks.
        index_t chunk = this->current_;
        while (chunk < this->chunks_.size() && this->chunks_[chunk].size - this->chunks_[chunk].used < n)
            ++chunk;
        if (chunk == this->chunks_.size())
            this->add_chunk(n);

        Chunk& c = this->chunks_[chunk];
        std::byte* ptr = c.data.get() + c.used;
        this->records_.push_back(Record{ptr, chunk, c.used, n, false});
        c.used         += n;
        this->current_  = chunk;

        this->bytes_in_use_ += n;
        if (this->bytes_in_use_ > this->high_water_) this->high_water_ = this->bytes_in_use_;
        return ptr;
    }

   /***********************************************************************************************
    * void release(std::byte* ptr)
    *
    * Releases a block. Blocks at the top of the stack are reclaimed immediately; a block released
    * out of order is reclaimed once all blocks allocated after it have been released.
    *
    * Parameters:
    *     ptr          Pointer returned by allocate()
    ***********************************************************************************************/
    inline void TemporaryArena::release(std::byte* ptr)
    {
        // Find the block (normally the top of the stack)
        index_t k = this->records_.size();
        while (k > 0 && this->records_[k-1].ptr != ptr) --k;
        if (k == 0 || this->records_[k-1].released)
            throw std::logic_error("TemporaryArena::release(std::byte*): pointer is not an allocated block");

        Record& record = this->records_[k-1];
        record.released      = true;
        this->bytes_in_use_ -= record.bytes;

        // Reclaim released blocks from the top of the stack
        while (!this->records_.empty() && this->records_.back().released)
        {
            const Record& top = this->records_.back();
            this->chunks_[top.chunk].used = top.offset_before;
            this->current_ = top.chunk;
            this->records_.pop_back();
        }
        if (!this->records_.empty()) this->current_ = this->records_.back().chunk;

        // When empty, merge chunks so the same work fits in one chunk next time
        if (this->records_.empty()) this->consolidate();
    }

   /***********************************************************************************************
    * void reserve(index_t bytes)
    *
    * Ensures that at least 'bytes' bytes can be allocated without heap allocation. When the arena
    * is empty its memory is consolidated into a single chunk of at least 'bytes' bytes.
    *
    * Parameters:
    *     bytes        Number of bytes to make available
    ***********************************************************************************************/
    inline void TemporaryArena::reserve(index_t bytes)
    {
        if (this->records_.empty())
        {
            if (this->capacity() >= bytes && this->chunks_.size() <= 1) return;
            const index_t size = (bytes > this->capacity()) ? bytes : this->capacity();
            this->chunks_.clear();
            this->chunks_.push_back(Chunk{std::make_unique<std::byte[]>(size), size, 0});
            this->current_ = 0;
        }
        else
        {
            // Blocks are in use: make sure a chunk after the top has room for 'bytes'
            for (index_t k = this->current_; k < this->chunks_.size(); ++k)
                if (this->chunks_[k].size - this->chunks_[k].used >= bytes) return;
            this->add_chunk(bytes);
        }
    }

   /***********************************************************************************************
    * void release_memory()
    *
    * Frees all memory held by the arena. The arena must be empty.
    ***********************************************************************************************/
    inline void TemporaryArena::release_memory()
    {
        if (!this->records_.empty())
            throw std::logic_error("TemporaryArena::release_memory(): temporaries are still in use");
        this->chunks_.clear();
        this->chunks_.shrink_to_fit();
        this->current_ = 0;
    }

   /***********************************************************************************************
    * Statistics
    ***********************************************************************************************/
    inline index_t TemporaryArena::bytes_in_use() const
    {return this->bytes_in_use_;}

    inline index_t TemporaryArena::capacity() const
    {
        index_t total = 0;
        for (const Chunk& c : this->chunks_) total += c.size;
        return total;
    }

    inline index_t TemporaryArena::high_water_mark() const
    {return this->high_water_;}

    inline void TemporaryArena::reset_high_water_mark()
    {this->high_water_ = this->bytes_in_use_;}

    inline index_t TemporaryArena::blocks_in_use() const
    {return this->records_.size();}
}
#endif
