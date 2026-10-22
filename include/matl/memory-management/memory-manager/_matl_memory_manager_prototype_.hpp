/**************************************************************************************************
 * File: _matl_memory_manager_prototype_.hpp
 *
 * Description: Prototype for the memory manager class.
 **************************************************************************************************/
namespace matl
{

   /***********************************************************************************************
    * class MemoryManager
    *
    * Per-thread memory manager (a thread_local singleton). It owns the thread's TemporaryArena,
    * from which the temporary matrices created during expression evaluation are allocated (e.g.
    * the result of A*B in 'D = A*B + C'), and holds the settings that control it.
    *
    * Matrices themselves allocate their memory on the heap (see allocate/deallocate) or use
    * externally owned memory.
    ***********************************************************************************************/
    class MemoryManager
    {
      private:
        // Memory manager enabled or disabled
        static inline thread_local bool ENABLED = defaults::ENABLE_MEMORY_MANAGER;

        // Stack-ordered memory for evaluation temporaries
        TemporaryArena temporaries_;

      protected:
        // Constructor
        MemoryManager() = default;

        // Destructor (releases the arena's memory when the thread exits)
        ~MemoryManager() = default;

      public:

        // Not cloneable
        MemoryManager(const MemoryManager&) = delete;

        // Not assignable
        MemoryManager& operator=(const MemoryManager&) = delete;

        // Get memory manager
        static MemoryManager& get_memory_manager();

        // Check if memory manager is enabled
        bool enabled() const;

        // Enable memory manager
        void enable();

        // Disable memory manager
        void disable();

        // Restore default settings
        void restore_defaults();

        // Memory for evaluation temporaries (this thread's arena)
        static TemporaryArena& temporary_arena();

        // Are evaluation temporaries allocated from the arena? (MemoryManager enabled and
        // CACHE_TEMPORARY_MATRICES set). Otherwise they are allocated on the heap.
        static bool use_temporary_arena();


        /**************** Settings *****************/

        // Allocate the temporary matrices created during expression evaluation from a recycled,
        // stack-ordered arena (see TemporaryArena) rather than the heap
        static inline thread_local bool CACHE_TEMPORARY_MATRICES = defaults::ENABLE_MEMORY_MANAGER;
    };
}
