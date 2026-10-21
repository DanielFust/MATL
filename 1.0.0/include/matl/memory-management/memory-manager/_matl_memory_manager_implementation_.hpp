/**************************************************************************************************
 * File: _matl_memory_manager_implementation_.hpp
 *
 * Description: Implementation for the memory manager class.
 **************************************************************************************************/
namespace matl
{
   /***********************************************************************************************
    * MemoryManager& get_memory_manager()
    *
    * Many of the settings and methods of the MemoryManager are static to avoid the syntax
    * MemoryManager::get_memory_manager().method() which is too cumbersome for the purpose of MATL
    * which is to provide a simple and intuitive interface for matrix operations. Static settings
    * are initialized from the defaults in matl.hpp; users may change them before or after the
    * singleton is created.
    *
    * Returns:
    *    Reference to the memory manager singleton for this thread (created on first use and
    *    destroyed when the thread exits)
    ***********************************************************************************************/
    inline MemoryManager& MemoryManager::get_memory_manager()
    {
        static thread_local MemoryManager memory_manager;
        return memory_manager;
    }


   /***********************************************************************************************
    * bool enabled()
    *
    * Check if memory manager is enabled.
    *
    * Returns:
    *     True if memory manager is enabled, false otherwise
    ***********************************************************************************************/
    inline bool MemoryManager::enabled() const
    {
        return MemoryManager::ENABLED;
    }


   /***********************************************************************************************
    * void enable()
    *
    * Enable memory manager.
    ***********************************************************************************************/
    inline void MemoryManager::enable()
    {
        MemoryManager::ENABLED = true;
    }

   /***********************************************************************************************
    * void disable()
    *
    * Disable memory manager. Evaluation temporaries are then allocated on the heap.
    ***********************************************************************************************/
    inline void MemoryManager::disable()
    {
        MemoryManager::ENABLED = false;
    }


   /***********************************************************************************************
    * void restore_defaults()
    *
    * Restore default settings for memory manager.
    ***********************************************************************************************/
    inline void MemoryManager::restore_defaults()
    {
        MemoryManager::ENABLED                  = defaults::ENABLE_MEMORY_MANAGER;
        MemoryManager::CACHE_TEMPORARY_MATRICES = defaults::ENABLE_MEMORY_MANAGER;
    }


   /***********************************************************************************************
    * TemporaryArena& temporary_arena()
    *
    * Returns:
    *     This thread's arena for evaluation temporaries
    ***********************************************************************************************/
    inline TemporaryArena& MemoryManager::temporary_arena()
    {return get_memory_manager().temporaries_;}

   /***********************************************************************************************
    * bool use_temporary_arena()
    *
    * Returns:
    *     'true' if evaluation temporaries are allocated from the arena (the MemoryManager is
    *     enabled and CACHE_TEMPORARY_MATRICES is set), 'false' if they are allocated on the heap
    ***********************************************************************************************/
    inline bool MemoryManager::use_temporary_arena()
    {return MemoryManager::ENABLED && MemoryManager::CACHE_TEMPORARY_MATRICES;}
}
