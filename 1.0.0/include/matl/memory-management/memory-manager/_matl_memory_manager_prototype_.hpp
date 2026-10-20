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
    * Singleton memory manager class for managing memory allocation and deallocation for the 
    * matrices class. The memory manager is responsible for assigning preallocated scrath memory
    * to matrices and expressions and reclaiming said scratch memory when objects are destroyed.
    ***********************************************************************************************/
    class MemoryManager
    {
      private:
        // Singleton instance
        static inline thread_local MemoryManager* memory_manager = nullptr;

        // Memory manager enabled or disabled
        static inline thread_local bool ENABLED = defaults::ENABLE_MEMORY_MANAGER;

        // Scratch memory containers. A std::deque is used (rather than std::vector) because
        // ScratchContainers and ScratchMemoryBorrowers hold raw pointers to the containers.
        // std::deque::emplace_back does not invalidate references to existing elements.
        std::deque<ScratchContainer> scratch_;

      protected:  
        // Constructor
        MemoryManager();

        // Destructor
        ~MemoryManager();

        // Split scratch array into two arrays
        void split(ScratchContainer& scratch_array, index_t size);
    
      public:
        
        // Not cloneable
        MemoryManager(const MemoryManager&) = delete;

        // Not assignable
        MemoryManager& operator=(const MemoryManager&) = delete;

        // Get memory manager
        static MemoryManager& get_memory_manager();

        // Get scratch memory containers
        std::deque<ScratchContainer>& containers();

        // Merge unused contiguous scratch arrays
        void merge_contiguous_unused();

        // Allocate scratch memory
        ScratchContainer& allocate_scratch(index_t size);

        // Deallocate scratch memory
        //void deallocate();

        // Add preallocated scratch memory and return its container
        template<class T>
        ScratchContainer& add_scratch_memory(T* scratch, index_t size);

        // Remove preallocated scratch memory
        template<class T>
        void remove_scratch_memory(T* scratch);

        // Reserve scratch memory container
        ScratchContainer& reserve_scratch_memory(index_t size);

        // Allocate RowMajorMatrix in scratch memory
        template<class T>
        void assign_memory(RowMajorMatrix<T>& matrix, index_t rows, index_t cols);

        // Free scratch memory
        void free_memory(ScratchContainer& scratch_array);

        // Delete all scratch memory
        void delete_scratch_memory();

        // Assign borrower to scratch memory container
        void assign_borrower(ScratchContainer*, ScratchMemoryBorrower*);

        // Check if memory manager is enabled
        bool enabled() const;

        // Enable memory manager
        void enable();

        // Disable memory manager
        void disable();

        // Restore default settings
        void restore_defaults();


       /*******************************************************************************************
        * Procedures for allocating and deallocating matrices are defined as static member 
        * functions of the MemoryManager class, as it has 'friend' access to the private data
        * of the matrix classes.
        * 
        * For logical and concise syntax, 'allocate' and 'deallocate' template functions are
        * defined in the 'matl' namespace that map to the static member functions of the
        * MemoryManager class. This is intended to loosely resemble Fortran90 allocate/deallocate
        * vernacular and provide greater clarity to the user when the MemoryManager is disabled
        * or ALLOCATE_MATRICES_IN_SCRATCH is set to 'false'. 
        *******************************************************************************************/

        // RowMajorMatrix (Allocate)
        template<class T>
        static void allocate(RowMajorMatrix<T>&, index_t, index_t);

        // RowMajorMatrix (Deallocate)
        template<class T>
        static void deallocate(RowMajorMatrix<T>&);

        // FixedSizeRowMajorMatrix (checks dimensions and returns)
        template<class T, index_t M, index_t N>
        static void allocate(FixedSizeRowMajorMatrix<T, M, N>&, index_t, index_t);

        // FixedSizeRowMajorMatrix (Deallocate)
        //template<class T, index_t M, index_t N>
        //static void deallocate(FixedSizeRowMajorMatrix<T, M, N>&);


        


        /**************** Settings *****************/

        // Cache temporary matrices
        static inline thread_local bool    CACHE_TEMPORARY_MATRICES     = defaults::ENABLE_MEMORY_MANAGER;

        // Allocate matrices in scratch memory
        static inline thread_local bool    ALLOCATE_MATRICES_IN_SCRATCH = defaults::ALLOCATE_MATRICES_IN_SCRATCH;

        // The MemoryManager will allocate scratch memory in multiples of this size
        static inline thread_local index_t SCRATCH_BLOCK_SIZE           = defaults::SCRATCH_BLOCK_SIZE;
    }; 



   /***********************************************************************************************
    * Allocating and Deallocating Matrices
    * 
    * 'allocate' and 'deallocate' functions in the 'matl' namespace are the preferred general 
    * purpose procedures for allocation/deallocation of matrices. If the MemoryManager is enabled,
    * and 'ALLOCATE_MATRICES_IN_SCRATCH' is set to 'true' assignment of memory to matrices
    * is deferred to the MemoryManager. Otherwise, the memory is allocated on the heap and owned
    * by the matrix objects themselves.
    * 
    * Note: These functions map to the static member functions of the MemoryManager class which has
    *       'friend' access to the private data of the matrix classes, however the allocation
    *       and deallocation of memory is not universally the conceptual responsibility of the 
    *       MemoryManager. When the MemoryManager is disabled or 'ALLOCATE_MATRICES_IN_SCRATCH' is 
    *       set to 'false' the user-facing syntax should not falsely suggest that the MemoryManager 
    *       is using cached memory for the matrices.
    ***********************************************************************************************/

    // Matrix allocation
    template<class T>
    void allocate(RowMajorMatrix<T>&, index_t, index_t);

    // Matrix deallocation
    template<class T>
    void deallocate(RowMajorMatrix<T>&);

}
