/**************************************************************************************************
 * File: _matl_memory_manager_implementation_.hpp
 * 
 * Description: Implementation for the memory manager class.
 **************************************************************************************************/
namespace matl
{
   /***********************************************************************************************
    * Constructor
    * 
    * Non-static defaults should be set here. Static defaults are set externally in matl.hpp
    * Users may choose to alter the settings of the MemoryManager before the singleton is 
    * initialized and the constructor must not override these settings.
    * 
    * Many of the settings and methods of the MemoryManager are static to avoid the syntax
    * MemoryManager::get_memory_manager().method() which is too cumbersome for the purpose of MATL
    * which is to provide a simple and intuitive interface for matrix operations.
    ***********************************************************************************************/
    inline MemoryManager::MemoryManager() {}

   /***********************************************************************************************
    * Destructor
    ***********************************************************************************************/
    inline MemoryManager::~MemoryManager()
    {
        // Deallocate scratch memory (externally owned memory is not deleted)
        get_memory_manager().scratch_.clear();
    } 
    

   /***********************************************************************************************
    * MemoryManager& get_memory_manager()
    * 
    * Returns:
    *    Reference to the memory manager singleton for this thread
    ***********************************************************************************************/
    inline MemoryManager& MemoryManager::get_memory_manager()
    {
        if (memory_manager == nullptr)
            memory_manager = new MemoryManager();
        return *memory_manager;
    } 

   /***********************************************************************************************
    * std::deque<ScratchContainer>& containers()
    * 
    * Returns:
    *   Reference to the scratch memory containers
    ***********************************************************************************************/
    inline std::deque<ScratchContainer>& MemoryManager::containers()
    {return get_memory_manager().scratch_;}  

   /***********************************************************************************************
    * void split(ScratchContainer&, index_t) 
    * 
    * Splits scratch array into two arrays so that excess memory can be used for other purposes.
    * 
    * Parameters:
    *     scratch_array        ScratchContainer object to split
    *     size                 Size of the first scratch array after split
    ***********************************************************************************************/
    inline void MemoryManager::split(ScratchContainer& scratch, index_t size)
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::split..." << std::endl;
        #endif

        index_t           excess      = scratch.size() - size;
        ScratchContainer* new_scratch = nullptr;

        #ifdef MATL_DEBUG
        // Check if split is possible
        if (excess <= 0) throw std::runtime_error("Cannot split scratch array. Size is too small.");
        #endif

        // If there is an unused next neighbor, donate excess memory
        if (scratch.next_neighbor() != nullptr)
        {
            ScratchContainer& next = *scratch.next_neighbor();
            if (next.bytes_used()==0)
            {
                // Donate excess memory to next neighbor
                next.size_ += excess;
                next.data_ =  &scratch.data_[size];
                next.prev_neighbor_ = &scratch;

                // Update scratch array
                scratch.size_          = size;
                scratch.next_neighbor_ = &next;

                return;
            }
        }

        // Check if an unallocated scratch array is available
        for (ScratchContainer& scratch : get_memory_manager().containers())
        {
            if (scratch.size() == 0)
            {
                scratch.nullify();
                new_scratch = &scratch;
                break;
            }
        }
        // If no unallocated scratch array is available, create a new one
        if (new_scratch == nullptr)
        {
            // Create new scratch array from excess memory
            std::byte* excess_data = &scratch.data_[size];
            get_memory_manager().containers().emplace_back(ScratchContainer(excess_data, excess));
            new_scratch = &get_memory_manager().scratch_.back();
        }

        // If next neighbor exists, update it
        if (scratch.next_neighbor_ != nullptr)
        {
            scratch.next_neighbor_->prev_neighbor_ = new_scratch;
        }

        // Update scratch array
        scratch.size_          = size;
        scratch.next_neighbor_ = new_scratch;

        // Update new scratch array (data must be set here, as a reused empty container has none)
        new_scratch->data_          = &scratch.data_[size];
        new_scratch->prev_neighbor_ = &scratch;
        new_scratch->external_mem_  = scratch.external_mem();
        new_scratch->size_          = excess;
        new_scratch->bytes_used_    = 0;
    }



   /***********************************************************************************************
    * void merge_contiguous_unused()
    * 
    * Combines all scratch arrays which are contiguous in memory and unused. This method should
    * not deallocate memory or modify the number of ScratchContainer objects owned by the memory 
    * manager as this could invalidate iterators.
    ***********************************************************************************************/
    inline void MemoryManager::merge_contiguous_unused()
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::merge_contiguous_unused..." << std::endl;
        #endif
        // Iterate through scratch arrays
        for (ScratchContainer& scratch : get_memory_manager().scratch_)
        {
            // Check if scratch array is unused
            if ((scratch.size() > 0)  && (scratch.bytes_used() == 0))
                scratch.merge_contiguous_unused();
        }
    }
    
    
   /***********************************************************************************************
    * void allocate(index_t)
    * 
    * Allocates scratch memory for the memory manager.
    * 
    * Parameters:
    *    size        Size (in bytes) of scratch memory to allocate
    ***********************************************************************************************/
    inline ScratchContainer& MemoryManager::allocate_scratch(index_t size)
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::allocate_scratch..." << std::endl;
        #endif
        // First check if any containers are empty
        for (ScratchContainer& container : get_memory_manager().scratch_)
        {
            if (container.size() == 0)
            {
                // Allocate memory
                container.nullify();
                container.data_ = new std::byte[size];
                container.size_ = size;
                container.external_mem_ = false;
        
                return container;
            }
        }

        // If no empty containers are found, create a new one
        get_memory_manager().scratch_.emplace_back(ScratchContainer(size));    
        return get_memory_manager().scratch_.back();
    }
    
   /***********************************************************************************************
    *  void add_scratch_memory(T*, index_t)
    * 
    * Add preallocated scratch memory to the memory manager.
    * 
    * Parameters:
    *     scratch        Pointer to preallocated memory
    *     size           Size (number of elements) of preallocated memory
    * Returns:
    *     Reference to the scratch memory container
    ***********************************************************************************************/
    template<class T>
    ScratchContainer& MemoryManager::add_scratch_memory(T* scratch, index_t size)
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::add_scratch_memory..." << std::endl;
        #endif
        index_t nbytes = size * sizeof(T);
        // First check if there are any empty containers
        for (ScratchContainer& container : get_memory_manager().scratch_)
        {
            if (container.size() == 0)
            {
                // Add scratch memory (externally owned, so it is never deleted by the MemoryManager)
                container.nullify();
                container.data_         = reinterpret_cast<std::byte*>(scratch);
                container.size_         = nbytes;
                container.external_mem_ = true;
                return container;
            }
        }
        // If there are no empty containers to hold the scratch memory, create one
        get_memory_manager().scratch_.emplace_back(ScratchContainer(reinterpret_cast<std::byte*>(scratch), nbytes));
        return get_memory_manager().scratch_.back();
    } 
    
   /***********************************************************************************************
    * void remove_scratch_memory(T*)
    * 
    * Remove preallocated scratch memory from the memory manager.
    * 
    * Parameters:
    *     scratch        Pointer to preallocated memory
    ***********************************************************************************************/
    template<class T>
    void MemoryManager::remove_scratch_memory(T* scratch)
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::remove_scratch_memory..." << std::endl;
        #endif

        std::byte* address = reinterpret_cast<std::byte*>(scratch);

        // Find scratch array by address
        for (ScratchContainer& container : get_memory_manager().scratch_)
        {
            if (container.data() == address)
            {
                // For safety, only remove the scratch memory if the pointer addresses
                // the zeroth element of the raw byte array
                if (container.prev_neighbor() != nullptr) 
                {
                    throw std::runtime_error("Cannot remove scratch memory. Pointer does not point to the start of the allocated memory block.");
                    return;
                }

                // Iterate through contiguous scratch arrays and nullify their borrowers
                for (ScratchContainer* contig = &container; contig != nullptr; contig = contig->next_neighbor())
                {
                    if (contig->borrower() != nullptr)
                    {
                        contig->borrower()->nullify_borrowed_memory();
                        contig->borrower_ = nullptr;
                    }
                }

                // Nullify all containers using this scratch memory
                container.nullify_neighbors();

                // If memory is owned by the memory manager, deallocate it
                if (!container.external_mem()) delete[] container.data();

                // Nullify scratch array
                container.nullify();

                return;
            }
        }
    }


   /***********************************************************************************************
    * ScratchContainer& reserve_scratch_memory(index_t)
    * 
    * Reserve scratch memory for the memory manager.
    * 
    * Parameters:
    *     size        Size (in bytes) of scratch memory to reserve 
    * Returns:
    *     Reference to the reserved scratch memory container
    ***********************************************************************************************/
    inline ScratchContainer& MemoryManager::reserve_scratch_memory(index_t size)
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::reserve_scratch_memory..." << std::endl;
        #endif

        // New allocation size (if necessary)
        index_t nbytes;
        index_t block_size = MemoryManager::SCRATCH_BLOCK_SIZE;

        // Check if any containers are large enough
        #ifdef MATL_DEBUG
        std::cout << "Requested size: " << size << std::endl;
        #endif
        for (ScratchContainer& container : get_memory_manager().containers())
        {
            #ifdef MATL_DEBUG
            std::cout << "Container size: " << container.size() << std::endl;
            #endif
            if (container.size()-container.bytes_used() >= size)
            {
                #ifdef MATL_DEBUG
                std::cout << "Viable container found! bytes used: " << container.bytes_used() << std::endl;
                #endif
                // Split container if necessary
                if (container.bytes_used() != 0)
                {
                    //get_memory_manager().split(container, size);
                    get_memory_manager().split(container, container.bytes_used());
                    return *container.next_neighbor();
                }
                return container;
            }
        }

        // Allocate new scratch array
        nbytes = (1 + size/block_size)*block_size;
        get_memory_manager().scratch_.emplace_back(ScratchContainer(nbytes));
    
        return get_memory_manager().scratch_.back();
    } 


   /***********************************************************************************************
    * void assign_memory(RowMajorMatrix<T>&, index_t, index_t)
    * 
    * Assign scratch memory to a RowMajorMatrix object.
    * 
    * Parameters:
    *     matrix        RowMajorMatrix object to assign memory to 
    *     rows          Number of rows in matrix
    *     cols          Number of columns in matrix
    ***********************************************************************************************/
    template<class T>
    void MemoryManager::assign_memory(RowMajorMatrix<T>& matrix, index_t rows, index_t cols)
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::assign_memory..." << std::endl;
        #endif
        MemoryManager& manager = MemoryManager::get_memory_manager();
        
        // Calculate size of matrix
        index_t nbytes = rows * cols * sizeof(T);

        // Get scratch memory of appropriate size
        ScratchContainer& scratch = manager.reserve_scratch_memory(nbytes);

        // Set borrower and bytes used
        manager.assign_borrower(&scratch, &matrix);
        scratch.bytes_used_ = nbytes;
        
        // Assign memory to matrix
        matrix.data_    = reinterpret_cast<T*>(scratch.data());
        matrix.rows_    = rows;
        matrix.cols_    = cols;
        matrix.external_mem_ = false; // From perspective of matrix, memory is owned by the MemoryManager even if it was preallocated
    } 


   /***********************************************************************************************
    * void free_memory(ScratchContainer&)
    * 
    * Free scratch memory from the memory manager. Memory will be merged with contiguous unused
    * scratch arrays. This method should not deallocate memory or modify the number of ScratchContainer
    * objects owned by the memory manager as this could invalidate iterators.
    * 
    * Parameters:
    *     scratch_array        ScratchContainer object to free memory from
    ***********************************************************************************************/
    inline void MemoryManager::free_memory(ScratchContainer& scratch)
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::free_memory..." << std::endl;
        #endif
        // Set bytes used to zero
        scratch.bytes_used_ = 0;
        
        // Nullify the borrower so memory is no longer accessed
        if (scratch.borrower() != nullptr)
        {
            scratch.borrower()->nullify_borrowed_memory();
            scratch.borrower_           = nullptr;
        }
        
        // Merge contiguous unused scratch arrays
        scratch.merge_contiguous_unused();
    } 

   /***********************************************************************************************
    * void delete_scratch_memory()
    * 
    * Delete all scratch memory from the memory manager. If any memory is being borrowed, it will
    * be copied to privately owned memory.
    ***********************************************************************************************/
    inline void MemoryManager::delete_scratch_memory()
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::delete_scratch_memory..." << std::endl;
        #endif
        // First nullify borrowed memory so it can no longer be accessed then and free scratch 
        // memory
        for (ScratchContainer& scratch : get_memory_manager().scratch_)
        {
            // Copy memory to privately owned memory
            if (scratch.borrower() != nullptr) 
            {
                scratch.borrower()->nullify_borrowed_memory();
                scratch.borrower_ = nullptr;
            }

            // Free memory and merge arrays
            get_memory_manager().free_memory(scratch);
        }

        // All contiguous scratch arrays should be merged, so scratch arrays owned by the memory
        // manager can be safely deleted
        for (ScratchContainer& scratch : get_memory_manager().scratch_)
        {
            // If there are no neighbors and memory is owned by the memory manager, deallocate memory
            if ((scratch.data() != nullptr) && !scratch.external_mem()) delete[] scratch.data();

            // Nullify scratch array
            scratch.nullify();
        }

        // Restore default number of scratch arrays
        get_memory_manager().scratch_.resize(defaults::DEFAULT_NUM_SCRATCH_CONTAINERS);
    }  

   /***********************************************************************************************
    * void assign_borrower(ScratchContainer&, ScratchMemoryBorrower*)
    * 
    * Assigns a ScratchContainer-ScratchMemoryBorrower pairing. 
    * 
    * Parameters:
    *     container        ScratchContainer object to assign borrower to
    *     borrower         ScratchMemoryBorrower object to assign to ScratchContainer
    ***********************************************************************************************/
    inline void MemoryManager::assign_borrower(ScratchContainer*      container, 
                                        ScratchMemoryBorrower* borrower)
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::assign_borrower..." << std::endl;
        #endif
        // If the container is null, return (throw error?)
        if (container == nullptr) 
        {
            //#ifdef MATL_DEBUG
            throw std::runtime_error("MemoryManager::assign_borrower(ScratchContainer&, ScratchMemoryBorrower*):  Container is null.");
            //#endif
            return;
        }

        // If the borrower is null, return (throw error?)
        if (borrower == nullptr) 
        {
            //#ifdef MATL_DEBUG
            throw std::runtime_error("MemoryManager::assign_borrower(ScratchContainer&, ScratchMemoryBorrower*):  Borrower is null.");
            //#endif
            return;
        }

        // If the pairing already exists, return
        if ((container->borrower() == borrower) && (borrower->scratch_ == container)) return;

        // If the scratch array is already being used, unpair it from its borrower
        if (container->borrower() != nullptr)
            throw std::runtime_error("MemoryManager::assign_borrower(ScratchContainer&, ScratchMemoryBorrower*):  Container is already is use.");

        // If the borrower is already using scratch memory, just exchange data between the 
        // containers and update the pointers.
        if (borrower->get_scratch_container() != nullptr) // (Maybe this should be disallowed)
        {    
            #ifdef MATL_DEBUG
            std::cout << "borrower->get_scratch_container() != nullptr" << std::endl;
            #endif
            ScratchContainer& old_container = *borrower->get_scratch_container();

            index_t    old_bytes_used = old_container.bytes_used();
            index_t    old_size       = old_container.size();
            std::byte* old_data       = old_container.data();
            ScratchContainer* old_prev_neighbor = old_container.prev_neighbor();
            ScratchContainer* old_next_neighbor = old_container.next_neighbor();
            bool old_external_mem = old_container.external_mem();

            // Transfer data to old container
            old_container.bytes_used_    = container->bytes_used();
            old_container.size_          = container->size();
            old_container.data_          = container->data();
            old_container.prev_neighbor_ = container->prev_neighbor();
            old_container.next_neighbor_ = container->next_neighbor();
            old_container.external_mem_  = container->external_mem();

            // Transfer data to new container
            container->bytes_used_       = old_bytes_used;
            container->size_             = old_size;
            container->data_             = old_data;
            container->prev_neighbor_    = old_prev_neighbor;
            container->next_neighbor_    = old_next_neighbor;
            container->external_mem_     = old_external_mem;

            // Set container-borrower relationships 
            // Modify data directly to avoid recursive calls to assign_borrower
            old_container.borrower_ = nullptr;
            container->borrower_    = borrower;
            borrower->scratch_      = container;

            return;
        }
        // If the prospective borrower is using external or private heap memory, throw error
        // (Possibly consider copying data to scratch memory in the future?)
        if (borrower->external_mem() || (borrower->bytes_used() > 0 && borrower->scratch_ == nullptr && borrower->data() != nullptr))
        {
            throw std::runtime_error("MemoryManager::assign_borrower(ScratchContainer&, ScratchMemoryBorrower*):  Borrower is using external or private memory.");
            return;
        }

        // Set pointers to borrower and container
        container->borrower_ = borrower;
        borrower->scratch_   = container;
    } 




   /***********************************************************************************************
    * bool enabled() const
    * 
    * Check if memory manager is enabled.
    * 
    * Returns:
    *     True if memory manager is enabled, false otherwise
    ***********************************************************************************************/
    inline bool MemoryManager::enabled() const
    {
        return get_memory_manager().ENABLED;
    } 


   /***********************************************************************************************
    * void enable()
    *
    * Enable memory manager.
    ***********************************************************************************************/
    inline void MemoryManager::enable()
    {
        get_memory_manager().ENABLED = true;
    }
    
   /***********************************************************************************************
    * void disable()
    * 
    * Disable memory manager.
    ***********************************************************************************************/
    inline void MemoryManager::disable()
    {
        get_memory_manager().ENABLED = false;
    }  


   /***********************************************************************************************
    * void restore_defaults()
    * 
    * Restore default settings for memory manager.
    ***********************************************************************************************/
    inline void MemoryManager::restore_defaults()
    {
        MemoryManager::ENABLED = defaults::ENABLE_MEMORY_MANAGER;

        // Only grow the container list. Shrinking would destroy containers that may still be
        // lending memory or referenced as neighbors.
        if (this->scratch_.size() < defaults::DEFAULT_NUM_SCRATCH_CONTAINERS)
            this->scratch_.resize(defaults::DEFAULT_NUM_SCRATCH_CONTAINERS);
        MemoryManager::SCRATCH_BLOCK_SIZE = defaults::SCRATCH_BLOCK_SIZE;
    }  




   /***********************************************************************************************
    *                                 Matrix Allocations
    ***********************************************************************************************/  



   /***********************************************************************************************
    * void allocate(RowMajorMatrix<T>&, index_t, index_t)
    * 
    * Allocate memory for a RowMajorMatrix object. If the MemoryManager is enabled, memory will be
    * allocated from scratch memory. If the MemoryManager is disabled, memory will be allocated
    * from the heap with a pointer owned by the matrix.
    * 
    * Parameters:
    *     matrix        RowMajorMatrix object to allocate memory for
    *     rows          Number of rows in matrix
    *     cols          Number of columns in matrix
    ***********************************************************************************************/  
    template<class T>
    void MemoryManager::allocate(RowMajorMatrix<T>& matrix, index_t rows, index_t cols)
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::allocate..." << std::endl;
        #endif
        // If memory manager is enabled and ALLOCATE_MATRICES_IN_SCRATCH is set to 'true', 
        // assign scratch memory to the matrix
        if (get_memory_manager().enabled() && MemoryManager::ALLOCATE_MATRICES_IN_SCRATCH)
            get_memory_manager().assign_memory(matrix, rows, cols);
        // Otherwise, allocate privately-owned memory from the heap
        else
        {
            matrix.data_ = new T[rows*cols];
            matrix.rows_ = rows;
            matrix.cols_ = cols;
            static_cast<ScratchMemoryBorrower&>(matrix).external_mem_ = false;
            static_cast<ScratchMemoryBorrower&>(matrix).scratch_      = nullptr;
        }
    }


   /***********************************************************************************************
    * void allocate(FixedSizeRowMajorMatrix<T>&, index_t, index_t) 
    * 
    * The 'allocate' routine for 'FixedSizeRowMajorMatrix' objects is present to support more
    * seemless polymorphism when using mixed matrix types. Because the size of the matrix is 
    * constant and determined at compile time, it cannot be altered. This function checks that
    * the requested size is the same as the size of the matrix and returns if it is. If the
    * requested size is different, an error is thrown.
    *
    * Parameters:
    *     matrix        FixedSizeRowMajorMatrix object
    *     rows          Number of rows in matrix
    *     cols          Number of columns in matrix
    ***********************************************************************************************/
    template<class T, index_t M, index_t N>
    void MemoryManager::allocate(FixedSizeRowMajorMatrix<T, M, N>& /*matrix*/, index_t rows, index_t cols)
    {
        // Check if requested size is the same as the size of the matrix
        if (rows != M || cols != N)
            throw std::runtime_error("MemoryManager::allocate(FixedSizeRowMajorMatrix<T, M, N>&, index_t, index_t):  Requested size does not match fixed size of matrix.");
        
        // If requested size is the same as the size of the matrix, return
        return;
    }



   /***********************************************************************************************
    *                                 Matrix Deallocations
    ***********************************************************************************************/

   /***********************************************************************************************
    * void deallocate(RowMajorMatrix<T>&)
    * 
    * Deallocate memory for a RowMajorMatrix object.
    * 
    * Parameters:
    *     matrix        RowMajorMatrix object to deallocate memory for
    ***********************************************************************************************/
    template<class T>
    void MemoryManager::deallocate(RowMajorMatrix<T>& matrix)
    {
        #ifdef MATL_DEBUG
        std::cout << "In MemoryManager::deallocate..." << std::endl;
        #endif
        // If memory is owned by the memory manager, free memory
        if (matrix.get_scratch_container() != nullptr)
        {
            ScratchContainer& scratch = *matrix.get_scratch_container();

            // Erases borrower and container pointers, resets bytes used, and triggers
            // the matrix to nullify its data pointer and set its size to 0 
            get_memory_manager().free_memory(scratch);
        }
        // If memory is externally owned reset matrix, do not deallocate memory
        else if (matrix.external_mem_){} // Fall through
        // If memory is privately owned, deallocate memory
        else
            delete[] matrix.data_;

        // Reset matrix data
        matrix.data_ = nullptr;
        matrix.rows_ = 0;
        matrix.cols_ = 0;

        // Reset borrower properties
        static_cast<ScratchMemoryBorrower&>(matrix).external_mem_ = false;
        static_cast<ScratchMemoryBorrower&>(matrix).scratch_      = nullptr;
    } 




   /***********************************************************************************************
    * void allocate(_MATRIX_&, index_t, index_t) 
    * 
    * General purpose call for allocation of a matrix object.
    * 
    * Parameters:
    *     matrix        Matrix object to allocate memory for
    *     rows          Number of rows in matrix
    *     cols          Number of columns in matrix
    ***********************************************************************************************/
   template<class T>
   void allocate(RowMajorMatrix<T>& matrix, index_t rows, index_t cols)
   {MemoryManager::allocate(matrix, rows, cols);}

   /***********************************************************************************************
    * void deallocate(_MATRIX_&) 
    * 
    * General purpose call for deallocation of a matrix object.
    * 
    * Parameters:
    *     matrix        Matrix object to deallocate memory for
    ***********************************************************************************************/
    template<class T>
    void deallocate(RowMajorMatrix<T>& matrix)
    {MemoryManager::deallocate(matrix);}
}