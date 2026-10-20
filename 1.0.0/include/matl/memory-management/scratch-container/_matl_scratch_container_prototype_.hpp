/**************************************************************************************************
 * File: _matl_scratch_container_prototype_.hpp
 * 
 * Description: Prototype for the scratch array class.
 **************************************************************************************************/
#ifndef _MATL_SCRATCH_CONTAINER_PROTOTYPE__HPP
#define _MATL_SCRATCH_CONTAINER_PROTOTYPE__HPP
 namespace matl
{
   /***********************************************************************************************
    * class ScratchContainer
    * 
    * ScratchContainer objects are used to store temporary data for matrix expressions owned by the 
    * MATL MemoryManager. Scratch containers can be used to recycle memory used in a compound
    * matrix expression or to preallocate memory for matrices and matrix expressions.
    ***********************************************************************************************/
    class ScratchContainer
    {
      private:  
        // Scratch memory array
        std::byte* data_;

        // Scratch memory size
        index_t size_;

        // Scratch memory in use
        index_t bytes_used_;

        // Pointer to previous contiguous memory scratch array
        ScratchContainer* prev_neighbor_;

        // Pointer to next contiguous memory scratch array
        ScratchContainer* next_neighbor_;

        // Is the memory externally owned
        bool external_mem_;

        // Pointer to matrix or object using array
        ScratchMemoryBorrower* borrower_;

       /***********************************************************
        * The nullification functions do not deallocate or free
        * memory. The memory MUST already be accounted for by the
        * MemoryManager or memory leaks will occur.
        ***********************************************************/
       
        // Nullifies the scratch array pointers and sets data to unused state
        void nullify();

        // Nullify all neighbors using contiguous memory
        void nullify_neighbors();

        // Nullify following neighbors using contiguous memory
        void nullify_next_neighbors();

        // Nullify previous neighbors using contiguous memory
        void nullify_prev_neighbors();

        // Point neighbors and borrower at this object's address (used after a move)
        void repoint_relationships();

      public:

        // Default Constructor
        ScratchContainer();

        // Sized-Constructor
        ScratchContainer(index_t size);

        // Preallocated Constructor
        ScratchContainer(std::byte* data, index_t size);

        // Copy Constructor
        ScratchContainer(const ScratchContainer& scratch_array);

        // Move Constructor
        ScratchContainer(ScratchContainer&& scratch_array) noexcept;

        // Destructor
        ~ScratchContainer();

        // Assignment Operator
        ScratchContainer& operator=(const ScratchContainer& scratch_array);

        // Move Assignment Operator
        ScratchContainer& operator=(ScratchContainer&& scratch_array) noexcept;

        // Get the size of the scratch array
        index_t size() const;

        // Get the number of bytes in use
        index_t bytes_used() const;

        // Is the memory externally owned
        bool external_mem() const;

        // Get the data pointer
        std::byte* data() const;

        // Get the previous contiguous memory scratch array
        ScratchContainer* prev_neighbor();

        // Get the next contiguous memory scratch array
        ScratchContainer* next_neighbor();

        // Get addess of the object using the scratch array
        ScratchMemoryBorrower* borrower() const;

        // Set the previous contiguous memory scratch array
        void set_prev_neighbor(ScratchContainer* prev_neighbor);

        // Set the next contiguous memory scratch array
        void set_next_neighbor(ScratchContainer* next_neighbor);

        // Set the object using the scratch array
        void set_borrower(ScratchMemoryBorrower* borrower);

        // Set the object using the scratch array, and number of required elements
        void set_borrower(ScratchMemoryBorrower* borrower, index_t size);

        // Merge contiguous unused scratch arrays
        void merge_contiguous_unused();

        // Merge previous unused scratch arrays
        void merge_prev_unused();

        // Merge next unused scratch arrays
        void merge_next_unused();

        // Allow MemoryManager to access private members
        friend class MemoryManager;

        // Allow ScratchMemoryBorrower to access private members
        friend class ScratchMemoryBorrower;
    };
}
#endif
 