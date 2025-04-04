#ifndef _MATL_SCRATCH_MEMORY_BORROWER_IMPLEMENTATION__HPP
#define _MATL_SCRATCH_MEMORY_BORROWER_IMPLEMENTATION__HPP

namespace matl
{
   /***********************************************************************************************
    * Default constructor 
    ***********************************************************************************************/ 
    ScratchMemoryBorrower::ScratchMemoryBorrower() : scratch_(nullptr), external_mem_(false) 
    {};

   /***********************************************************************************************
    * External memory flag constructor 
    ***********************************************************************************************/
    ScratchMemoryBorrower::ScratchMemoryBorrower(bool external_mem) : scratch_(nullptr), 
                                                                      external_mem_(external_mem) 
    {};

   /***********************************************************************************************
    * Move constructor
    * 
    * This method should be called by the Derived class during their move constructor
    * to ensure that pointers to the scratch memory are properly assigned.
    *   
    * Parameters:
    *   other.......other scratch memory borrower
    ***********************************************************************************************/
    ScratchMemoryBorrower::ScratchMemoryBorrower(ScratchMemoryBorrower&& other) : 
        scratch_(other.scratch_),
        external_mem_(other.external_mem_)
    {
        if (other.external_mem())
            throw std::runtime_error("Cannot move a ScratchMemoryBorrower that is externally owned");
        if (other.scratch_ == nullptr)
            throw std::runtime_error("Cannot move a ScratchMemoryBorrower that has no scratch memory");
       
        // Update the container's pointer
        ScratchContainer* scratch = other.scratch_;
        scratch->borrower_ = this;

        // Nullify the other borrower's pointer
        other.scratch_      = nullptr;
        other.external_mem_ = false;
    } 

   /***********************************************************************************************
    * Move assignment
    * 
    * This method should be called by the Derived class during their move assignment operator
    * to ensure that pointers to the scratch memory are properly assigned.
    * 
    * Parameters: 
    *   other.......other scratch memory borrower
    ***********************************************************************************************/ 
    ScratchMemoryBorrower& ScratchMemoryBorrower::operator=(ScratchMemoryBorrower&& other)
    {
        #ifdef MATL_DEBUG
        std::cout << "ScratchMemoryBorrower::operator=(ScratchMemoryBorrower&&) called" << std::endl;
        #endif
        // Check for self assignment
        if (this == &other) return *this;


        if (other.external_mem())
            throw std::runtime_error("Cannot move a ScratchMemoryBorrower that is externally owned");
        if (other.scratch_ == nullptr)
            throw std::runtime_error("Cannot move a ScratchMemoryBorrower that has no scratch memory");
        
        // Set pointer to container
        ScratchContainer* scratch = other.scratch_;
        this->scratch_            = scratch;
        this->external_mem_       = other.external_mem_; // should be 'false'

        // Update the container's pointer
        scratch->borrower_ = this;

        // Nullify the other borrower's pointer
        other.scratch_      = nullptr;
        other.external_mem_ = false; 

        return *this;
    }    

   /***********************************************************************************************
    * ScratchContainer* get_scratch_container()
    * 
    * Returns:
    *     Pointer to the scratch memory container this object is borrowing from
    ***********************************************************************************************/
   ScratchContainer* ScratchMemoryBorrower::get_scratch_container()
   {return this->scratch_;};
   
 
   /***********************************************************************************************
    * void set_scratch_container(ScratchContainer*)
    * 
    * Invokes the MemoryManager to assigns a pairing with the ScratchContainer
    ***********************************************************************************************/
    void ScratchMemoryBorrower::set_scratch_container(ScratchContainer* scratch)
    {MemoryManager::get_memory_manager().assign_borrower(scratch, this);};

   /***********************************************************************************************
    * void nullify_borrowed_memory
    * 
    * Nullify borrowed memory (This should nullify the borrower's pointer to raw array and reset 
    * any size data to 0 to prevent undefined behavior). This method is called whenever the
    * lending ScratchContainer is deallocated, destroyed or the memory is freed.
    ***********************************************************************************************/
    void ScratchMemoryBorrower::nullify_borrowed_memory()
    {this->scratch_ = nullptr;}

   /***********************************************************************************************
    * void set_external_mem(bool)
    * 
    * Sets the external memory flag.
    * 
    * Parameters:
    *    external_mem........external memory flag value
    ***********************************************************************************************/
    void ScratchMemoryBorrower::set_external_mem(bool external_mem)
    {this->external_mem_ = external_mem;};

   /***********************************************************************************************
    * bool external_mem() const
    * 
    * Returns:
    *     True if memory is externally owned, false otherwise
    ***********************************************************************************************/
    bool ScratchMemoryBorrower::external_mem() const
    {return this->external_mem_;};

   /***********************************************************************************************
    * void copy_memory_to_private()
    * 
    * Copies the data from shared scratch memory to newly allocated array owned by the borrower
    ***********************************************************************************************/ 
    void ScratchMemoryBorrower::copy_memory_to_private()
    {this->scratch_ = nullptr;};
}
#endif