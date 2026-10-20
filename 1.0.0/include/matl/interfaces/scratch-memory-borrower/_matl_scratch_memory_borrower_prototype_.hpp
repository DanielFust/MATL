/**************************************************************************************************
 * File: _matl_scratch_memory_borrower_.hpp
 * 
 * Description: virtual base class for MATL objects that utilize scratch memory borrower class.
 **************************************************************************************************/
#ifndef _MATL_SCRATCH_MEMORY_BORROWER__HPP
#define _MATL_SCRATCH_MEMORY_BORROWER__HPP
namespace matl
{
   /***********************************************************************************************
    * class ScratchMemoryBorrower
    * 
    * Interface base class for MATL objects that utilize scratch memory borrower class.
    ***********************************************************************************************/  
    class ScratchMemoryBorrower
    {
      protected:
        /****************** Data *****************/
        
        // Scratch Memory Container
        ScratchContainer* scratch_;

        // Is the memory externally owned (not by the MemoryManager)
        bool external_mem_;

        /*************** Constructors ************/
        
        // Default Constructor
        ScratchMemoryBorrower();

        // External Memory Flag Constructor
        ScratchMemoryBorrower(bool external_mem);

        // Virtual destructor (derived matrix types are polymorphic through this interface)
        virtual ~ScratchMemoryBorrower() = default;

        /************ Move Semantics *************/
        // The move contructor and move assignment operator should be called by the derived class
        // during their move constructor and move assignment operator to ensure that pointers to the
        // scratch memory are properly assigned.

        // Move Constructor
        ScratchMemoryBorrower(ScratchMemoryBorrower&&);

        // Move Assignment Operator
        ScratchMemoryBorrower& operator=(ScratchMemoryBorrower&&);

        /***************** Methods ***************/

        // Get Scratch Memory Container
        ScratchContainer* get_scratch_container();
        
        // Set Scratch Memory Container
        virtual void set_scratch_container(ScratchContainer*);

        // Nullify borrowed memory (This should nullify the borrower's pointer to raw array and reset 
        // any size data to 0 to prevent undefined behavior). This method is called whenever the
        // lending ScratchContainer is deallocated, destroyed or the memory is freed. 
        virtual void nullify_borrowed_memory();

        // Get raw data pointer 
        virtual std::byte* data() const = 0;

        // Set external memory flag
        virtual void set_external_mem(bool);

      public:

        // Get the number of bytes in use
        virtual index_t bytes_used() const = 0;

        // Is the memory externally owned
        bool external_mem() const;

        // Copy memory to privately owned memory
       /* 
          This method should be overridden by derived classes and should
          terminate by nullifying the pointer to the scratch memory container. 
        */  
        virtual void copy_memory_to_private() = 0;
    
        // Only the MemoryManager, ScratchContainer and derived classes should call  non-const
        // ScratchMemoryBorrower methods 
        friend class MemoryManager;
        friend class ScratchContainer;
    };
}
#endif