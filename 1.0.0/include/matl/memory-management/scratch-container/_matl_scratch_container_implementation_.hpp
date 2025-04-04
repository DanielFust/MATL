/**************************************************************************************************
 * File: _matl_scratch_container_implementation_.hpp
 * 
 * Description: Implementation for the scratch array class.
 **************************************************************************************************/
#ifndef _MATL_SCRATCH_CONTAINER_IMPLEMENTATION__HPP
#define _MATL_SCRATCH_CONTAINER_IMPLEMENTATION__HPP
namespace matl
{
   /***********************************************************************************************
    * Default Constructor
    ***********************************************************************************************/
    ScratchContainer::ScratchContainer() : size_(0), 
                                           bytes_used_(0), 
                                           borrower_(nullptr), 
                                           prev_neighbor_(nullptr), 
                                           next_neighbor_(nullptr),
                                           external_mem_(false)
    {}
    

   /***********************************************************************************************
    * Sized Constructor
    ***********************************************************************************************/ 
    ScratchContainer::ScratchContainer(index_t size) : size_(size), 
                                                       bytes_used_(0), 
                                                       borrower_(nullptr), 
                                                       prev_neighbor_(nullptr), 
                                                       next_neighbor_(nullptr),
                                                       external_mem_(false)
    {
        // Allocate memory
        data_ = new std::byte[size_];
    }

   /***********************************************************************************************
    * Copy Constructor
    * 
    * Creates a shallow copy. Only one ScratchContainer object should own the memory at a time.
    * 
    * Parameters:
    *     other        ScratchContainer object to copy
    ***********************************************************************************************/
    ScratchContainer::ScratchContainer(const ScratchContainer& other) : size_(other.size_), 
                                                                        bytes_used_(0), 
                                                                        borrower_(nullptr), 
                                                                        prev_neighbor_(nullptr), 
                                                                        next_neighbor_(nullptr),
                                                                        external_mem_(false)
    {
        // Allocate memory
        data_ = new std::byte[size_];
    }  


   /***********************************************************************************************
    * Move Constructor
    * 
    * Parameters:
    *     other        ScratchContainer object to move 
    ***********************************************************************************************/
    ScratchContainer::ScratchContainer(ScratchContainer&& other) : size_(other.size_), 
                                                                   bytes_used_(other.bytes_used_), 
                                                                   borrower_(other.borrower_), 
                                                                   prev_neighbor_(other.prev_neighbor_), 
                                                                   next_neighbor_(other.next_neighbor_),
                                                                   external_mem_(other.external_mem_)
    {
        // Move data
        this->data_ = other.data_;

        // Reset source
        other.nullify();
    }

   /***********************************************************************************************
    * Preallocated Constructor
    * 
    * Construct a scratch array container for a preallocated memory chunk. If memory is externally 
    * owned, this must be set outside of the constructor.
    * 
    * Parameters:
    *     data         Pointer to preallocated memory
    *     size         Size of preallocated memory
    ***********************************************************************************************/
    ScratchContainer::ScratchContainer(std::byte* data, index_t size) : data_(data), 
                                                                        size_(size), 
                                                                        bytes_used_(0), 
                                                                        borrower_(nullptr), 
                                                                        prev_neighbor_(nullptr), 
                                                                        next_neighbor_(nullptr),
                                                                        external_mem_(true)
    {} 

   /***********************************************************************************************
    * Destructor
    ***********************************************************************************************/ 
    ScratchContainer::~ScratchContainer()
    {
        // If the borrowing object has not been deleted, notify it so it can copy the data
        if (borrower_ != nullptr) 
        {
            this->borrower_->copy_memory_to_private();
            #ifdef MATL_DEBUG
            std::cout << "ScratchContainer destructor: copying data from scratch array to borrower." << std::endl;
            #endif
        }

        // If there are previous and next neighbors, notify them
        if (this->prev_neighbor_ != nullptr && this->next_neighbor_ != nullptr)
        {
            // Donate memory to previous neighbor
            this->prev_neighbor_->size_ += this->size_;

            // Update relationships
            this->prev_neighbor_->next_neighbor_ = this->next_neighbor_;
            this->next_neighbor_->prev_neighbor_ = this->prev_neighbor_;
            return;
        }
        // Notify previous neighbor
        else if (this->prev_neighbor_ != nullptr) 
        {
            // Donate memory
            this->prev_neighbor_->size_ += this->size_;

            // Update relationship
            this->prev_neighbor_->next_neighbor_ = nullptr;
            return;
        }
        // Notify next neighbor
        else if (this->next_neighbor_ != nullptr)
        {
            // If the next neighbor is lending memory, the data must be copied to the new start
            // location of the neighbor's array
            if (this->next_neighbor_->borrower_ != nullptr)
            {
                // Copy data from next neighbor to this array
                for (index_t i = 0; i < this->next_neighbor_->bytes_used(); i++)
                    this->data_[i] = this->next_neighbor_->data_[i];    
            }
            
            // Donate memory and set data pointer
            this->next_neighbor_->size_ += this->size_;
            this->next_neighbor_->data_  = this->data_;

            // Update relationship
            this->next_neighbor_->prev_neighbor_ = nullptr;
            return;
        }
    
        // The memory can be deallocated only if it is not externally owned and has no neighbors
        // Neighbor condition should already be met if the destructor has reached this point.
        if ( (this->data_ != nullptr) 
             && !this->external_mem()
             && (this->prev_neighbor_ == nullptr)
             && (this->next_neighbor_ == nullptr) ) 
        {delete[] this->data_;}     
    }


   /***********************************************************************************************
    * Assignment Operator
    * 
    * Parameters:
    *     other        ScratchContainer object to copy 
    * Returns:
    *     Reference to this object
    ***********************************************************************************************/
    ScratchContainer& ScratchContainer::operator=(const ScratchContainer& other)
    {
        // Copy info
        this->size_          = other.size_;
        this->bytes_used_    = other.bytes_used_;
        this->borrower_      = other.borrower_;

        // Newly allocated array has no contiguous neighbors
        this->prev_neighbor_ = nullptr;
        this->next_neighbor_ = nullptr;
        this->external_mem_  = false;

        // Allocate and copy data
        if (other.data_ != nullptr)
        {
            this->data_ = new std::byte[size_];
            for (index_t i = 0; i < size_; i++) data_[i] = other.data_[i];
        }
        else
        {
            this->data_ = nullptr;
        }

        return *this;
    }
    
   /***********************************************************************************************
    * Move Assignment Operator
    * 
    * Parameters:
    *     other        ScratchContainer object to move 
    * Returns:
    *     Reference to this object
    ***********************************************************************************************/
    ScratchContainer& ScratchContainer::operator=(ScratchContainer&& other)
    {
        // Move info
        this->size_          = other.size_;
        this->bytes_used_    = other.bytes_used_;
        this->borrower_      = other.borrower_;
        this->prev_neighbor_ = other.prev_neighbor_;
        this->next_neighbor_ = other.next_neighbor_;
        this->external_mem_  = other.external_mem_;

        // Move data
        this->data_ = other.data_;

        // Reset source
        other.nullify();

        return *this;
    } 

   /***********************************************************************************************
    * size()
    * 
    * Returns:
    *     size of scratch array
    ***********************************************************************************************/ 
    index_t ScratchContainer::size() const
    {return this->size_;}

   /***********************************************************************************************
    * bytes_used()
    * 
    * Returns:
    *     The number of bytes in use
    ***********************************************************************************************/ 
    index_t ScratchContainer::bytes_used() const
    {return this->bytes_used_;}

   /***********************************************************************************************
    * bool external_mem()
    * 
    * Returns:
    *     'True' if the memory is externally owned, 'False' otherwise
    ***********************************************************************************************/ 
    bool ScratchContainer::external_mem() const
    {return this->external_mem_;} 

   /***********************************************************************************************
    * data()
    * 
    * Returns:
    *     Pointer to the raw byte array
    ***********************************************************************************************/ 
    std::byte* ScratchContainer::data() const
    {return this->data_;}


   /***********************************************************************************************
    * prev_neighbor()
    * 
    * Returns:
    *     Pointer to container for array that is contiguous in memory preceding data_ pointer
    ***********************************************************************************************/
    ScratchContainer* ScratchContainer::prev_neighbor()
    {return this->prev_neighbor_;}

   /***********************************************************************************************
    * next_neighbor()
    * 
    * Returns:
    *     Pointer to container for array that is contiguous in memory following data_ pointer
    ***********************************************************************************************/
    ScratchContainer* ScratchContainer::next_neighbor()
    {return this->next_neighbor_;}


   /***********************************************************************************************
    * void* borrower()
    * 
    * Returns:
    *     Pointer to the object matrix using the scratch memory
    ***********************************************************************************************/
    ScratchMemoryBorrower* ScratchContainer::borrower() const
    {return this->borrower_;}
    
   


   /***********************************************************************************************
    * void set_prev_neighbor(ScratchContainer*) 
    * 
    * Sets pointer to ScratchContainer with contiguous data preceding this array's data
    ***********************************************************************************************/
    void ScratchContainer::set_prev_neighbor(ScratchContainer* prev)
    {this->prev_neighbor_ = prev;}

   /***********************************************************************************************
    * void set_next_neighbor(ScratchContainer*) 
    * 
    * Sets pointer to ScratchContainer with contiguous data following this array's data
    ***********************************************************************************************/
    void ScratchContainer::set_next_neighbor(ScratchContainer* next)
    {this->prev_neighbor_ = next;} 

   /***********************************************************************************************
    * void set_borrower(ScratchMemoryBorrower* borrower) 
    * 
    * Sets pointer to object using the scratch memory.
    ***********************************************************************************************/
    void ScratchContainer::set_borrower(ScratchMemoryBorrower* borrower) 
    {
        MemoryManager& manager = MemoryManager::get_memory_manager();
        manager.assign_borrower(this, borrower);   
    }


   /***********************************************************************************************
    * void set_borrower(ScratchMemoryBorrower* borrower, index_t size)
    * 
    * Sets pointer to object using the scratch memory and the size of the memory in use and number
    * of bytes used by the borrower.
    * 
    * Parameters:
    *    borrower        Pointer to object using the scratch memory 
    *    size            Number of bytes required by the borrower
    ***********************************************************************************************/
    void ScratchContainer::set_borrower(ScratchMemoryBorrower* borrower, index_t size)
    {
        // Check if ScratchContainer is large enough
        if (size > this->size_)
        {
            // Throw error if ScratchContainer is too small
            throw std::runtime_error("ScratchContainer::set_borrower(ScratchMemoryBorrower*, index_t):  ScratchContainer is too small for borrower.");
            return;
        }

        // First assume no bytes used
        this->set_borrower(borrower);

        // Set bytes used
        this->bytes_used_ = size;       
    } 


   /***********************************************************************************************
    * void nullify()
    * 
    * Resets the scratch array to an empty state. This does not deallocate or free memory.
    * 
    * Note: This method should be called judiciously as it invalidates pointers without dealing 
    *       with the memory they point to, which can lead to memory leaks if not handled properly.
    *       ONLY CALL THIS METHOD IF THE MEMORY IS ALREADY ACCOUNTED FOR!
    ***********************************************************************************************/
    void ScratchContainer::nullify()
    {   
        // Reset data
        this->data_ = nullptr;

        // Reset size
        this->size_ = 0;

        // Reset bytes used
        this->bytes_used_ = 0;

        // Reset borrower
        this->borrower_ = nullptr;

        // Reset neighbors
        this->prev_neighbor_ = nullptr;
        this->next_neighbor_ = nullptr;

        // Reset external memory flag
        this->external_mem_ = false;
    }  


   /***********************************************************************************************
    * void merge_contiguous_unused()
    * 
    * Combines all scratch arrays which are contiguous in memory with this array and unused.
    ***********************************************************************************************/
    void ScratchContainer::merge_contiguous_unused()
    {
        #ifdef MATL_DEBUG
        std::cout << "In ScratchContainer::merge_contiguous_unused..." << std::endl;
        #endif
        if (this->bytes_used()==0) 
        {
            this->merge_prev_unused();
            this->merge_next_unused();
        }
    }


   /***********************************************************************************************
    * void merge_prev_unused()
    * 
    * Combines all preceding scratch arrays which are contiguous in memory with this array and 
    * unused.
    ***********************************************************************************************/ 
    void ScratchContainer::merge_prev_unused()
    {
        #ifdef MATL_DEBUG
        std::cout << "In ScratchContainer::merge_prev_unused..." << std::endl;
        #endif
        // Previous neighbor
        ScratchContainer* prev_ptr = this->prev_neighbor();
        if (prev_ptr != nullptr)
        {
            ScratchContainer& prev = *prev_ptr;
            // Trigger recursive merge if previous neighbor is unused
            if (prev.bytes_used() == 0)
            {
                // Recursive merge of previous neighbors
                prev.merge_prev_unused();

                // Join contiguous unused scratch arrays
                this->size_ += prev.size();
                this->data_ =  prev.data_;

                // Reset previous neighbor
                prev.nullify();
            }//if (prev.bytes_used() == 0)
        }//if (prev_ptr != nullptr)
    }

   /***********************************************************************************************
    * void merge_next_unused()
    * 
    * Combines all following scratch arrays which are contiguous in memory with this array and 
    * unused.
    ***********************************************************************************************/ 
    void ScratchContainer::merge_next_unused()
    {
        #ifdef MATL_DEBUG
        std::cout << "In ScratchContainer::merge_next_unused..." << std::endl;
        #endif
        // Next neighbor
        ScratchContainer* next_ptr = this->next_neighbor();
        if (next_ptr != nullptr)
        {
            ScratchContainer& next = *next_ptr;
            // Trigger recursive merge if next neighbor is unused
            if (next.bytes_used() == 0)
            {
                // Recursive merge of next neighbors
                next.merge_next_unused();

                // Join contiguous unused scratch arrays (Only size needs to be updated here)
                this->size_ += next.size();

                // Reset next neighbor
                next.nullify();
            }//if (next.bytes_used() == 0)
        }//if (next_ptr != nullptr)
    }

   /***********************************************************************************************
    * void nullify_neighbors()
    * 
    * Nullifies all contiguous neighbors
    ***********************************************************************************************/
    void ScratchContainer::nullify_neighbors()
    {
        this->nullify_prev_neighbors();
        this->nullify_next_neighbors();
    }
    
   /***********************************************************************************************
    * void nullify_next_neighbors()
    *  
    * Nullifies all contiguous neighbors following this scratch array
    * **********************************************************************************************/
    void ScratchContainer::nullify_next_neighbors()
    {
        // Next neighbor
        ScratchContainer* next_ptr = this->next_neighbor();
        if (next_ptr != nullptr)
        {
            ScratchContainer& next = *next_ptr;
            // Recursive nullify of next neighbors
            next.nullify_next_neighbors();

            // Reset next neighbor
            next.nullify();
        }//if (next_ptr != nullptr)
    }

   /***********************************************************************************************
    * void nullify_prev_neighbors()
    * 
    * Nullifies all contiguous neighbors preceding this scratch array
    ***********************************************************************************************/
    void ScratchContainer::nullify_prev_neighbors()
    {
        // Previous neighbor
        ScratchContainer* prev_ptr = this->prev_neighbor();
        if (prev_ptr != nullptr)
        {
            ScratchContainer& prev = *prev_ptr;
            // Recursive nullify of previous neighbors
            prev.nullify_prev_neighbors();

            // Reset previous neighbor
            prev.nullify();
        }//if (prev_ptr != nullptr)
    } 

}
#endif