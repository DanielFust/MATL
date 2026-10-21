/**************************************************************************************************
 * File: _matl_evaluation_temporary_.hpp
 *
 * Description: RAII holders for the temporary matrices created while an expression is evaluated.
 **************************************************************************************************/
#ifndef MATL_EVALUATION_TEMPORARY_HPP
#define MATL_EVALUATION_TEMPORARY_HPP
namespace matl::detail
{
   /***********************************************************************************************
    * class TemporaryBlock
    *
    * A block of raw memory for one evaluation temporary. It is taken from this thread's
    * TemporaryArena when MemoryManager::use_temporary_arena() is true, and from the heap
    * otherwise, and returned to the same place when the TemporaryBlock is destroyed.
    ***********************************************************************************************/
    class TemporaryBlock
    {
      private:
        std::byte* data_;
        bool       from_arena_;

      public:
        explicit TemporaryBlock(index_t bytes) :
            data_(nullptr),
            from_arena_(MemoryManager::use_temporary_arena())
        {
            this->data_ = this->from_arena_ ? MemoryManager::temporary_arena().allocate(bytes)
                                            : new std::byte[bytes == 0 ? 1 : bytes];
        }

        ~TemporaryBlock()
        {
            if (this->from_arena_) MemoryManager::temporary_arena().release(this->data_);
            else                   delete[] this->data_;
        }

        TemporaryBlock(TemporaryBlock const&) = delete;
        TemporaryBlock& operator=(TemporaryBlock const&) = delete;

        std::byte* data() const {return this->data_;}
    };


   /***********************************************************************************************
    * class EvaluationTemporary<MAT>
    *
    * A temporary matrix of type MAT, with the given dimensions, that lives for the duration of an
    * evaluation. The general case allocates the matrix with 'allocate' (fixed-size matrices need
    * no allocation).
    ***********************************************************************************************/
    template<typename MAT>
    class EvaluationTemporary
    {
      private:
        MAT mat_;

      public:
        EvaluationTemporary(index_t rows, index_t cols)
        {allocate(this->mat_, rows, cols);}

        EvaluationTemporary(EvaluationTemporary const&) = delete;
        EvaluationTemporary& operator=(EvaluationTemporary const&) = delete;

        MAT&       matrix()       {return this->mat_;}
        MAT const& matrix() const {return this->mat_;}
    };


   /***********************************************************************************************
    * class EvaluationTemporary<RowMajorMatrix<T>>
    *
    * Dynamic temporaries take their memory from a TemporaryBlock (arena or heap) and view it as a
    * RowMajorMatrix with external memory. Used for element types that may live in raw memory
    * (trivially copyable, with alignment no stricter than the arena's); other element types use
    * the general case.
    ***********************************************************************************************/
    template<class T>
    requires (std::is_trivially_copyable_v<T> && alignof(T) <= TemporaryArena::ALIGNMENT)
    class EvaluationTemporary< RowMajorMatrix<T> >
    {
      private:
        TemporaryBlock    block_;
        RowMajorMatrix<T> mat_;

      public:
        EvaluationTemporary(index_t rows, index_t cols) :
            block_(rows * cols * sizeof(T)),
            mat_(rows, cols, reinterpret_cast<T*>(block_.data()))
        {}

        EvaluationTemporary(EvaluationTemporary const&) = delete;
        EvaluationTemporary& operator=(EvaluationTemporary const&) = delete;

        RowMajorMatrix<T>&       matrix()       {return this->mat_;}
        RowMajorMatrix<T> const& matrix() const {return this->mat_;}
    };
}
#endif
