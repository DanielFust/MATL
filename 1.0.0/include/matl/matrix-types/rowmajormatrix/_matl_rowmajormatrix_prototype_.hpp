namespace matl
{
   /***********************************************************************************************
    * class RowMajorMatrix
    * 
    * Row major matrix class  
    * 
    * Template Arguments:
    *     T     Matrix data type     
    ***********************************************************************************************/ 
    template<class T>
    class RowMajorMatrix : public _MATRIX_< RowMajorMatrix<T> >, 
                           public ScratchMemoryBorrower
    {
      private:
        index_t rows_;
        index_t cols_;
        T*      data_;

      public:
        
        typedef T dtype;

        /********************* Constructors *********************/
        
        // Default constructor
        RowMajorMatrix();

        // Sized constructor
        RowMajorMatrix(index_t, index_t);

        // Preallocated constructor
        RowMajorMatrix(index_t, index_t, T*);

        // Copy constructor (General matrix)
        template<typename MAT>
        RowMajorMatrix(_MATRIX_<MAT> const&);

        // Copy constructor (Same type)
        RowMajorMatrix(RowMajorMatrix<T> const&);

        // Move constructor
        RowMajorMatrix(RowMajorMatrix<T>&&);

        // Destructor
        ~RowMajorMatrix();

        /************* Assignment and Copying ******************/

        // Expression execution trigger
        template<typename EXPR>
        RowMajorMatrix<T>& operator= (_EXPRESSION_<EXPR> const&);

        #if 0
        // Expression execution trigger (matrix multiplication specialization)
        template<typename MAO1, typename MAO2>
        RowMajorMatrix<T>& operator= (_MATRIX_TIMES_MATRIX_<MAO1,MAO2> const&);
        #endif

        // RowMajorMatrix move assignment
        RowMajorMatrix<T>& operator= (RowMajorMatrix<T>&&);

        // Copy assignment/constructor for RowMajorMatrix
        template<class T2>
        RowMajorMatrix<T>& operator= (RowMajorMatrix<T2>&);

        // Copy assignment/constructor for General Matrix
        template<typename MAT>
        RowMajorMatrix<T>& operator= (_MATRIX_<MAT>&);

        // Scalar assignement shortcut
        template<class T2, enable_if_scalar<T2> >
        RowMajorMatrix<T>& operator= (const T2);

        // Addition assignment (by scalar)
        template<class T2, enable_if_scalar<T2> = true >
        RowMajorMatrix<T>& operator+= (const T2);

        // Subtraction assignment (by scalar)
        template<class T2, enable_if_scalar<T2> = true >
        RowMajorMatrix<T>& operator-= (const T2);

        // Multiplication assignment (by scalar)
        template<class T2, enable_if_scalar<T2> = true >
        RowMajorMatrix<T>& operator*= (const T2);

        // Division assignment (by scalar)
        template<class T2, enable_if_scalar<T2> = true >
        RowMajorMatrix<T>& operator/= (const T2);

        // Addition assignment (by matrix)
        template<typename MAT>
        RowMajorMatrix<T>& operator+= (_MATRIX_<MAT>&);

        // Subtraction assignment (by matrix)
        template<typename MAT>
        RowMajorMatrix<T>& operator-= (_MATRIX_<MAT>&);

        // Multiplication assignment (by matrix)
        template<typename MAT>
        RowMajorMatrix<T>& operator*= (_MATRIX_<MAT>&);

        // Addition assignment (by matrix expression)
        template<typename EXPR>
        RowMajorMatrix<T>& operator+= (_EXPRESSION_<EXPR> const&);

        // Subtraction assignment (by matrix expression)
        template<typename EXPR>
        RowMajorMatrix<T>& operator-= (_EXPRESSION_<EXPR> const&);

        // Multiplication assignment (by matrix expression)
        template<typename EXPR>
        RowMajorMatrix<T>& operator*= (_EXPRESSION_<EXPR> const&);

        // General matrix copy
        template<typename MAT>
        void copy(const _MATRIX_<MAT>&);

        // RowMajorMatrix copy
        template<class T2>
        void copy(const RowMajorMatrix<T2>&);

        /********************* Data Access *********************/

        // Single index getter
        dtype operator() (index_t) const;

        // Double index getter
        dtype operator() (index_t, index_t) const;

        // Single index accessor
        dtype& operator() (index_t);

        // Double index accessor
        dtype& operator() (index_t, index_t);

        /********************** Utilities **********************/

        // Data writing shortcut
        template<class T2, enable_if_scalar<T2> = true>
        void operator= (const T2);

        // Row Getter
        index_t rows() const;

        // Coloumn Getter
        index_t cols() const;

        // Size Getter
        index_t size() const;

        // Resize (This might be antithecal to the design)
        //void resize(index_t, index_t);

        // Print
        void print() const;

        // Print with newline
        void println() const;

        /*********** ScratchMemoryBorrower interface ***********/

      protected: // Only for use by MemoryManager

        // Get Scratch Memory Container (probably doesn't need an override)
        ScratchContainer* get_scratch_container();

        // Set Scratch Memory Container (probably doesn't need an override)
        void set_scratch_container(ScratchContainer*) override;

        // Copy memory to privately owned memory 
        void copy_memory_to_private() override;

        // Nullify borrowed memory and reset matrix
        void nullify_borrowed_memory() override;

        // Get raw data pointer
        std::byte* data() const override;

        // Set external memory flag
        void set_external_mem(bool) override;

      public: // const ScratchMemoryBorrower methods  

        // Get the number of bytes in use
        index_t bytes_used() const override;

        // Is the memory externally owned (Not managed by MemoryManager)
        bool external_mem() const;

        // Give privileged access to memory manager
        friend class MemoryManager;
    };
}
