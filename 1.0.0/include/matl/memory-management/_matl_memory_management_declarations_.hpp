namespace matl
{
    // Memory Manager
    class MemoryManager;

    // Stack-ordered memory for evaluation temporaries
    class TemporaryArena;

    // RAII holders for evaluation temporaries
    class TemporaryBlock;
    template<typename MAT> class EvaluationTemporary;
}
