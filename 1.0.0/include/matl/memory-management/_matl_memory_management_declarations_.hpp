namespace matl
{
    // Memory Manager
    class MemoryManager;

    // Scratch Memory Container
    class ScratchContainer;

    // Stack-ordered memory for evaluation temporaries
    class TemporaryArena;

    // RAII holders for evaluation temporaries
    class TemporaryBlock;
    template<typename MAT> class EvaluationTemporary;
}
