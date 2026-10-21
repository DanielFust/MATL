
// ScratchContainer must be complete before the MemoryManager prototype, which stores
// ScratchContainers in a std::deque (std::deque does not support incomplete element types)
#include "./scratch-container/_matl_scratch_container_prototype_.hpp"
#include "./temporary-arena/_matl_temporary_arena_prototype_.hpp"
#include "./memory-manager/_matl_memory_manager_prototype_.hpp"
