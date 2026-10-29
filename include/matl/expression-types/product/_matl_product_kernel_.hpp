/**************************************************************************************************
 * File: _matl_product_kernel_.hpp
 *
 * Description: Cache-blocked matrix product kernel.
 *
 * The product C = A*B (A is MxK, B is KxN) is computed in the structure of GotoBLAS/BLIS:
 *
 *   for each block of NC columns of B and C
 *     for each block of KC rows of B (columns of A)
 *       pack the KC x NC block of B into micro-panels of NR columns
 *       for each block of MC rows of A and C
 *         pack the MC x KC block of A into micro-panels of MR rows
 *         for each MR x NR tile of C in the block
 *           micro-kernel: accumulate the tile in registers over the KC products
 *
 * Packing copies the blocks into contiguous buffers, in the order the micro-kernel reads them, so
 * the innermost loops read memory sequentially from cache. With OpenMP, large products run on
 * several threads: the row blocks of C are divided among them (each packs its own blocks of A,
 * and all share the packed block of B), or, when there are fewer row blocks than threads, the
 * column panels of each block. Each element of C is computed by one thread in the same order, so
 * the result does not depend on the number of threads. The packed block of A stays in the L2
 * cache while it is multiplied with every micro-panel of B, and each micro-panel of B stays in the
 * L1 cache while it is multiplied with every micro-panel of A. The packing buffers are evaluation
 * temporaries (see TemporaryBlock), so they are recycled and make no heap allocations.
 *
 * Each element of C is the sum of its K products added in increasing k, as in the simple
 * i-j-k loop: a tile's partial sums are loaded back from C when the next block of KC products is
 * added to them. Without fused multiply-add (FMA) instructions, the kernel therefore gives exactly
 * the result of the simple loop. When the target has FMA (e.g. -march=native on recent x86, or
 * ARM), GCC and Clang fuse each multiply and add by default, which rounds once instead of twice;
 * the kernel then does one FMA per product in increasing k, but the simple loop may be fused
 * differently, so the two can differ in the last bits. -ffp-contract=off disables fusion.
 **************************************************************************************************/
#ifndef MATL_PRODUCT_KERNEL_HPP
#define MATL_PRODUCT_KERNEL_HPP

// Promise to the compiler that pointers do not alias (lets it keep values in registers)
#if defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER)
#define MATL_RESTRICT __restrict
#else
#define MATL_RESTRICT
#endif

namespace matl::detail
{
   /***********************************************************************************************
    * SIMD vectors
    *
    * With GCC and Clang, the micro-kernel computes with SIMD vectors (the 'vector_size'
    * extension) of the target's register width: 32 bytes with AVX, otherwise 16 bytes (SSE2 on
    * x86-64, NEON on ARM). float and double use SIMD vectors; other element types, and other
    * compilers, use scalar code.
    ***********************************************************************************************/
    #if defined(__AVX__)
    inline constexpr std::size_t simd_width = 32;
    #else
    inline constexpr std::size_t simd_width = 16;
    #endif

    template<typename T>
    inline constexpr bool use_simd_v =
    #if defined(__GNUC__) || defined(__clang__)
        std::is_same_v<T, float> || std::is_same_v<T, double>;
    #else
        false;
    #endif

   /***********************************************************************************************
    * struct product_blocking<T>
    *
    * Block sizes of the product kernel for element type T, in elements. MR x NR is the tile of C
    * held in registers (each row of the tile is NR / LANES SIMD vectors); an MR x KC micro-panel
    * of A and a KC x NR micro-panel of B fit in the L1 cache, an MC x KC block of A in the L2
    * cache, and a KC x NC block of B in the L3 cache.
    ***********************************************************************************************/
    template<typename T>
    struct product_blocking
    {
        static constexpr index_t LANES = use_simd_v<T> ? simd_width / sizeof(T) : 1;
        static constexpr index_t MR = 4;
        static constexpr index_t NR = std::max<index_t>(4, LANES);
        static constexpr index_t MC = 64;
        static constexpr index_t KC = 256;
        static constexpr index_t NC = 2048;

        // Products with fewer multiply-adds (M*N*K) than this use the simple loop: packing costs
        // more than it saves for small matrices (measured crossover: about 10 x 10 x 10)
        static constexpr index_t MIN_BLOCKED_WORK = 10 * 10 * 10;
    };

    // Smallest multiple of 'multiple' that is at least n
    constexpr index_t round_up(index_t n, index_t multiple)
    {return ((n + multiple - 1) / multiple) * multiple;}


   /***********************************************************************************************
    * void pack_lhs(packed, lhs, i0, mc, k0, kc)
    *
    * Packs rows [i0, i0+mc) and columns [k0, k0+kc) of the left operand (0-based) into micro-panels
    * of MR rows. Within a micro-panel, the MR elements of each column are consecutive. Rows past
    * the end of the block are zero.
    ***********************************************************************************************/
    template<typename T, typename LHS>
    void pack_lhs(T* MATL_RESTRICT packed, LHS const& lhs, index_t i0, index_t mc, index_t k0, index_t kc)
    {
        constexpr index_t MR = product_blocking<T>::MR;
        for (index_t ir = 0; ir < mc; ir += MR)
        {
            const index_t mr = std::min(MR, mc - ir);
            for (index_t k = 0; k < kc; ++k)
            {
                for (index_t r = 0; r < mr; ++r) packed[r] = static_cast<T>(lhs(i0 + ir + r + 1, k0 + k + 1));
                for (index_t r = mr; r < MR; ++r) packed[r] = T(0);
                packed += MR;
            }
        }
    }

   /***********************************************************************************************
    * void pack_rhs(packed, rhs, k0, kc, j0, nc)
    *
    * Packs rows [k0, k0+kc) and columns [j0, j0+nc) of the right operand (0-based) into
    * micro-panels of NR columns. Within a micro-panel, the NR elements of each row are
    * consecutive. Columns past the end of the block are zero.
    ***********************************************************************************************/
    template<typename T, typename RHS>
    void pack_rhs(T* MATL_RESTRICT packed, RHS const& rhs, index_t k0, index_t kc, index_t j0, index_t nc)
    {
        constexpr index_t NR = product_blocking<T>::NR;
        for (index_t jr = 0; jr < nc; jr += NR)
        {
            const index_t nr = std::min(NR, nc - jr);
            for (index_t k = 0; k < kc; ++k)
            {
                for (index_t c = 0; c < nr; ++c) packed[c] = static_cast<T>(rhs(k0 + k + 1, j0 + jr + c + 1));
                for (index_t c = nr; c < NR; ++c) packed[c] = T(0);
                packed += NR;
            }
        }
    }

   /***********************************************************************************************
    * void micro_kernel(kc, a, b, acc)
    *
    * acc += (MR x kc micro-panel a) * (kc x NR micro-panel b). The MR x NR accumulators are kept in
    * registers; the fixed loop bounds let the compiler unroll and vectorize the update.
    ***********************************************************************************************/
    template<typename T>
    inline void micro_kernel(index_t kc, T const* MATL_RESTRICT a, T const* MATL_RESTRICT b,
                             T (&acc)[product_blocking<T>::MR][product_blocking<T>::NR])
    {
        constexpr index_t MR = product_blocking<T>::MR;
        constexpr index_t NR = product_blocking<T>::NR;

        #if defined(__GNUC__) || defined(__clang__)
        if constexpr (use_simd_v<T>)
        {
            // Each row of the tile is V SIMD vectors (in registers). Memory is read and written
            // through an unaligned vector type, which may alias T.
            constexpr index_t LANES = product_blocking<T>::LANES;
            constexpr index_t V     = NR / LANES;
            typedef T simd           __attribute__((vector_size(simd_width)));
            typedef T simd_unaligned __attribute__((vector_size(simd_width), aligned(alignof(T))));

            simd sum[MR][V];
            for (index_t r = 0; r < MR; ++r)
                for (index_t v = 0; v < V; ++v)
                    sum[r][v] = *reinterpret_cast<simd_unaligned const*>(&acc[r][v * LANES]);

            // (The loops over the tile are unrolled completely, so the tile stays in registers)
            for (index_t k = 0; k < kc; ++k)
            {
                simd bk[V];
                #pragma GCC unroll 16
                for (index_t v = 0; v < V; ++v)
                    bk[v] = *reinterpret_cast<simd_unaligned const*>(b + v * LANES);
                #pragma GCC unroll 16
                for (index_t r = 0; r < MR; ++r)
                    #pragma GCC unroll 16
                    for (index_t v = 0; v < V; ++v)
                        sum[r][v] += a[r] * bk[v];
                a += MR;
                b += NR;
            }

            for (index_t r = 0; r < MR; ++r)
                for (index_t v = 0; v < V; ++v)
                    *reinterpret_cast<simd_unaligned*>(&acc[r][v * LANES]) = sum[r][v];
        }
        else
        #endif
        {
            // Local accumulators (not aliased by the panels, so they stay in registers)
            T sum[MR][NR];
            for (index_t r = 0; r < MR; ++r)
                for (index_t c = 0; c < NR; ++c)
                    sum[r][c] = acc[r][c];

            for (index_t k = 0; k < kc; ++k)
            {
                for (index_t r = 0; r < MR; ++r)
                    for (index_t c = 0; c < NR; ++c)
                        sum[r][c] += a[r] * b[c];
                a += MR;
                b += NR;
            }

            for (index_t r = 0; r < MR; ++r)
                for (index_t c = 0; c < NR; ++c)
                    acc[r][c] = sum[r][c];
        }
    }

   /***********************************************************************************************
    * void blocked_product<T, SUBTRACT>(target, lhs, rhs, M, N, K)
    *
    * target = lhs * rhs, or with SUBTRACT, target -= lhs * rhs (without a temporary: the packed
    * blocks of lhs are negated, and the sums start from the target's elements). lhs is M x K and
    * rhs is K x N. The operands are read through operator()(i,j) (1-based) only while packing. T is
    * the element type of the product and of the target, and must be an arithmetic type (the
    * packing buffers are raw memory).
    *
    * Parameters:
    *     target       Target matrix (M x N; must not share elements with the operands)
    *     lhs, rhs     Operands (matrices or evaluated operands, see product_operand_t)
    *     M, N, K      Dimensions
    ***********************************************************************************************/
    template<typename T, bool SUBTRACT = false, WritableMatrix MAT, typename LHS, typename RHS>
    void blocked_product(MAT& target, LHS const& lhs, RHS const& rhs, index_t M, index_t N, index_t K)
    {
        constexpr index_t MR = product_blocking<T>::MR;
        constexpr index_t NR = product_blocking<T>::NR;
        constexpr index_t MC = product_blocking<T>::MC;
        constexpr index_t KC = product_blocking<T>::KC;
        constexpr index_t NC = product_blocking<T>::NC;

        // Several threads for large products (see _matl_parallel_.hpp); every element of the
        // target is computed by one thread, in the same order as serially
        const bool parallel = use_threads(M * N * K, PRODUCT_PARALLEL_THRESHOLD);
        const index_t row_blocks = (M + MC - 1) / MC;

        // Packing buffers, sized for the largest blocks of this product (each thread that packs
        // blocks of lhs has its own)
        const index_t kc_max = std::min(K, KC);
        const index_t lhs_block_bytes = kc_max * round_up(std::min(M, MC), MR) * sizeof(T);
        TemporaryBlock packed_rhs_block(kc_max * round_up(std::min(N, NC), NR) * sizeof(T));
        TemporaryBlock packed_lhs_block(lhs_block_bytes);
        T* const packed_rhs = reinterpret_cast<T*>(packed_rhs_block.data());
        T* const packed_lhs = reinterpret_cast<T*>(packed_lhs_block.data());

        // Packs the block of lhs at rows i0.., columns k0.. (negated with SUBTRACT)
        const auto pack_block = [&lhs](T* packed, index_t i0, index_t mc, index_t k0, index_t kc)
        {
            pack_lhs(packed, lhs, i0, mc, k0, kc);
            if constexpr (SUBTRACT)
                for (index_t p = 0; p < kc * round_up(mc, MR); ++p) packed[p] = -packed[p];
        };

        // The micro-tiles of target rows i0..i0+mc-1, in the column panels [jr_first, jr_last)
        // (in units of NR) of the column block starting at j0
        const auto compute = [&](T const* packed, index_t i0, index_t mc, index_t j0, index_t nc, index_t k0,
                                 index_t kc, index_t jr_first, index_t jr_last)
        {
            for (index_t jr = jr_first * NR; jr < std::min(nc, jr_last * NR); jr += NR)
            {
                const index_t nr = std::min(NR, nc - jr);
                for (index_t ir = 0; ir < mc; ir += MR)
                {
                    const index_t mr = std::min(MR, mc - ir);

                    // Start from 0 for the first KC products, otherwise from the partial sums
                    // stored in C (so the products are added in increasing k); with SUBTRACT,
                    // always from C
                    T acc[MR][NR] = {};
                    if (SUBTRACT || k0 > 0)
                        for (index_t r = 0; r < mr; ++r)
                            for (index_t c = 0; c < nr; ++c)
                                acc[r][c] = target(i0 + ir + r + 1, j0 + jr + c + 1);

                    micro_kernel(kc, packed + ir * kc, packed_rhs + jr * kc, acc);

                    for (index_t r = 0; r < mr; ++r)
                        for (index_t c = 0; c < nr; ++c)
                            target(i0 + ir + r + 1, j0 + jr + c + 1) = acc[r][c];
                }
            }
        };

        for (index_t j0 = 0; j0 < N; j0 += NC)
        {
            const index_t nc = std::min(NC, N - j0);
            const index_t panels = (nc + NR - 1) / NR;
            for (index_t k0 = 0; k0 < K; k0 += KC)
            {
                const index_t kc = std::min(KC, K - k0);
                pack_rhs(packed_rhs, rhs, k0, kc, j0, nc);

                if (parallel && row_blocks >= thread_count())
                {
                    // Enough row blocks: each thread packs and computes its own row blocks
                    parallel_ranges_if(true, row_blocks, [&](index_t first, index_t last)
                    {
                        TemporaryBlock own_block(lhs_block_bytes);
                        T* const own = reinterpret_cast<T*>(own_block.data());
                        for (index_t b = first; b < last; ++b)
                        {
                            const index_t i0 = b * MC, mc = std::min(MC, M - i0);
                            pack_block(own, i0, mc, k0, kc);
                            compute(own, i0, mc, j0, nc, k0, kc, 0, panels);
                        }
                    }, equal_split{row_blocks});
                }
                else
                {
                    // Few row blocks (e.g. a 64-row product): the threads share each packed row
                    // block, and divide its column panels
                    for (index_t i0 = 0; i0 < M; i0 += MC)
                    {
                        const index_t mc = std::min(MC, M - i0);
                        pack_block(packed_lhs, i0, mc, k0, kc);
                        parallel_ranges_if(parallel, panels, [&](index_t first, index_t last)
                        {
                            compute(packed_lhs, i0, mc, j0, nc, k0, kc, first, last);
                        }, equal_split{panels});
                    }
                }
            }
        }
    }
}
#endif
