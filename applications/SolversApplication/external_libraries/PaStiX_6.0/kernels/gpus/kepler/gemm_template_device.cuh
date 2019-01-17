/*
    -- MAGMA (version 2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Jakub Kurzak
       @author Stan Tomov
       @author Mark Gates
       @author Azzam Haidar
       @author Ahmad Abdelfattah

       See [zcds]gemm_fermi.cu for description of related files.
*/

#ifndef GEMM_TEMPLATE_DEVICE_CUH
#define GEMM_TEMPLATE_DEVICE_CUH

///////////////////////////////////////////////////////////////////////////////////////////////////
// op<trans>( x ) returns x or conj(x).
template< const int conjugate, typename T >
__host__ __device__ static inline
T op( T& x )
{
    if (conjugate == 1) {
        return conj(x);
    } else {
        return x;
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int THR_M, const int THR_N, const int CONJA, const int CONJB>
static __device__
void gemm_template_device_nn(
    int M, int N, int K,
    const T* __restrict__ A, int LDA,
    const T* __restrict__ B, int LDB,
    T*       __restrict__ C, int LDC,
    T alpha, T beta )
{
#if (__CUDA_ARCH__ >= 200)
    int idx = threadIdx.x;  // thread's m dimension
    int idy = threadIdx.y;  // thread's n dimension

    int idt = DIM_X * idy + idx;    // thread's global number

    int idxA = idt % DIM_XA;    // idx within A
    int idyA = idt / DIM_XA;    // idy within A

    int idxB = idt % DIM_XB;    // idx within B
    int idyB = idt / DIM_XB;    // idy within B

    int blx = blockIdx.x;   // block's m dimension
    int bly = blockIdx.y;   // block's n dimension

    __shared__ T sA[BLK_K][BLK_M+1];      // +1 only required if A is transposed
    __shared__ T sB[BLK_N][BLK_K+1];      // +1 always required

    // Registers for the innermost loop
    T rC[THR_N][THR_M];
    T rA[THR_M];
    T rB[THR_N];

    T ra[BLK_K/DIM_YA][BLK_M/DIM_XA];
    T rb[BLK_N/DIM_YB][BLK_K/DIM_XB];

    #ifdef TEXTURE_1D
        int coord_A = offsetA + blx*BLK_M     + idyA*LDA + idxA;
        int coord_B = offsetB + bly*BLK_N*LDB + idyB*LDB + idxB;
    #else
        const T *offs_dA = A + blx*BLK_M     + idyA*LDA + idxA;
        ptrdiff_t boundA = (LDA*(K-1) + M) - ( blx*BLK_M  + idyA*LDA + idxA ) -1;

        const T *offs_dB = B + bly*BLK_N*LDB + idyB*LDB + idxB;
        ptrdiff_t boundB = (LDB*(N-1) + K) - ( bly*BLK_N*LDB + idyB*LDB + idxB ) -1;
    #endif

    int m, n, k, kk;


    // Zero C
    #pragma unroll
    for (n = 0; n < THR_N; n++)
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rC[n][m] = make_FloatingPoint(0.0, 0.0);

    #pragma unroll
    for (n = 0; n < BLK_K; n += DIM_YA)
        #pragma unroll
        for (m = 0; m < BLK_M; m += DIM_XA)
            sA[n+idyA][m+idxA] = fetch(A, m, n, boundA);

    #pragma unroll
    for (n = 0; n < BLK_N; n += DIM_YB)
        #pragma unroll
        for (m = 0; m < BLK_K; m += DIM_XB)
            sB[n+idyB][m+idxB] = fetch(B, m, n, boundB);

    __syncthreads();

    for (kk = 0; kk < K-BLK_K; kk += BLK_K)
    {
        #ifdef TEXTURE_1D
            coord_A += BLK_K*LDA;
            coord_B += BLK_K;
        #else
            offs_dA += BLK_K*LDA;
            boundA  -= BLK_K*LDA;

            offs_dB += BLK_K;
            boundB  -= BLK_K;
        #endif

        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_M/DIM_XA; m++)
                ra[n][m] = fetch(A, m*DIM_XA, n*DIM_YA, boundA);

        #pragma unroll
        for (n = 0; n < BLK_N/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XB; m++)
                rb[n][m] = fetch(B, m*DIM_XB, n*DIM_YB, boundB);

        // Multiply
        #pragma unroll
        for (k = 0; k < BLK_K; k++)
        {
            // Load A shmem->regs
            #pragma unroll
            for (m = 0; m < THR_M; m++)
                rA[m] = sA[k][m*DIM_X+idx];

            // Load B shmem->regs
            #pragma unroll
            for (n = 0; n < THR_N; n++)
                rB[n] = sB[n*DIM_Y+idy][k];

            // Compute
            #pragma unroll
            for (n = 0; n < THR_N; n++) {
                #pragma unroll
                for (m = 0; m < THR_M; m++) {
                    fma(op<CONJA>(rA[m]), op<CONJB>(rB[n]), rC[n][m]);
                }
            }
        }

        __syncthreads();

        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_M/DIM_XA; m++)
                sA[n*DIM_YA+idyA][m*DIM_XA+idxA] = ra[n][m];

        #pragma unroll
        for (n = 0; n < BLK_N/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XB; m++)
                sB[n*DIM_YB+idyB][m*DIM_XB+idxB] = rb[n][m];

        __syncthreads();
    }

    // Multiply last full (BLK_K) or partial block of
    // columns of op(A) and rows of op(B).
    // It's okay that m,n exceed matrix bounds as all work is in registers
    // or shared memory, and out-of-bounds rC[n][m] will not be saved later.
    kk = K - kk;
    #pragma unroll
    for (k = 0; k < kk; k++)
    {
        // Load A shmem->regs
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rA[m] = sA[k][m*DIM_X+idx];

        // Load B shmem->regs
        #pragma unroll
        for (n = 0; n < THR_N; n++)
            rB[n] = sB[n*DIM_Y+idy][k];

        // Compute
        #pragma unroll
        for (n = 0; n < THR_N; n++) {
            #pragma unroll
            for (m = 0; m < THR_M; m++) {
                fma(op<CONJA>(rA[m]), op<CONJB>(rB[n]), rC[n][m]);
            }
        }
    }

    // Store C regs->dev
    #pragma unroll
    for (n = 0; n < THR_N; n++) {
        int coord_dCn = bly*BLK_N + n*DIM_Y + idy;
        #pragma unroll
        for (m = 0; m < THR_M; m++) {
            int coord_dCm = blx*BLK_M + m*DIM_X + idx;
            if (coord_dCm < M && coord_dCn < N) {
                int offsC = coord_dCn*LDC + coord_dCm;

                T &regC = rC[n][m];
                T &memC = C[offsC];

                memC = add(mul(alpha, regC), mul(beta, memC));
            }
        }
    }
#endif /* (__CUDA_ARCH__ >= 200) */
}


///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int THR_M, const int THR_N, const int CONJA, const int CONJB>
static __device__
void gemm_template_device_nt(
    int M, int N, int K,
    const T* __restrict__ A, int LDA,
    const T* __restrict__ B, int LDB,
    T*       __restrict__ C, int LDC,
    T alpha, T beta )
{
#if (__CUDA_ARCH__ >= 200)
    int idx = threadIdx.x;  // thread's m dimension
    int idy = threadIdx.y;  // thread's n dimension

    int idt = DIM_X * idy + idx;    // thread's global number

    int idxA = idt % DIM_XA;    // idx within A
    int idyA = idt / DIM_XA;    // idy within A

    int idxB = idt % DIM_XB;    // idx within B
    int idyB = idt / DIM_XB;    // idy within B

    int blx = blockIdx.x;   // block's m dimension
    int bly = blockIdx.y;   // block's n dimension

    __shared__ T sA[BLK_K][BLK_M+1];      // +1 only required if A is transposed
    __shared__ T sB[BLK_N][BLK_K+1];      // +1 always required

    // Registers for the innermost loop
    T rC[THR_N][THR_M];
    T rA[THR_M];
    T rB[THR_N];

    T ra[BLK_K/DIM_YA][BLK_M/DIM_XA];
    T rb[BLK_K/DIM_YB][BLK_N/DIM_XB];

    #ifdef TEXTURE_1D
        int coord_A = offsetA + blx*BLK_M     + idyA*LDA + idxA;
        int coord_B = offsetB + bly*BLK_N     + idyB*LDB + idxB;
    #else
        const T *offs_dA = A + blx*BLK_M     + idyA*LDA + idxA;
        ptrdiff_t boundA = (LDA*(K-1) + M) - ( blx*BLK_M  + idyA*LDA + idxA ) -1;

        const T *offs_dB = B + bly*BLK_N     + idyB*LDB + idxB;
        ptrdiff_t boundB = (LDB*(K-1) + N) - ( bly*BLK_N     + idyB*LDB + idxB ) -1;
    #endif

    int m, n, k, kk;


    // Zero C
    #pragma unroll
    for (n = 0; n < THR_N; n++)
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rC[n][m] = make_FloatingPoint(0.0, 0.0);

     #pragma unroll
     for (n = 0; n < BLK_K; n += DIM_YA)
         #pragma unroll
         for (m = 0; m < BLK_M; m += DIM_XA)
             sA[n+idyA][m+idxA] = fetch(A, m, n, boundA);

    // Load B dev->shmem
    #pragma unroll
    for (n = 0; n < BLK_K; n += DIM_YB)
        #pragma unroll
        for (m = 0; m < BLK_N; m += DIM_XB)
            sB[m+idxB][n+idyB] = fetch(B, m, n, boundB);

    __syncthreads();

    for (kk = 0; kk < K-BLK_K; kk += BLK_K)
    {
        #ifdef TEXTURE_1D
            coord_A += BLK_K*LDA;
            coord_B += BLK_K*LDB;
        #else
            offs_dA += BLK_K*LDA;
            boundA  -= BLK_K*LDA;

            offs_dB += BLK_K*LDB;
            boundB  -= BLK_K*LDB;
        #endif

        // Load A dev->regs
        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_M/DIM_XA; m++)
                ra[n][m] = fetch(A, m*DIM_XA, n*DIM_YA, boundA);

        // Load B dev->regs
        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_N/DIM_XB; m++)
                rb[n][m] = fetch(B, m*DIM_XB, n*DIM_YB, boundB);

        // Multiply
        #pragma unroll
        for (k = 0; k < BLK_K; k++)
        {
            // Load A shmem->regs
            #pragma unroll
            for (m = 0; m < THR_M; m++)
                rA[m] = sA[k][m*DIM_X+idx];

            // Load B shmem->regs
            #pragma unroll
            for (n = 0; n < THR_N; n++)
                rB[n] = sB[n*DIM_Y+idy][k];

            // Compute
            #pragma unroll
            for (n = 0; n < THR_N; n++) {
                #pragma unroll
                for (m = 0; m < THR_M; m++) {
                    fma(op<CONJA>(rA[m]), op<CONJB>(rB[n]), rC[n][m]);
                }
            }
        }

        __syncthreads();

        // Load A regs->shmem
        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_M/DIM_XA; m++)
                sA[n*DIM_YA+idyA][m*DIM_XA+idxA] = ra[n][m];

        // Load B regs->shmem
        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_N/DIM_XB; m++)
                sB[m*DIM_XB+idxB][n*DIM_YB+idyB] = rb[n][m];
        __syncthreads();
    }

    // Multiply last full (BLK_K) or partial block of
    // columns of op(A) and rows of op(B).
    // It's okay that m,n exceed matrix bounds as all work is in registers
    // or shared memory, and out-of-bounds rC[n][m] will not be saved later.
    kk = K - kk;
    #pragma unroll
    for (k = 0; k < kk; k++)
    {
        // Load A shmem->regs
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rA[m] = sA[k][m*DIM_X+idx];

        // Load B shmem->regs
        #pragma unroll
        for (n = 0; n < THR_N; n++)
            rB[n] = sB[n*DIM_Y+idy][k];

        // Compute
        #pragma unroll
        for (n = 0; n < THR_N; n++) {
            #pragma unroll
            for (m = 0; m < THR_M; m++) {
                fma(op<CONJA>(rA[m]), op<CONJB>(rB[n]), rC[n][m]);
            }
        }
    }

    // Store C regs->dev
    #pragma unroll
    for (n = 0; n < THR_N; n++) {
        int coord_dCn = bly*BLK_N + n*DIM_Y + idy;
        #pragma unroll
        for (m = 0; m < THR_M; m++) {
            int coord_dCm = blx*BLK_M + m*DIM_X + idx;
            if (coord_dCm < M && coord_dCn < N) {
                int offsC = coord_dCn*LDC + coord_dCm;

                T &regC = rC[n][m];
                T &memC = C[offsC];

                memC = add(mul(alpha, regC), mul(beta, memC));
            }
        }
    }
#endif /* (__CUDA_ARCH__ >= 200) */
}


///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int THR_M, const int THR_N, const int CONJA, const int CONJB>
static __device__
void gemm_template_device_tn(
    int M, int N, int K,
    const T* __restrict__ A, int LDA,
    const T* __restrict__ B, int LDB,
    T*       __restrict__ C, int LDC,
    T alpha, T beta )
{
#if (__CUDA_ARCH__ >= 200)
    int idx = threadIdx.x;  // thread's m dimension
    int idy = threadIdx.y;  // thread's n dimension

    int idt = DIM_X * idy + idx;    // thread's global number

    int idxA = idt % DIM_XA;    // idx within A
    int idyA = idt / DIM_XA;    // idy within A

    int idxB = idt % DIM_XB;    // idx within B
    int idyB = idt / DIM_XB;    // idy within B

    int blx = blockIdx.x;   // block's m dimension
    int bly = blockIdx.y;   // block's n dimension

    __shared__ T sA[BLK_K][BLK_M+1];      // +1 only required if A is transposed
    __shared__ T sB[BLK_N][BLK_K+1];      // +1 always required

    // Registers for the innermost loop
    T rC[THR_N][THR_M];
    T rA[THR_M];
    T rB[THR_N];

    // Registers for the dev->shmem copy
    T ra[BLK_M/DIM_YA][BLK_K/DIM_XA];
    T rb[BLK_N/DIM_YB][BLK_K/DIM_XB];

    #ifdef TEXTURE_1D
        int coord_A = offsetA + blx*BLK_M*LDA + idyA*LDA + idxA;
        int coord_B = offsetB + bly*BLK_N*LDB + idyB*LDB + idxB;
    #else
        // bound is the correction to offs_d in order to not get out of memory bound
        // so bound could be negative value since offs_d could be out of bound
            const T *offs_dA = A + blx*BLK_M*LDA + idyA*LDA + idxA;
            ptrdiff_t boundA = (LDA*(M-1) + K) - ( blx*BLK_M*LDA + idyA*LDA + idxA ) -1;

            const T *offs_dB = B + bly*BLK_N*LDB + idyB*LDB + idxB;
            ptrdiff_t boundB = (LDB*(N-1) + K) - ( bly*BLK_N*LDB + idyB*LDB + idxB ) -1;
    #endif

    int m, n, k, kk;

    // Zero C
    #pragma unroll
    for (n = 0; n < THR_N; n++)
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rC[n][m] = make_FloatingPoint(0.0, 0.0);

    // Load A dev->shmem
    #pragma unroll
    for (n = 0; n < BLK_M; n += DIM_YA)
        #pragma unroll
        for (m = 0; m < BLK_K; m += DIM_XA)
            sA[m+idxA][n+idyA] = fetch(A, m, n, boundA);

    #pragma unroll
    for (n = 0; n < BLK_N; n += DIM_YB)
        #pragma unroll
        for (m = 0; m < BLK_K; m += DIM_XB)
            sB[n+idyB][m+idxB] = fetch(B, m, n, boundB);

    __syncthreads();

    for (kk = 0; kk < K-BLK_K; kk += BLK_K)
    {
        #ifdef TEXTURE_1D
            coord_A += BLK_K;
            coord_B += BLK_K;
        #else
            offs_dA += BLK_K;
            boundA  -= BLK_K;

            offs_dB += BLK_K;
            boundB  -= BLK_K;
        #endif

        // Load A dev->regs
        #pragma unroll
        for (n = 0; n < BLK_M/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XA; m++)
                ra[n][m] = fetch(A, m*DIM_XA, n*DIM_YA, boundA);

        // Load B dev->regs
        #pragma unroll
        for (n = 0; n < BLK_N/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XB; m++)
                rb[n][m] = fetch(B, m*DIM_XB, n*DIM_YB, boundB);

        // Multiply
        #pragma unroll
        for (k = 0; k < BLK_K; k++)
        {
            // Load A shmem->regs
            #pragma unroll
            for (m = 0; m < THR_M; m++)
                rA[m] = sA[k][m*DIM_X+idx];

            // Load B shmem->regs
            #pragma unroll
            for (n = 0; n < THR_N; n++)
                rB[n] = sB[n*DIM_Y+idy][k];

            // Compute
            #pragma unroll
            for (n = 0; n < THR_N; n++) {
                #pragma unroll
                for (m = 0; m < THR_M; m++) {
                    fma(op<CONJA>(rA[m]), op<CONJB>(rB[n]), rC[n][m]);
                }
            }
        }

        __syncthreads();

        // Load A regs->shmem
        #pragma unroll
        for (n = 0; n < BLK_M/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XA; m++)
                sA[m*DIM_XA+idxA][n*DIM_YA+idyA] = ra[n][m];

        // Load B regs->shmem
        #pragma unroll
        for (n = 0; n < BLK_N/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XB; m++)
                sB[n*DIM_YB+idyB][m*DIM_XB+idxB] = rb[n][m];

        __syncthreads();
    }

    // Multiply last full (BLK_K) or partial block of
    // columns of op(A) and rows of op(B).
    // It's okay that m,n exceed matrix bounds as all work is in registers
    // or shared memory, and out-of-bounds rC[n][m] will not be saved later.
    kk = K - kk;
    #pragma unroll
    for (k = 0; k < kk; k++)
    {
        // Load A shmem->regs
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rA[m] = sA[k][m*DIM_X+idx];

        // Load B shmem->regs
        #pragma unroll
        for (n = 0; n < THR_N; n++)
            rB[n] = sB[n*DIM_Y+idy][k];

        // Compute
        #pragma unroll
        for (n = 0; n < THR_N; n++) {
            #pragma unroll
            for (m = 0; m < THR_M; m++) {
                fma(op<CONJA>(rA[m]), op<CONJB>(rB[n]), rC[n][m]);
            }
        }
    }

    // Store C regs->dev
    #pragma unroll
    for (n = 0; n < THR_N; n++) {
        int coord_dCn = bly*BLK_N + n*DIM_Y + idy;
        #pragma unroll
        for (m = 0; m < THR_M; m++) {
            int coord_dCm = blx*BLK_M + m*DIM_X + idx;
            if (coord_dCm < M && coord_dCn < N) {
                int offsC = coord_dCn*LDC + coord_dCm;

                T &regC = rC[n][m];
                T &memC = C[offsC];

                memC = add(mul(alpha, regC), mul(beta, memC));
            }
        }
    }
#endif /* (__CUDA_ARCH__ >= 200) */
}


///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, const int DIM_X, const int DIM_Y, const int BLK_M, const int BLK_N, const int BLK_K,
         const int DIM_XA, const int DIM_YA, const int DIM_XB, const int DIM_YB,
         const int THR_M, const int THR_N, const int CONJA, const int CONJB>
static __device__
void gemm_template_device_tt(
    int M, int N, int K,
    const T* __restrict__ A, int LDA,
    const T* __restrict__ B, int LDB,
    T*       __restrict__ C, int LDC,
    T alpha, T beta )
{
#if (__CUDA_ARCH__ >= 200)
    int idx = threadIdx.x;  // thread's m dimension
    int idy = threadIdx.y;  // thread's n dimension

    int idt = DIM_X * idy + idx;    // thread's global number

    int idxA = idt % DIM_XA;    // idx within A
    int idyA = idt / DIM_XA;    // idy within A

    int idxB = idt % DIM_XB;    // idx within B
    int idyB = idt / DIM_XB;    // idy within B

    int blx = blockIdx.x;   // block's m dimension
    int bly = blockIdx.y;   // block's n dimension

    __shared__ T sA[BLK_K][BLK_M+1];      // +1 only required if A is transposed
    __shared__ T sB[BLK_N][BLK_K+1];      // +1 always required

    // Registers for the innermost loop
    T rC[THR_N][THR_M];
    T rA[THR_M];
    T rB[THR_N];

    // Registers for the dev->shmem copy
    T ra[BLK_M/DIM_YA][BLK_K/DIM_XA];
    T rb[BLK_K/DIM_YB][BLK_N/DIM_XB];

    #ifdef TEXTURE_1D
        int coord_A = offsetA + blx*BLK_M*LDA + idyA*LDA + idxA;
        int coord_B = offsetB + bly*BLK_N     + idyB*LDB + idxB;
    #else
        // bound is the correction to offs_d in order to not get out of memory bound
        // so bound could be negative value since offs_d could be out of bound
        const T *offs_dA = A + blx*BLK_M*LDA + idyA*LDA + idxA;
        ptrdiff_t boundA = (LDA*(M-1) + K) - ( blx*BLK_M*LDA + idyA*LDA + idxA ) -1;

        const T *offs_dB = B + bly*BLK_N     + idyB*LDB + idxB;
        ptrdiff_t boundB = (LDB*(K-1) + N) - ( bly*BLK_N     + idyB*LDB + idxB ) -1;
    #endif

    int m, n, k, kk;

    // Zero C
    #pragma unroll
    for (n = 0; n < THR_N; n++)
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rC[n][m] = make_FloatingPoint(0.0, 0.0);

    // Load A dev->shmem
    #pragma unroll
    for (n = 0; n < BLK_M; n += DIM_YA)
        #pragma unroll
        for (m = 0; m < BLK_K; m += DIM_XA)
            sA[m+idxA][n+idyA] = fetch(A, m, n, boundA);

    // Load B dev->shmem
    #pragma unroll
    for (n = 0; n < BLK_K; n += DIM_YB)
        #pragma unroll
        for (m = 0; m < BLK_N; m += DIM_XB)
            sB[m+idxB][n+idyB] = fetch(B, m, n, boundB);

    __syncthreads();

    for (kk = 0; kk < K-BLK_K; kk += BLK_K)
    {
        #ifdef TEXTURE_1D
            coord_A += BLK_K;
            coord_B += BLK_K*LDB;
        #else
            offs_dA += BLK_K;
            boundA  -= BLK_K;

            offs_dB += BLK_K*LDB;
            boundB  -= BLK_K*LDB;
        #endif

        // Load A dev->regs
        #pragma unroll
        for (n = 0; n < BLK_M/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XA; m++)
                ra[n][m] = fetch(A, m*DIM_XA, n*DIM_YA, boundA);

        // Load B dev->regs
        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_N/DIM_XB; m++)
                rb[n][m] = fetch(B, m*DIM_XB, n*DIM_YB, boundB);

        // Multiply
        #pragma unroll
        for (k = 0; k < BLK_K; k++)
        {
            // Load A shmem->regs
            #pragma unroll
            for (m = 0; m < THR_M; m++)
                rA[m] = sA[k][m*DIM_X+idx];

            // Load B shmem->regs
            #pragma unroll
            for (n = 0; n < THR_N; n++)
                rB[n] = sB[n*DIM_Y+idy][k];

            // Compute
            #pragma unroll
            for (n = 0; n < THR_N; n++) {
                #pragma unroll
                for (m = 0; m < THR_M; m++) {
                    fma(op<CONJA>(rA[m]), op<CONJB>(rB[n]), rC[n][m]);
                }
            }
        }

        __syncthreads();

        // Load A regs->shmem
        #pragma unroll
        for (n = 0; n < BLK_M/DIM_YA; n++)
            #pragma unroll
            for (m = 0; m < BLK_K/DIM_XA; m++)
                sA[m*DIM_XA+idxA][n*DIM_YA+idyA] = ra[n][m];

        // Load B regs->shmem
        #pragma unroll
        for (n = 0; n < BLK_K/DIM_YB; n++)
            #pragma unroll
            for (m = 0; m < BLK_N/DIM_XB; m++)
                sB[m*DIM_XB+idxB][n*DIM_YB+idyB] = rb[n][m];

        __syncthreads();
    }

    // Multiply last full (BLK_K) or partial block of
    // columns of op(A) and rows of op(B).
    // It's okay that m,n exceed matrix bounds as all work is in registers
    // or shared memory, and out-of-bounds rC[n][m] will not be saved later.
    kk = K - kk;
    #pragma unroll
    for (k = 0; k < kk; k++)
    {
        // Load A shmem->regs
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rA[m] = sA[k][m*DIM_X+idx];

        // Load B shmem->regs
        #pragma unroll
        for (n = 0; n < THR_N; n++)
            rB[n] = sB[n*DIM_Y+idy][k];

        // Compute
        #pragma unroll
        for (n = 0; n < THR_N; n++) {
            #pragma unroll
            for (m = 0; m < THR_M; m++) {
                fma(op<CONJA>(rA[m]), op<CONJB>(rB[n]), rC[n][m]);
            }
        }
    }

    // Store C regs->dev
    #pragma unroll
    for (n = 0; n < THR_N; n++) {
        int coord_dCn = bly*BLK_N + n*DIM_Y + idy;
        #pragma unroll
        for (m = 0; m < THR_M; m++) {
            int coord_dCm = blx*BLK_M + m*DIM_X + idx;
            if (coord_dCm < M && coord_dCn < N) {
                int offsC = coord_dCn*LDC + coord_dCm;

                T &regC = rC[n][m];
                T &memC = C[offsC];

                memC = add(mul(alpha, regC), mul(beta, memC));
            }
        }
    }
#endif /* (__CUDA_ARCH__ >= 200) */
}
///////////////////////////////////////////////////////////////////////////////////////////////////

#endif //GEMM_TEMPLATE_DEVICE_CUH
