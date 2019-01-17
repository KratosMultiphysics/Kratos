//on garde tous les tests sur M et N

///////////////////////////////////////////////////////////////////////////////////////////////////

// size of work for a thread
#define THR_M ( BLK_M / DIM_X )
#define THR_N ( BLK_N / DIM_Y )

///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef KERNEL_BOTTOM
#  define KERNEL_BOTTOM_NAME _bottom
#else
#  define KERNEL_BOTTOM_NAME
#endif

#ifdef KERNEL_RIGHT
#  define KERNEL_RIGHT_NAME _right
#else
#  define KERNEL_RIGHT_NAME
#endif

#ifdef KERNEL_LDLT
#  define KERNEL_SHORT_NAME gemdm
#else
#  define KERNEL_SHORT_NAME gemm
#endif

#define CONCAT_NAME2(a,b,c,d) a ## b ## c ## d
#define CONCAT_NAME3(a,b,c,d) CONCAT_NAME2(a,b,c,d)
#define CONCAT_NAME(a) CONCAT_NAME3(KERNEL_SHORT_NAME, KERNEL_BOTTOM_NAME, KERNEL_RIGHT_NAME, a)
#ifndef GENERATE_SM_VERSION_KERNEL_NAME
#  error GENERATE_SM_VERSION_KERNEL_NAME must be defined
#endif
#ifndef version
#  error version must be defined
#endif
#if   (version == trans_nn)
#  define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(CONCAT_NAME(_nn))
#elif (version == trans_nt)
#  define TRANS_B
#  define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(CONCAT_NAME(_nt))
#elif (version == trans_nc)
#  define TRANS_B
#  define CONJ_B
#  define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(CONCAT_NAME(_nc))
#elif (version == trans_tn)
#  define TRANS_A
#  define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(CONCAT_NAME(_tn))
#elif (version == trans_tt)
#  define TRANS_A
#  define TRANS_B
#  define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(CONCAT_NAME(_tt))
#elif (version == trans_tc)
#  define TRANS_A
#  define TRANS_B
#  define CONJ_B
#  define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(CONCAT_NAME(_tc))
#elif (version == trans_cn)
#  define TRANS_A
#  define CONJ_A
#  define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(CONCAT_NAME(_cn))
#elif (version == trans_ct)
#  define TRANS_A
#  define CONJ_A
#  define TRANS_B
#  define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(CONCAT_NAME(_ct))
#elif (version == trans_cc)
#  define TRANS_A
#  define CONJ_A
#  define TRANS_B
#  define CONJ_B
#  define kernel_name GENERATE_SM_VERSION_KERNEL_NAME(CONCAT_NAME(_cc))
#endif
#ifndef kernel_name
#  error "kernel_name must be defined"
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__
void     kernel_name (int M, int N, int K,
                      FloatingPoint_t alpha,
                      const FloatingPoint_t *A, int LDA,
#ifdef KERNEL_LDLT
                      const FloatingPoint_t *D, int LDD,
#endif
                      const FloatingPoint_t *B, int LDB,
                      FloatingPoint_t beta,
                      FloatingPoint_t       *C, int LDC,
                      int offsetA,
#ifdef KERNEL_LDLT
                      int offsetD,
#endif
                      int offsetB,
                      int blocknbr, const int *blocktab,
                      int fblocknbr, const int *fblocktab)
{
    int offset[THR_M+1];

    int idx = threadIdx.x;  // thread's m dimension
    int idy = threadIdx.y;  // thread's n dimension

    int idt = DIM_X * idy + idx;    // thread's global number

    int idxA = idt % DIM_XA;    // idx within A
    int idyA = idt / DIM_XA;    // idy within A

    int idxB = idt % DIM_XB;    // idx within B
    int idyB = idt / DIM_XB;    // idy within B
#ifdef KERNEL_BOTTOM
    int blx = M/BLK_M;      // block's m dimension
#else
    int blx = blockIdx.x;   // block's m dimension
#endif
#ifdef KERNEL_RIGHT
    int bly = N/BLK_N;      // block's n dimension
#else
    int bly = blockIdx.y;   // block's n dimension
#endif

    __shared__ FloatingPoint_t sA[BLK_K][BLK_M+1];      // +1 only required if A is transposed
#ifdef KERNEL_LDLT
    __shared__ FloatingPoint_t sD[BLK_K];
#endif
    __shared__ FloatingPoint_t sB[BLK_N][BLK_K+1];      // +1 always required

    // Registers for the innermost loop
    FloatingPoint_t rC[THR_N][THR_M];
    FloatingPoint_t rA[THR_M];
#ifdef KERNEL_LDLT
    FloatingPoint_t rD;
#endif
    FloatingPoint_t rB[THR_N];

#ifdef TRANS_A
    const FloatingPoint_t *offs_dA = A + blx*BLK_M*LDA + idyA*LDA+idxA;
#else
    const FloatingPoint_t *offs_dA = A + blx*BLK_M     + idyA*LDA+idxA;
#endif
#ifdef KERNEL_LDLT
    const FloatingPoint_t *offs_dD = D + idyA*LDD + idyA;
#endif
#ifdef TRANS_B
    const FloatingPoint_t *offs_dB = B + bly*BLK_N     + idyB*LDB+idxB;
#else
    const FloatingPoint_t *offs_dB = B + bly*BLK_N*LDB + idyB*LDB+idxB;
#endif

    int m, n, k, kk;
    int coordm, coordn;
#ifdef KERNEL_LDLT
    int coordd;
#endif
#if (defined KERNEL_BOTTOM || defined KERNEL_RIGHT)
    int kk_aux=-1*BLK_K;
#endif

    // Zero C
#pragma unroll
    for (n = 0; n < THR_N; n++)
#pragma unroll
        for (m = 0; m < THR_M; m++)
            rC[n][m] = make_FloatingPoint(0.0, 0.0);

    for (kk = 0; kk < K-BLK_K; kk += BLK_K)
    {
#if (defined KERNEL_BOTTOM || defined KERNEL_RIGHT)
        kk_aux=kk;
#endif
        // Load A dev->shmem
#ifdef TRANS_A
#pragma unroll
        for (n = 0; n < BLK_M; n += DIM_YA){
#pragma unroll
            for (m = 0; m < BLK_K; m += DIM_XA){
                // TODO : implement this correctly, not require in PaStiX yet
#  if (defined KERNEL_BOTTOM || defined KERNEL_RIGHT)
                coordn = n + blx*BLK_M + idyA;
                if(coordn < M)
#  endif
                    sA[m+idxA][n+idyA] = fetch(A, m, n);
#  if (defined KERNEL_BOTTOM || defined KERNEL_RIGHT)
                else
                    sA[m+idxA][n+idyA] = make_FloatingPoint(0.0,0.0);
#  endif
            }
        }
#  ifdef KERNEL_LDLT
        // Load D dev->shmem
#pragma unroll
        for (n = 0; n < BLK_K; n += DIM_XA)
        {
            sD[n+idyA] = fetch(D, n, n);

        }
#  endif /* KERNEL_LDLT */
#else /* TRANS_A */
        /* Not TRANS_A : It's the normal case in PaStiX */
#pragma unroll
        for (n = 0; n < BLK_K; n += DIM_YA){
#pragma unroll
            for (m = 0; m < BLK_M; m += DIM_XA){
#  ifdef KERNEL_BOTTOM
                coordm = m + blx*BLK_M + idxA;
                //coordn = n + idyA + kk;
                if(coordm < M/* && coordn < K*/)
#  endif
                    sA[n+idyA][m+idxA] = fetch(A, m, n);
#  ifdef KERNEL_BOTTOM
                else
                    sA[n+idyA][m+idxA] = make_FloatingPoint(0.0,0.0);
#  endif
            }
        }
#  ifdef KERNEL_LDLT
        // Load D dev->shmem
#pragma unroll
        for (n = 0; n < BLK_K; n += DIM_YA)
        {
            sD[n+idyA] = fetch(D, n, n);
        }
#  endif
#endif /* TRANS_A */

        // Load B dev->shmem
#ifdef TRANS_B
        /* TRANS_B : It's the normal case in PaStiX */

#pragma unroll
        for (n = 0; n < BLK_K; n += DIM_YB){
#pragma unroll
            for (m = 0; m < BLK_N; m += DIM_XB){
#  ifdef KERNEL_RIGHT
                coordm = m + bly*BLK_N + idxB;
                if(coordm < N)
#  endif
                    sB[m+idxB][n+idyB] = fetch(B, m, n);
#  ifdef KERNEL_RIGHT
                else
                    sB[m+idxB][n+idyB] = make_FloatingPoint(0.0,0.0);
#  endif
            }
        }
#else /* TRANS_B */
        /* TODO : Check, not TRANS_B : It's NOT the case in PaStiX */
#pragma unroll
        for (n = 0; n < BLK_N; n += DIM_YB){
#pragma unroll
            for (m = 0; m < BLK_K; m += DIM_XB){
#  ifdef KERNEL_RIGHT
                coordn = n + bly*BLK_N + idyB;
                if(coordn < N)
#  endif
                    sB[n+idyB][m+idxB] = fetch(B, m, n);
#  ifdef KERNEL_RIGHT
                else
                    sB[n+idyB][m+idxB] = make_FloatingPoint(0.0,0.0);
#  endif
            }
        }
#endif

        __syncthreads();

        // Multiply
#pragma unroll
        for (k = 0; k < BLK_K; k++)
        {
#ifdef KERNEL_LDLT
            rD = sD[k];
#endif
            // Load A shmem->regs
#pragma unroll
            for (m = 0; m < THR_M; m++)
            {
#ifdef KERNEL_LDLT
                rA[m] = mul(sA[k][m*DIM_X+idx],rD);
#else
                rA[m] = sA[k][m*DIM_X+idx];
#endif
            }
            // Load B shmem->regs
#pragma unroll
            for (n = 0; n < THR_N; n++)
                rB[n] = sB[n*DIM_Y+idy][k];

            // Compute
#pragma unroll
            for (n = 0; n < THR_N; n++)
#pragma unroll
                for (m = 0; m < THR_M; m++) {
#ifdef CONJ_A
#ifdef CONJ_B
                    fma(conj(rA[m]), conj(rB[n]), rC[n][m]);
#else
                    fma(conj(rA[m]), rB[n], rC[n][m]);
#endif
#else
#ifdef CONJ_B
                    fma(rA[m], conj(rB[n]), rC[n][m]);
#else
                    fma(rA[m], rB[n], rC[n][m]);
#endif
#endif
                }
        }
        __syncthreads();


        //maj offset
#ifdef TRANS_A
        offs_dA += BLK_K;
#else
        offs_dA += BLK_K*LDA;
#endif
#ifdef KERNEL_LDLT
        offs_dD += BLK_K*LDD + BLK_K;
#endif
#ifdef TRANS_B
        offs_dB += BLK_K*LDB;
#else
        offs_dB += BLK_K;
#endif

        __syncthreads();

    }
    /////////////////////////////////////////////////////////////////////



#  if (defined KERNEL_BOTTOM || defined KERNEL_RIGHT)
    kk_aux+=BLK_K;
#  endif
    // Load A dev->shmem
#ifdef TRANS_A
#pragma unroll
    for (n = 0; n < BLK_M; n += DIM_YA){
#pragma unroll
        for (m = 0; m < BLK_K; m += DIM_XA){
#  ifdef KERNEL_BOTTOM
            coordn = n + blx*BLK_M + idyA;
#  endif
#  if (defined KERNEL_BOTTOM || defined KERNEL_RIGHT)
            coordm = m + idxA + kk_aux;
#  else
            coordm = m + idxA + kk;
#  endif
            if(coordm < K
#  if (defined KERNEL_BOTTOM)
               && coordn < M
#  endif
               )
                sA[m+idxA][n+idyA] = fetch(A, m, n);
            else
                sA[m+idxA][n+idyA] = make_FloatingPoint(0.0,0.0);
        }
    }
#  ifdef KERNEL_LDLT
    // Load D dev->shmem
#pragma unroll
    for (n = 0; n < BLK_K; n += DIM_XA)
    {
        coordd = idyA + kk + n;
        if(coordd<K)
            sD[n+idyA] = fetch(D, n, n);
        else
            sD[n+idyA] = make_FloatingPoint(1.0,0.0);
    }
#endif
#else /* not TRANS_A */
#pragma unroll
    for (n = 0; n < BLK_K; n += DIM_YA){
#pragma unroll
        for (m = 0; m < BLK_M; m += DIM_XA){
#  ifdef KERNEL_BOTTOM
            coordm = m + blx*BLK_M + idxA;
#  endif
#  if (defined KERNEL_BOTTOM || defined KERNEL_RIGHT)
            coordn = n + idyA + kk_aux;
#  else
            coordn = n + idyA + kk;
#  endif
            if(
#  if (defined KERNEL_BOTTOM)
                coordm < M &&
#  endif
                coordn < K)
                sA[n+idyA][m+idxA] = fetch(A, m, n);
            else
                sA[n+idyA][m+idxA] = make_FloatingPoint(0.0,0.0);
        }
    }
#ifdef KERNEL_LDLT
    // Load D dev->shmem
#pragma unroll
    for (n = 0; n < BLK_K; n += DIM_YA)
    {
        coordd = idyA + kk + n;
        if(coordd<K)
            sD[n+idyA] = fetch(D, n, n);
        else
            sD[n+idyA] = make_FloatingPoint(1.0,0.0);
    }
#endif
#endif

    // Load B dev->shmem
#ifdef TRANS_B
#pragma unroll
    for (n = 0; n < BLK_K; n += DIM_YB){
#pragma unroll
        for (m = 0; m < BLK_N; m += DIM_XB){
#  ifdef KERNEL_RIGHT
            coordm = m + bly*BLK_N + idxB;
#  endif
#  if (defined KERNEL_BOTTOM || defined KERNEL_RIGHT)
            coordn = n + idyB + kk_aux;
#  else
            coordn = n + idyB + kk;
#  endif
            if(
#  ifdef KERNEL_RIGHT
                coordm < N &&
#  endif
                coordn < K)
                sB[m+idxB][n+idyB] = fetch(B, m, n);
            else
                sB[m+idxB][n+idyB] = make_FloatingPoint(0.0,0.0);
        }
    }
#else
#pragma unroll
    for (n = 0; n < BLK_N; n += DIM_YB){
#pragma unroll
        for (m = 0; m < BLK_K; m += DIM_XB){
#  if (defined KERNEL_BOTTOM || defined KERNEL_RIGHT)
            coordm = m + idxB + kk_aux;
#  else
            coordm = m + idxB + kk;
#  endif
#  ifdef KERNEL_RIGHT
            coordn = n + bly*BLK_N + idyB;
#  endif
            if(coordm < K
#  ifdef KERNEL_RIGHT
               && coordn < N
#  endif
               )
                sB[n+idyB][m+idxB] = fetch(B, m, n);
            else
                sB[n+idyB][m+idxB] = make_FloatingPoint(0.0,0.0);
        }
    }
#endif

    __syncthreads();

    // Multiply
#pragma unroll
    for (k = 0; k < BLK_K; k++)
    {
#ifdef KERNEL_LDLT
        rD = sD[k];
#endif
        // Load A shmem->regs
#pragma unroll
        for (m = 0; m < THR_M; m++) {
#ifdef KERNEL_LDLT
            rA[m] = mul(sA[k][m*DIM_X+idx],rD);
#else
            rA[m] = sA[k][m*DIM_X+idx];
#endif
        }
        // Load B shmem->regs
#pragma unroll
        for (n = 0; n < THR_N; n++)
            rB[n] = sB[n*DIM_Y+idy][k];

        // Compute
#pragma unroll
        for (n = 0; n < THR_N; n++)
#pragma unroll
            for (m = 0; m < THR_M; m++) {
#ifdef CONJ_A
#ifdef CONJ_B
                fma(conj(rA[m]), conj(rB[n]), rC[n][m]);
#else
                fma(conj(rA[m]), rB[n], rC[n][m]);
#endif
#else
#ifdef CONJ_B
                fma(rA[m], conj(rB[n]), rC[n][m]);
#else
                fma(rA[m], rB[n], rC[n][m]);
#endif
#endif
            }
    }

    __syncthreads();
#if (defined KERNEL_BOTTOM || defined KERNEL_RIGHT)
    //maj offset
#  ifdef TRANS_A
    offs_dA += BLK_K;
#  else
    offs_dA += BLK_K*LDA;
#  endif
    // WARNING: THIS WASN'T HERE BEFORE GENERIC STENCIL
#  ifdef KERNEL_LDLT
    offs_dD += BLK_K*LDD + BLK_K;
#  endif
#  ifdef TRANS_B
    offs_dB += BLK_K*LDB;
#  else
    offs_dB += BLK_K;
#  endif
#endif
    __syncthreads();
    ////////////////////////////////////////////////////////////////////////

    {
#define FROWNUM(tab, b) tab[2*b]
#define LROWNUM(tab, b) tab[2*b+1]
#define BLOCKSIZE(tab, b) LROWNUM(tab, b) - FROWNUM(tab, b) + 1
        int blocknum = 0, fblocknum = 0;
        size_t totalblocksize = 0;
        size_t blocksize = BLOCKSIZE(blocktab, blocknum);
        int    rownum;

        offset[0] = 0;
        for (m = 0; m < THR_M; m++) {
            int coord_dCm = blx*BLK_M + m*DIM_X+idx;

            if (coord_dCm < M) {

                /*
                 * We should keep blocknum < blocknbr
                 */
                while( totalblocksize + blocksize < coord_dCm + 1)
                {
                    totalblocksize += blocksize;
                    blocknum++;
                    blocksize = BLOCKSIZE(blocktab, blocknum);
                }

                /* Global row index */
                rownum = coord_dCm - totalblocksize + FROWNUM(blocktab, blocknum);

                while (LROWNUM(fblocktab, fblocknum) < rownum) {
                    offset[m] += BLOCKSIZE(fblocktab, fblocknum);
                    fblocknum++;
                }
                offset[m+1] = offset[m];
                offset[m] += rownum - FROWNUM(fblocktab, fblocknum);
            }
        }
        __syncthreads();
#undef FROWNUM
#undef LROWNUM
    }



    // Store C regs->dev
#pragma unroll
    for (n = 0; n < THR_N; n++) {
        int coord_dCn = bly*BLK_N + n*DIM_Y+idy;
#pragma unroll
        for (m = 0; m < THR_M; m++) {
            int coord_dCm = blx*BLK_M + m*DIM_X+idx;
            if (coord_dCm < M && coord_dCn < N) {
                int offsC = coord_dCn*LDC + offset[m]; /*coord_dCm;*/

                FloatingPoint_t &regC = rC[n][m];
                FloatingPoint_t &memC = C[offsC];

                memC = add(mul(alpha, regC), mul(beta, memC));
            }
        }
    }
    (void)coordm; (void)coordn;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#undef TRANS_A
#undef TRANS_B
#undef CONJ_A
#undef CONJ_B
#undef KERNEL_BOTTOM_NAME
#undef KERNEL_RIGHT_NAME
#undef KERNEL_SHORT_NAME
#undef CONCAT_NAME
#undef CONCAT_NAME2
#undef CONCAT_NAME3
#undef THR_M
#undef THR_N

#undef kernel_name
