!
! @file fmultilap.f90
!
! Fortran 90 multiple rhs example.
!
! This example tries to solve a problem with multiple matrices and
! rhs. It is possible to use multiple pastix instances to solve the
! problems in parallel, and it is possible to force multiple
! sequential solves in parallel with OpenMP instead of mumti-threaded
! pastix solves.
!
! To launch the testing: ./fmultilap < test_config.in
! with test_config.in of the following format:
!
! ---------------------------------------------------------------
! A single comment line to describe the file
! 1                      Enable/disable the verbose mode
! 0                      Enable/disable the check and correct
! 0                      Enable/disable the multi-threaded solves
! 5                      Total number of threads
! 1                      Number of PaStiX instances
! 2                      Number of outermost iterations
! 2                      Number of distinct matrices per iteration
! 10                     Number of righ-hand-side per matrix
! 2                      Number of solve step per problem
! 10                     First dimension of each laplacian matrix
! 10                     Second dimension of each laplacian matrix
! 10                     Third dimension of each laplacian matrix
! ---------------------------------------------------------------
!
! @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.1
! @author Andrea Piacentini
! @author Mathieu Faverge
! @date 2018-07-16
!
program fmultilap
  use iso_c_binding
  use pastix_enums
  use spmf
  use pastixf

  implicit none

  type multilap_param
     logical                    :: checkmat    ! Enable/disable the check and correct on the spm matrix
     logical                    :: multirhs    ! Enable/disable multi-rhs solves
     logical                    :: output      ! Enable output information
     integer                    :: dim1        ! First dimension of the Laplacian matrix
     integer                    :: dim2        ! Second dimension of the Laplacian matrix
     integer                    :: dim3        ! Third dimension of the Laplacian matrix
     integer(kind=pastix_int_t) :: n           ! Size of the matrix
     integer                    :: nnz         ! Number of non zeroes in the matrix
     integer                    :: nb_outit    ! Number of outer iteration
     integer                    :: nb_mat      ! Number of matrices to factorize per iteration
     integer                    :: nb_rhs      ! Number of right hand side per matrix
     integer                    :: nb_solve    ! Number of solve steps per RHS
     integer                    :: nb_sys      ! Number of systems (each couple composed of 1 matrix with 1 subset of rhs)
     integer                    :: nb_thrd     ! Number of thread available
     integer                    :: nb_fact_omp ! Number of PaStiX instances for the factorization
     integer                    :: nb_fact_thd ! Number of threads per PaStiX instances for the factorization
     integer                    :: nb_solv_omp ! Number of PaStiX instances for the solve
     integer                    :: nb_solv_thd ! Number of threads per PaStiX instances for the solve
  end type multilap_param

  type rhs_subset
     integer(kind=pastix_int_t)                            :: idmat ! Index of the associated matrix
     integer(kind=pastix_int_t)                            :: nrhs  ! Number of RHS in the subset
     complex(kind=c_double_complex), dimension(:), pointer :: b0    ! Array of the initial b, to restart at each iteration
     complex(kind=c_double_complex), dimension(:), pointer :: b     ! Array of the b vector
     complex(kind=c_double_complex), dimension(:), pointer :: x     ! Array of the x vector
  end type rhs_subset

  type sys_lin
     integer(kind=pastix_int_t)                            :: pid         ! Index of the pastix instance factorizing the matrix
     integer(kind=pastix_int_t)                            :: nsys        ! Number of subsystems used to solve the rhs
     type(pastix_data_t),                          pointer :: pastix_data ! The pastix_data associated to the matrix
     integer(kind=pastix_int_t),     dimension(iparm_size) :: iparm       ! iparm array associated to the matrix
     real(kind=c_double),            dimension(dparm_size) :: dparm       ! dparm array associated to the matrix
     type(spmatrix_t),                             pointer :: spm         ! Pointer to the spm matrix
     integer(kind=pastix_int_t),     dimension(:), pointer :: rowptr      ! User matrix
     integer(kind=pastix_int_t),     dimension(:), pointer :: colptr      ! User matrix
     complex(kind=c_double_complex), dimension(:), pointer :: values      ! User matrix
     type(rhs_subset),               dimension(:), pointer :: rhs         ! Array of rhs subsets
  end type sys_lin

  type(multilap_param)                              :: params    ! Parameters of the example
  type(sys_lin), dimension(:), allocatable, target  :: sys_array ! Array of all the linear systems to solve
  integer(kind=pastix_int_t), dimension(iparm_size) :: iparm     ! Global iparm array used to initialize the systems
  real(kind=c_double), dimension(dparm_size)        :: dparm     ! Global dparm array used to initialize the systems
  type(sys_lin), pointer                            :: matrix
  type(rhs_subset), pointer                         :: rhs
  type(c_ptr)                                       :: x_ptr, b_ptr
  integer, dimension(:), allocatable                :: ila_thrmn, ila_thrmx, ila_thrsz
  integer(kind=c_int), dimension(:), pointer        :: bindtab
  integer(kind=pastix_int_t)                        :: th, im, ir, i, j, k
  integer(kind=pastix_int_t)                        :: size
  integer(c_int)                                    :: info, ginfo = 0
  !
  integer, parameter :: MULTILAP_ANALYZE_TIME = 1
  integer, parameter :: MULTILAP_FACT_TIME    = 2
  integer, parameter :: MULTILAP_FACT_FLOPS   = 3
  integer, parameter :: MULTILAP_SOLV_TIME    = 4
  integer, parameter :: MULTILAP_MAX_STAT     = MULTILAP_SOLV_TIME
  !
  real(kind=c_double),  allocatable, dimension(:,:)  :: dla_thread_stats
  real(kind=c_double),  dimension(MULTILAP_MAX_STAT) :: dla_final_stats = 0.0

  !
  ! Get the problem confirguration
  !
  call multilap_init( params )

  !
  ! Set some parameters in iparm
  !
  call pastixInitParam( iparm, dparm )

  iparm(IPARM_THREAD_NBR) = params%nb_fact_thd
  if ( iparm(IPARM_THREAD_NBR) == 1 ) then
     iparm(IPARM_SCHEDULER)  = PastixSchedSequential
  else
     iparm(IPARM_SCHEDULER)  = PastixSchedStatic
  end if
  iparm(IPARM_VERBOSE)       = PastixVerboseNot
  iparm(IPARM_FACTORIZATION) = PastixFactLLT

  ! Initialize the matrices structures
  allocate(sys_array(params%nb_mat))
  do im = 1, params%nb_mat
     !
     matrix => sys_array(im)
     !
     matrix%iparm(:) = iparm(:)
     matrix%dparm(:) = dparm(:)
     !

     ! Initialize the solve structures
     if ( params%multirhs ) then
        !
        ! multi-thread RHS is enabled
        ! Distribute the rhs solve in the same manner as the matrices factorizations
        !
        matrix%nsys = 1
        allocate(matrix%rhs( matrix%nsys ))

        rhs => matrix%rhs(1)
        rhs%nrhs = params%nb_rhs

        size = rhs%nrhs * params%n

        if ( params%output ) then
           ! Allocate a backup of b that will be destroyed by the check
           allocate(rhs%b0(size))
        end if
        allocate(rhs%b(size))
        allocate(rhs%x(size))

     else
        !
        ! multi-thread rhs is disabled
        ! Distribute the rhs solve among each sequential thread issued from the factorization
        !
        call split_parall( params%nb_rhs, params%nb_fact_thd )

        matrix%nsys = iparm(IPARM_THREAD_NBR)
        allocate(matrix%rhs( matrix%nsys ))

        do th = 1, iparm(IPARM_THREAD_NBR)

           rhs => matrix%rhs(th)

           rhs%nrhs = ila_thrsz(th)
           size = rhs%nrhs * params%n

           if ( rhs%nrhs .gt. 0 ) then
              if ( params%output ) then
                 ! Allocate a backup of b that will be destroyed by the check
                 allocate(rhs%b0(size))
              end if
              allocate(rhs%b(size))
              allocate(rhs%x(size))
           end if
        end do
     end if
  end do

  if (params%output) then
     write( 6, * )        '!--------------------------------------------------------------------!'
     write( 6, * )        '     Size of x = ', params%n
     write( 6, fmt=8888 ) 'Matrix', 'NRHS', 'Start', 'End'

     do im = 1, params%nb_mat
        matrix => sys_array(im)
        k = 1
        do j = 1, matrix%nsys
           rhs => matrix%rhs(j)
           write( 6, fmt=8889 ) im, rhs%nrhs, k, k + rhs%nrhs - 1
           k = k + rhs%nrhs
        end do
     end do
  end if


  !
  ! Initialization loop to create the PaStiX instances
  !
  call split_parall( params%nb_mat, params%nb_fact_omp )

  !
  !$OMP PARALLEL NUM_THREADS(params%nb_fact_omp) DEFAULT(NONE) &
  !$OMP SHARED(params, sys_array, k)            &
  !$OMP SHARED(ila_thrsz, ila_thrmn, ila_thrmx) &
  !$OMP SHARED(dla_thread_stats, iparm, dparm)  &
  !$OMP PRIVATE(th, im, ir, i, j )              &
  !$OMP PRIVATE(matrix, rhs, x_ptr, b_ptr)      &
  !$OMP PRIVATE(info, bindtab)

  !$OMP DO SCHEDULE(STATIC,1)
  init_loop: do th = 1, params%nb_fact_omp
     !
     ! Give each matrix a subset of cores
     !
     allocate( bindtab( iparm(IPARM_THREAD_NBR) ) )
     do i = 1, iparm(IPARM_THREAD_NBR)
        bindtab( i ) = int((th-1) * iparm(IPARM_THREAD_NBR) + (i - 1), c_int)
     end do

     do im = ila_thrmn(th), ila_thrmx(th)

        matrix => sys_array(im)
        matrix%pid = th

        ! 1- Initialize the parameters and the solver
        call pastixInitWithAffinity( matrix%pastix_data, 0, &
             & matrix%iparm, matrix%dparm, bindtab )

     end do

     deallocate(bindtab)
  end do init_loop
  !$OMP END DO
  !$OMP END PARALLEL

  !
  ! Outmost iteration
  !
  main_loop: do k = 1, params%nb_outit
     !
     write(6,*) '!====================================================================!'
     write(6,*) ' Outer iteration: ', k
     write(6,*) '!--------------------------------------------------------------------!'
     write(6,*) ' Nb of factorization performed in parallel      = ', params%nb_fact_omp
     write(6,*) ' Nb of threads used by PaStiX per factorization = ', params%nb_fact_thd

     !
     ! Distribute factorizations among parallel threads
     !
     call split_parall( params%nb_mat, params%nb_fact_omp )
     !
     allocate(dla_thread_stats(params%nb_mat, MULTILAP_MAX_STAT))
     dla_thread_stats(:,:) = 0.0

     !
     !$OMP PARALLEL NUM_THREADS(params%nb_fact_omp) DEFAULT(NONE) &
     !$OMP SHARED(params, sys_array, k)            &
     !$OMP SHARED(ila_thrsz, ila_thrmn, ila_thrmx) &
     !$OMP SHARED(dla_thread_stats, iparm, dparm)  &
     !$OMP PRIVATE(th, im, ir, i, j )              &
     !$OMP PRIVATE(matrix, rhs, x_ptr, b_ptr)      &
     !$OMP PRIVATE(info, bindtab)

     !$OMP DO SCHEDULE(STATIC,1)
     fact_loop: do th = 1, params%nb_fact_omp
        !
        ! Give each matrix a subset of cores
        !
        do im = ila_thrmn(th), ila_thrmx(th)

           matrix => sys_array(im)

           matrix%pid = th

           ! 0- Generate the SPM matrix
           call multilap_genOneMatrix( params, matrix, k, im )

           ! 2- Analyze the problem
           call pastix_task_analyze( matrix%pastix_data, matrix%spm, info )
           dla_thread_stats(th,MULTILAP_ANALYZE_TIME) = &
                & dla_thread_stats(th,MULTILAP_ANALYZE_TIME) + matrix%dparm(DPARM_ANALYZE_TIME)

           ! 3- Factorize the matrix
           call pastix_task_numfact( matrix%pastix_data, matrix%spm, info )

           dla_thread_stats(th,MULTILAP_FACT_TIME) = &
                & dla_thread_stats(th,MULTILAP_FACT_TIME) + matrix%dparm(DPARM_FACT_TIME)
           dla_thread_stats(th,MULTILAP_FACT_FLOPS) = &
                & dla_thread_stats(th,MULTILAP_FACT_FLOPS) + (matrix%dparm(DPARM_FACT_FLOPS)/(1024.0**3))

        end do

     end do fact_loop
     !$OMP END DO
     !$OMP END PARALLEL

     if ( params%output ) then
        write(6,*) '!--------------------------------------------------------------------!'
        write(6,*) '! Results per matrix'
        write(6,*) '!'
        do im = 1, params%nb_mat
           matrix => sys_array(im)

           write(6,*) ' Matrix ', im, ' done by instance ', matrix%pid
           write(6,*) ' Time for analysis      ', matrix%dparm(DPARM_ANALYZE_TIME)
           write(6,*) ' Pred Time for fact     ', matrix%dparm(DPARM_PRED_FACT_TIME)
           write(6,*) ' Time for factorization ', matrix%dparm(DPARM_FACT_TIME)
           write(6,*) ' GFlops/s for fact      ', matrix%dparm(DPARM_FACT_FLOPS)/(1024.0**3)
        end do

        write(6,*) '!--------------------------------------------------------------------!'
        write(6,*) '! Results per PaStiX instance'
        write(6,*) '!'
        do th = 1, params%nb_fact_omp
           dla_thread_stats(th,MULTILAP_FACT_FLOPS) = &
                & dla_thread_stats(th,MULTILAP_FACT_FLOPS) / ila_thrsz(th)
           write(6,*) ' Thread ',th
           write(6,*) ' Time for analysis      ', dla_thread_stats(th,MULTILAP_ANALYZE_TIME)
           write(6,*) ' Time for factorization ', dla_thread_stats(th,MULTILAP_FACT_TIME)
           write(6,*) ' GFlops/s for fact      ', dla_thread_stats(th,MULTILAP_FACT_FLOPS)
        end do

        write(6,*) '!--------------------------------------------------------------------!'
        write(6,*) ' Outer iterate ', k
        write(6,*) ' Time for analysis      ', maxval(dla_thread_stats(:,MULTILAP_ANALYZE_TIME))
        write(6,*) ' Time for factorization ', maxval(dla_thread_stats(:,MULTILAP_FACT_TIME))
        write(6,*) ' GFlops/s for fact      ', sum(dla_thread_stats(:,MULTILAP_FACT_FLOPS)) / params%nb_fact_omp
     end if

     dla_final_stats(MULTILAP_ANALYZE_TIME) = dla_final_stats(MULTILAP_ANALYZE_TIME) + &
          & maxval(dla_thread_stats(:,MULTILAP_ANALYZE_TIME))
     dla_final_stats(MULTILAP_FACT_TIME) = dla_final_stats(MULTILAP_FACT_TIME) + &
          & maxval(dla_thread_stats(:,MULTILAP_FACT_TIME))
     dla_final_stats(MULTILAP_FACT_FLOPS) = dla_final_stats(MULTILAP_FACT_FLOPS) + &
          & (sum(dla_thread_stats(:,MULTILAP_FACT_FLOPS)) / params%nb_fact_omp)
     !
     deallocate(dla_thread_stats)

     !
     ! Distribute solves among parallel threads
     !
     allocate(dla_thread_stats( params%nb_solv_omp, MULTILAP_MAX_STAT ))
     dla_thread_stats(:,:) = 0.0

     write(6,*) '!--------------------------------------------------------------------!'
     write(6,*) ' Nb of OpenMP threads enrolled for solution      = ', params%nb_solv_omp
     write(6,*) ' Nb of pthreads used in PaStiX for solution      = ', params%nb_solv_thd

     call split_parall( params%nb_sys, params%nb_solv_omp )

     do i = 1, params%nb_solve

        write(6,*) '!--------------------------------------------------------------------!'
        write(6,*) ' Solve iteration nr. ', i

        !
        !$OMP PARALLEL NUM_THREADS(params%nb_solv_omp) DEFAULT(NONE) &
        !$OMP SHARED(params, sys_array, k, i)         &
        !$OMP SHARED(ila_thrsz, ila_thrmn, ila_thrmx) &
        !$OMP SHARED(dla_thread_stats, iparm, dparm)  &
        !$OMP PRIVATE(th, im, ir, j )                 &
        !$OMP PRIVATE(matrix, rhs, x_ptr, b_ptr)      &
        !$OMP PRIVATE(info)

        !$OMP DO SCHEDULE(STATIC,1)
        solve_loop: do th = 1, params%nb_solv_omp
           solve_loop2: do j = ila_thrmn(th), ila_thrmx(th)

              if ( params%multirhs ) then
                 im = j
                 ir = 1
              else
                 im = ((j-1) / params%nb_fact_thd) + 1
                 ir = mod( j-1, params%nb_fact_thd ) + 1
              end if
              matrix => sys_array(im)
              rhs    => matrix%rhs( ir )

              if ( .not. (rhs%nrhs .gt. 0)) then
                 cycle
              endif

              !   2- The right hand side
              x_ptr = c_loc(rhs%x)
              b_ptr = c_loc(rhs%b)

              if (i==1) then
                 call spmGenRHS( SpmRhsRndX, rhs%nrhs, matrix%spm, &
                      & c_null_ptr, params%n, b_ptr, params%n, info )

                 if ( params%output ) then
                    ! Backup initial b that will be overwritten by check
                    rhs%b0(:) = rhs%b(:)
                 end if
              else
                 if ( params%output ) then
                    ! Restore the initial b overwritten by check
                    rhs%b(:) = rhs%b0(:)
                 end if
              end if

              ! Set the sequential scheduler to enable multiple solves in parallel with the same matrix
              if ( .not. params%multirhs ) then
                 matrix%iparm(IPARM_SCHEDULER) = PastixSchedSequential
              end if

              ! 3- Permute the b pointer
              call pastix_subtask_applyorder( matrix%pastix_data, SpmComplex64, PastixDirForward, &
                   &                          params%n, rhs%nrhs, b_ptr, params%n, info )

              rhs%x(:) = rhs%b(:)

              ! 4- Solve the problem
              call pastix_subtask_solve( matrix%pastix_data, rhs%nrhs, &
                   & x_ptr, matrix%spm%n, info )

              dla_thread_stats(th, MULTILAP_SOLV_TIME) =  &
                   & dla_thread_stats(th,MULTILAP_SOLV_TIME) + matrix%dparm(DPARM_SOLV_TIME)

              ! 5- Refine the solution
              call pastix_subtask_refine(            &
                   & matrix%pastix_data,             &
                   & matrix%spm%n, rhs%nrhs,         &
                   & b_ptr, matrix%spm%n,            &
                   & x_ptr, matrix%spm%n, info )

              ! 6- Apply the backward permutation on b and x
              call pastix_subtask_applyorder( matrix%pastix_data, SpmComplex64, PastixDirBackward, &
                   &                          params%n, rhs%nrhs, b_ptr, params%n, info )

              call pastix_subtask_applyorder( matrix%pastix_data, SpmComplex64, PastixDirBackward, &
                   &                          params%n, rhs%nrhs, x_ptr, params%n, info )

           end do solve_loop2
        end do solve_loop
        !$OMP END DO
        !$OMP END PARALLEL

        ! Restore the initial scheduler
        if ( .not. params%multirhs ) then
           do im = 1, params%nb_mat
              matrix => sys_array(im)
              matrix%iparm(IPARM_SCHEDULER) = iparm(IPARM_SCHEDULER)
           end do
        end if

        if (params%output) then
           do th = 1, params%nb_solv_omp
              do j = ila_thrmn(th), ila_thrmx(th)

                 if ( params%multirhs ) then
                    im = j
                    ir = 1
                 else
                    im = ((j-1) / params%nb_fact_thd) + 1
                    ir = mod( j-1, params%nb_fact_thd ) + 1
                 end if
                 matrix => sys_array(im)
                 rhs    => matrix%rhs( ir )

                 if ( .not. (rhs%nrhs .gt. 0)) then
                    cycle
                 endif

                 write(6,*) '!--------------------------------------------------------------------!'
                 write(6,*) ' Check results for system ', im, ir

                 !
                 ! Check the solution
                 !
                 call spmPrintInfo( matrix%spm )

                 b_ptr = c_loc(rhs%b)
                 x_ptr = c_loc(rhs%x)

                 call spmCheckAxb( matrix%dparm(DPARM_EPSILON_REFINEMENT), rhs%nrhs, &
                      & matrix%spm,    &
                      & c_null_ptr, params%n, &
                      & b_ptr,      params%n, &
                      & x_ptr,      params%n, info )

                 ginfo = ginfo + info
              end do
           end do

           write(6,*) '!--------------------------------------------------------------------!'
           do th = 1, params%nb_fact_omp
              write(6,*) ' Thread ', th
              write(6,*) ' Time for solution      ', dla_thread_stats(th,MULTILAP_SOLV_TIME)
           end do
        end if
     end do ! End of iteration loop

     if (params%output) then
        write(6,*) '!--------------------------------------------------------------------!'
        write(6,*) ' Outer iterate ', k
        write(6,*) ' Time for solution      ', maxval(dla_thread_stats(:,MULTILAP_SOLV_TIME))
     end if

     dla_final_stats(MULTILAP_SOLV_TIME) = dla_final_stats(MULTILAP_SOLV_TIME) + &
          & maxval(dla_thread_stats(:,MULTILAP_SOLV_TIME))

     deallocate(dla_thread_stats)

     ! Free memory
     do im=1, params%nb_mat
        matrix => sys_array(im)
        ! 6- Destroy the C data structure
        call spmExit( matrix%spm )
        deallocate( matrix%spm )
     end do
  end do main_loop

  ! Destroy the matrices structures
  do im = 1, params%nb_mat
     !
     matrix => sys_array(im)

     call pastixFinalize( matrix%pastix_data )

     do ir = 1, matrix%nsys
        rhs => matrix%rhs(ir)

        if ( rhs%nrhs .gt. 0 ) then
           if ( params%output ) then
              ! Allocate a backup of b that will be destroyed by the check
              deallocate(rhs%b0)
           end if
           deallocate(rhs%b)
           deallocate(rhs%x)
        end if
     end do
     deallocate( matrix%rhs )
  end do

  deallocate(sys_array)
  if (allocated(ila_thrsz)) deallocate(ila_thrsz)
  if (allocated(ila_thrmn)) deallocate(ila_thrmn)
  if (allocated(ila_thrmx)) deallocate(ila_thrmx)

  write(6,*) '!====================================================================!'

  write(6,*) ' Overall final statistics'
  write(6,*) ' Time for analysis      ', dla_final_stats(MULTILAP_ANALYZE_TIME)
  write(6,*) ' Time for factorization ', dla_final_stats(MULTILAP_FACT_TIME)
  write(6,*) ' GFlops/s for fact      ', dla_final_stats(MULTILAP_FACT_FLOPS) / params%nb_outit
  write(6,*) ' Time for solution      ', dla_final_stats(MULTILAP_SOLV_TIME)
  write(6,*) '!====================================================================!'

  call exit(ginfo)

8888 format( '  ', A6, ' ', A6, ' ', A6, ' ', A6 )
8889 format( '  ', I6, ' ', I6, ' ', I6, ' ', I6 )

contains

  !
  ! @brief Initialize the testing parameters based on the input configuration file
  !
  ! @param[output] params
  !
  subroutine multilap_init( params )
    type(multilap_param), intent(out), target :: params
    integer                                   :: val, dim1, dim2, dim3

    ! Read a dummy line.
    read( 5, * )

    ! Read if we enable/disable the verbose mode
    read( 5, * ) val
    params%output = .not. ( val .eq. 0 )

    ! Read if we enable/disable the check and correct
    read( 5, * ) val
    params%checkmat = .not. ( val .eq. 0 )

    ! Read if we enable/disable the multi-threaded solve
    read( 5, * ) val
    params%multirhs = .not. ( val .eq. 0 )

    ! Read the total number of threads
    read( 5, * ) params%nb_thrd

    if( params%nb_thrd < 1 ) then
       write( 6, fmt = 9999 ) 'nb_thrd', params%nb_thrd, params%nb_thrd, 1
       call exit(1)
    end if

    ! Read the number of pastix instances
    read( 5, * ) params%nb_fact_omp

    if( params%nb_fact_omp .lt. 1 ) then
       write( 6, fmt = 9999 ) 'nb_fact_omp', params%nb_fact_omp, 1
       call exit(1)
    else if( params%nb_fact_omp .gt. params%nb_thrd) then
       write( 6, fmt = 9998 ) 'nb_fact_omp', params%nb_fact_omp, params%nb_thrd
       call exit(1)
    endif

    params%nb_fact_thd = params%nb_thrd / params%nb_fact_omp
    if ( (params%nb_fact_thd * params%nb_fact_omp) .ne. params%nb_thrd ) then
       write( 6, * ) 'nb_thrd (', params%nb_thrd, ')must be a multiple of nb_pastix (', params%nb_fact_omp, ')'
       call exit(1)
    end if

    ! Read the number of outer most iterations
    read( 5, * ) params%nb_outit

    if( params%nb_outit .lt. 1 ) then
       write( 6, fmt = 9999 ) 'nb_outit', params%nb_outit, 1
       call exit(1)
    end if

    ! Read the number of matrices to factorize
    read( 5, * ) params%nb_mat

    if( params%nb_mat .lt. params%nb_fact_omp ) then
       write( 6, fmt = 9999 ) 'nb_mat', params%nb_mat, params%nb_fact_omp
       call exit(1)
    endif

    ! Read the number of RHS per iteration
    read( 5, * ) params%nb_rhs

    if( params%nb_rhs .lt. 1 ) then
       write( 6, fmt = 9999 ) 'nb_rhs', params%nb_rhs, 1
       call exit(1)
    end if

    ! Read the number of solve per iteration
    read( 5, * ) params%nb_solve

    if( params%nb_solve .lt. 1 ) then
       write( 6, fmt = 9999 ) 'nb_solve', params%nb_solve, 1
       call exit(1)
    end if

    ! Read the dimensions of the laplacian
    read( 5, * ) dim1

    if( dim1 .lt. 1 ) then
       write( 6, fmt = 9999 ) 'dim1', dim1, 1
       call exit(1)
    end if

    read( 5, * ) dim2

    if( dim2 .lt. 1 ) then
       write( 6, fmt = 9999 ) 'dim2', dim2, 1
       call exit(1)
    end if

    read( 5, * ) dim3

    if( dim3 .lt. 1 ) then
       write( 6, fmt = 9999 ) 'dim3', dim3, 1
       call exit(1)
    end if

    if ( params%multirhs ) then
       params%nb_sys      = params%nb_mat
       params%nb_solv_omp = params%nb_fact_omp
       params%nb_solv_thd = params%nb_fact_thd
    else
       params%nb_sys      = params%nb_mat * params%nb_fact_thd
       params%nb_solv_omp = params%nb_thrd
       params%nb_solv_thd = 1
    endif

    params%dim1 = dim1
    params%dim2 = dim2
    params%dim3 = dim3
    params%n    = dim1 * dim2 * dim3
    params%nnz  = (2*(dim1)-1) * dim2 * dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1)

    write(6,*) '!--------------------------------------------------------------------!'
    write(6,*) '!           Multiple Laplacian testing configuration                 !'
    write(6,*) '!--------------------------------------------------------------------!'
    write(6,*) ' Nb of threads                               = ', params%nb_thrd
    write(6,*) ' Nb of PaStiX instances                      = ', params%nb_fact_omp
    write(6,*) ' Nb of outer iterations                      = ', params%nb_outit
    write(6,*) ' Nb of distinct matrices                     = ', params%nb_mat
    write(6,*) ' Nb of RHS to solve per matrix               = ', params%nb_rhs
    write(6,*) ' Nb of solve phase to perform per RHS        = ', params%nb_solve
    write(6,*) ' Size of each matrix                         = ', params%n, &
         &     '(', dim1, ' x ', dim2, ' x ', dim3, ')'
    write(6,*) ' Nbr of non zero entries per matrix          = ', params%nnz

    if ( params%multirhs ) then
       write(6,*) ' The multirhs mode is enabled'
    else
       write(6,*) ' The multirhs mode is disabled'
    endif

9998 format( ' Invalid input value: ', A8, '=', I6, '; must be <=', I6 )
9999 format( ' Invalid input value: ', A8, '=', I6, '; must be >=', I6 )

  end subroutine multilap_init

  !
  ! Generate a single laplacian spm matrix with the parameters stores
  ! in the matrix description
  !
  subroutine multilap_genOneMatrix( params, matrix, ib_out, ib )
    type(multilap_param),       intent(in),    target :: params
    type(sys_lin),              intent(inout), target :: matrix
    integer(kind=pastix_int_t), intent(in)            :: ib, ib_out
    integer(kind=pastix_int_t)                        :: dim1, dim2, dim3
    integer(kind=pastix_int_t)                        :: i, j, k, l
    type(spmatrix_t), pointer                         :: spm2

    !
    ! Laplacian dimensions
    !
    dim1 = params%dim1
    dim2 = params%dim2
    dim3 = params%dim3

    !
    ! Create the spm out of the internal data
    !
    allocate( matrix%spm )
    call spmInit( matrix%spm )
    matrix%spm%mtxtype = SpmHermitian
    matrix%spm%flttype = SpmComplex64
    matrix%spm%fmttype = SpmIJV
    matrix%spm%n       = params%n
    matrix%spm%nnz     = params%nnz
    matrix%spm%dof     = 1

    call spmUpdateComputedFields( matrix%spm )
    call spmAlloc( matrix%spm )

    call c_f_pointer( matrix%spm%rowptr, matrix%rowptr, [params%nnz] )
    call c_f_pointer( matrix%spm%colptr, matrix%colptr, [params%nnz] )
    call c_f_pointer( matrix%spm%values, matrix%values, [params%nnz] )

    l = 1
    do i=1,dim1
       do j=1,dim2
          do k=1,dim3
             matrix%rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
             matrix%colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
             matrix%values(l) = (6., 0.)

             if (i == 1) then
                matrix%values(l) = matrix%values(l) - (1., 0.)
             end if
             if (i == dim1) then
                matrix%values(l) = matrix%values(l) - (1., 0.)
             end if
             if (j == 1) then
                matrix%values(l) = matrix%values(l) - (1., 0.)
             end if
             if (j == dim2) then
                matrix%values(l) = matrix%values(l) - (1., 0.)
             end if
             if (k == 1) then
                matrix%values(l) = matrix%values(l) - (1., 0.)
             end if
             if (k == dim3) then
                matrix%values(l) = matrix%values(l) - (1., 0.)
             end if

             matrix%values(l) = matrix%values(l) * 8.
             l = l + 1

             if (i < dim1) then
                matrix%rowptr(l) =  i    + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                matrix%colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                matrix%values(l) = -(1., 0.)
                l = l + 1
             end if
             if (j < dim2) then
                matrix%rowptr(l) = (i-1) + dim1 *  j    + dim1 * dim2 * (k-1) + 1
                matrix%colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                matrix%values(l) = -(1., 0.)
                l = l + 1
             end if
             if (k < dim3) then
                matrix%rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 *  k    + 1
                matrix%colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                matrix%values(l) = -(1., 0.)
                l = l + 1
             end if
          end do
       end do
    end do

    matrix%values(:) = matrix%values(:)*(dble(ib) * dble(ib_out) / 4.0)

    if (params%checkmat) then
       allocate( spm2 )
       call spmCheckAndCorrect( matrix%spm, spm2, info )
       if ( info .ne. 0 ) then
          call spmExit( matrix%spm )
          matrix%spm = spm2
       end if
       deallocate( spm2 )
    else
       call spmConvert( SpmCSC, matrix%spm, info )
    endif
  end subroutine multilap_genOneMatrix

  !
  subroutine split_parall(id_size, id_thr)

    integer, intent(in) :: id_size
    integer, intent(in) :: id_thr

    if (allocated(ila_thrsz)) deallocate(ila_thrsz)
    allocate(ila_thrsz(id_thr))
    if (allocated(ila_thrmn)) deallocate(ila_thrmn)
    allocate(ila_thrmn(id_thr))
    if (allocated(ila_thrmx)) deallocate(ila_thrmx)
    allocate(ila_thrmx(id_thr))

    ila_thrsz(:) = id_size / id_thr
    ila_thrsz(1:id_size-id_thr*ila_thrsz(1)) =&
         &  ila_thrsz(1:id_size-id_thr*ila_thrsz(1)) + 1

    ila_thrmn(1) = 1
    do i = 1, id_thr-1
       ila_thrmx(i)   = ila_thrmn(i) + ila_thrsz(i)-1
       ila_thrmn(i+1) = ila_thrmx(i) + 1
    end do
    ila_thrmx(id_thr) = id_size

  end subroutine split_parall

end program fmultilap


