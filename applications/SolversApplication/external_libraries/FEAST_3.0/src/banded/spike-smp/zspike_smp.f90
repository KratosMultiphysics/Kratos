
  !=========================================================================================
  !Copyright (c) 2009-2015, The Regents of the University of Massachusetts, Amherst.
  !E. Polizzi research lab
  !All rights reserved.
  !
  !Redistribution and use in source and binary forms, with or without modification, 
  !are permitted provided that the following conditions are met:
  !
  !1. Redistributions of source code must retain the above copyright notice, this list of conditions 
  !   and the following disclaimer.
  !2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
  !   and the following disclaimer in the documentation and/or other materials provided with the distribution.
  !3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
  !    products derived from this software without specific prior written permission.
  !
  !THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
  !BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
  !ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
  !EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
  !SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
  !LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
  !IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  !==========================================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! SPIKE-SMP - COMPLEX DOUBLE-PRECISIONS routines   !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! List of routines:

!!!!!! Documented !!!!!!!!

! zspike_gbtrf
! zspike_gbtrs
! zspike_gbtrsi
! zspike_gbsv
! zspike_tune


!!!! Auxiliary !!!!!!!!!

! zspike_gbtrf2
! zspike_gbtrs2
! zspike_gbtrst2
! zspike_itrefinement
! zspike_matmul
! zspike_multi
! zspike_invred
! zspike_prep_recn
! zspike_solve_recn
! zspike_solve_recn_transpose
! zspike_multi_transpose

! zspike_GBTRFk2
! zspike_GBTRSk2
! zspike_simple_solve
! zspike_simple_solve_transpose
! zspike_vector_norm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine zspike_gbtrf(spm,n,kl,ku,A,lda,work,info)
! Purpose
! -------
! This subroutine performs the SPIKE DS factorization
! NOTE : This subroutine naturally destroys A, so if you 
! might want to perform iterative refinement later, or
! use A for something else, you should save it before
! calling this subroutine.
! -------
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! A (in/out) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! work (in/out)
! This is a large work array, which contains parts of the reduced system which must be communicated to the solve stage
! The required size is (klu^2)*4*(#levels*#partitions)+(#partitions-1))
! The value 4*(#levels*#partitions)+(#partitions-1)) is placed in spm(10) by spikeinit.
! So, you may instead allocate work to be (klu*klu)*spm(10)
! 
! info (out)
! Information parameter
! info = 0 -> The banded LU factorization on each partition did not require boosting
! info = 1 -> The banded LU factorization on some partition required boosting. The answer will be approximate. Iterative refinement may be necessary.
! info = 2 -> There was an illegal value in your A matrix.
  !=====================================================================
  !  Braegan Spring - Eric Polizzi - 2015
  !=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer, dimension(64) :: spm_default
integer :: n,kl,ku,lda,klu,info,i
complex(kind=(kind(1.0d0))), dimension(*) :: work
complex(kind=(kind(1.0d0))), dimension(*) :: A
integer, dimension(:), pointer :: Ajmin
integer :: infoloc
integer(8),parameter :: fout=6
logical :: test


klu=max(kl,ku)

if (spm(1)==1) then
  call wwrite_n(fout)
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_s(fout, '*********** SPIKE-SMP -BEGIN ******************')
  call wwrite_n(fout) 
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout)
!!!!!!!!!!!! Print the Spike parameters which has been changed from default
  call wwrite_s(fout, 'List of input parameters spm(1:64)-- if different from default')
  call wwrite_n(fout)

  call spikeinit_default(spm_default)
  do i=1,19
      if ((i/=10).and.(spm(i)/=spm_default(i))) then
      call wwrite_s(fout, '   spm(')
      call wwrite_i(fout, i)
      call wwrite_s(fout, ')=')
      call wwrite_i(fout, spm(i))
      call wwrite_n(fout)
    endif
  enddo
call wwrite_s(fout, 'Size system')  
     call wwrite_t(fout) 
     call wwrite_i(fout,N)
     call wwrite_n(fout)
     call wwrite_s(fout, 'kl,ku      ')  
     call wwrite_t(fout) 
     call wwrite_i(fout,kl)
     call wwrite_t(fout) 
     call wwrite_i(fout,ku)
     call wwrite_n(fout)
     call wwrite_s(fout, '#Threads (Total) available           ')  
     call wwrite_i(fout,spm(25))
     call wwrite_n(fout)
     if (spm(25)/=spm(22)) then
      call wwrite_s(fout, '**ATTENTION** --- SPIKE cannot use all the threads for this problem') 
      call wwrite_n(fout)
      call wwrite_s(fout, '#Threads effectively used by SPIKE   ') 
      call wwrite_i(fout,spm(22))
      call wwrite_n(fout)
     end if
  call wwrite_i(fout, spm(20) )
  call wwrite_s(fout, " Partitions and " )
  call wwrite_i(fout, spm(22) )
  call wwrite_s(fout, " Threads" )
  call wwrite_n(fout)
  call wwrite_s(fout, "Partition ratios: " )
  call wwrite_d(fout, spm(4) * .1d0)
  call wwrite_t(fout)
  call wwrite_d(fout, spm(5) * .1d0)
  call wwrite_n(fout)
  call wwrite_s(fout, "#Partitions with 1 and 2 threads" )
  call wwrite_n(fout)
  call wwrite_i(fout, spm(24) )
  call wwrite_t(fout)
  call wwrite_i(fout, spm(23) )
  call wwrite_n(fout)

  call wallocate_1i(Ajmin,spm(20)+1,infoloc)
  call spikerl_calc_size_partitions(spm,Ajmin,n)

  call wwrite_s(fout, "Partition sizes: " )
  do i=1,spm(20)
    call wwrite_i(fout, Ajmin(i+1) - Ajmin(i))
    call wwrite_t(fout)
  enddo
  call wwrite_n(fout)

end if

call zspike_gbtrf2(spm,n,kl,ku,A,lda,work(1),work(2*klu*klu*spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),info)

if (spm(1)==1) call wdeallocate_1i(Ajmin)


end subroutine zspike_gbtrf




subroutine zspike_gbtrs(spm,trans,n,kl,ku,nrhs,A,lda,work,B,ldb)
! Purpose
! -------
! This subroutine performs the SPIKE solve
! -------
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The number of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! nrhs (in)
! The number of columns in B
! 
! A (in) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
!  
! work (in)
! This is a large work array, which contains parts of the reduced system which must be communicated to the solve stage
! The required size is (klu^2)*4*(#levels*#partitions)+(#partitions-1))
! The value 4*(#levels*#partitions)+(#partitions-1)) is placed in spm(10) by spikeinit.
! So, you may instead allocate work to be (klu*klu)*spm(10)
!
! B (in/out)
! The collection of vectors on which the solve is performed
!
! ldb (in)
! The leading dimension of B. Usually this is n, but if you feel clever you can try other things. 
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: n,kl,ku,lda,klu,ldb,nrhs
complex(kind=(kind(1.0d0))), dimension(*) :: work
complex(kind=(kind(1.0d0))), dimension(lda,*) :: A
complex(kind=(kind(1.0d0))), dimension(ldb,*) :: B
character :: trans
integer(8),parameter :: fout=6

klu=max(kl,ku)

if((trans .eq. 'n') .or. (trans .eq. 'N')) then
  call zspike_gbtrs2(spm,n,kl,ku,lda,nrhs,A,work(1),work(2*klu*klu* spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),B,ldb)
else
 ! call zspike_gbtrst2(spm,trans,n,kl,ku,lda,nrhs,A,work(1),work(2*klu*klu* spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),B,ldb)
  if((trans .eq. 't') .or. (trans .eq. 'T')) then
    call zspike_gbtrst2(spm,'T',n,kl,ku,lda,nrhs,A,work(1),work(2*klu*klu* spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),B,ldb)
  else
    call zspike_gbtrst2(spm,'C',n,kl,ku,lda,nrhs,A,work(1),work(2*klu*klu* spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),B,ldb)
  endif
endif


if((spm(1) .eq. 1) .and. (spm(31) .eq. 0)) then
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_s(fout, '************SPIKE-SMP END *********************')
  call wwrite_n(fout) 
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_n(fout)
endif

end subroutine zspike_gbtrs



subroutine zspike_gbtrsi(spm,trans,n,kl,ku,nrhs,C,ldc,A,lda,work,B,ldb)
! Purpose
! -------
! This subroutine performs the SPIKE solve with iterative refinement
! -------
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The number of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! nrhs (in)
! The number of columns in B
! 
! C (in)
! A copy of the origional version of A, used for iterative refinement
! 
! ldc (in) 
! The leading dimensions of C. ldc=kl+ku+1
!
! A (in) 
! The banded matrix on which the SPIKE factorization was performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
!  
! work (in)
! This is a large work array, which contains parts of the reduced system which must be communicated to the solve stage
! The required size is (klu^2)*4*(#levels*#partitions)+(#partitions-1))
! The value 4*(#levels*#partitions)+(#partitions-1)) is placed in spm(10) by spikeinit.
! So, you may instead allocate work to be (klu*klu)*spm(10)
!
! B (in/out)
! The collection of vectors on which the solve is performed
!
! ldb (in)
! The leading dimension of B. Usually this is n, but if you feel clever you can try other things. 
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: n,kl,ku,lda,klu,ldb,nrhs,ldc
complex(kind=(kind(1.0d0))), dimension(*) :: work
complex(kind=(kind(1.0d0))), dimension(lda,*) :: A
complex(kind=(kind(1.0d0))), dimension(ldb,*) :: B
complex(kind=(kind(1.0d0))), dimension(ldc,*) :: C
character :: trans
complex(kind=(kind(1.0d0))), dimension(:,:), pointer :: oB
integer :: i,p
integer :: infoloc
integer(8),parameter :: fout=6


call wallocate_2z(oB,ldb,nrhs,infoloc)
call ZLACPY( 'F',n,nrhs,B, N, oB, N )

klu=max(kl,ku)

!Spike solve
spm(31) = 1
call zspike_gbtrs(spm,trans,n,kl,ku,nrhs,A,lda,work,B,ldb)

!Applying iterative refinement
call zspike_itrefinement(spm,trans,n,kl,ku,nrhs,A,work,lda,C,ldc,B,oB,ldb)

call wdeallocate_2z(oB)
spm(31) = 0
if((spm(1) .eq. 1) .and. (spm(31) .eq. 0)) then
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_s(fout, '************SPIKE-SMP END *********************')
  call wwrite_n(fout) 
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_n(fout)
endif


end subroutine zspike_gbtrsi




subroutine zspike_gbsv(spm,N,KL,KU,NRHS,A,LDA,B,LDB,INFO)
! Purpose
! -------
! This is the do-it-all factorize and solve subroutine. 
! A matrix and a collection of vectors are entered, and the 
! problem  A^(-1)B => B is solved. 
! The residual value for B is checked, and simple iterative 
! refinement is performed if it is above the given tolerance.
!
! Overall, this is similar to the lapack Xgbsv, with the exception that
! the matrix A is returned to the initial state upon return.
! 
! Optimization note: 
! If the optimization flag spm(2)=2 is set, this subroutine 
! will attempt to find appropriate ratios for the partition sizes. 
! These ratios are dependent on the characteristics of the matrix used 
! (specificly the relationship between klu and nrhs) and the hardware 
! on which SPIKE is run. The hardware characteristics are embodied
! in the tuning constant K=spm(7).
! This subroutine may attempt to find K, but finding K requires a 
! large single threaded factorization and solve. To avoid this cost,
! zspike_gbsv will check if there is a non-zero value in spm(7), which 
! indicates that the constant has already been found.
! Because this constant should only depend on the the hardware, 
! one may manually set the tuning constant by setting spm(7) after calling
! spikeinit. 
! (zspike_tune will help with this process)
! This need not be run every time zspike_gbsv is called, or even every
! time the calling program is run (although communicating the constant
! between program runs is an exercise left to the reader [I don't want
! to muck up your environment variables]).  
! 
! -------
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! nrhs (in)
! The number of columns in B
! 
! A (in) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! B (in/out) 
! The set of vectors on which the matrix solve will be performed; solution on exit
! 
! ldb (in)
! The leading dimensions of B
! 
! info (out)
! Information parameter
! info = 0 -> The banded LU factorization on each partition did not require boosting
! info = 1 -> The banded LU factorization on some partition required boosting. The answer will be approximate. Iterative refinement may be necessary.
! info = 2 -> There was an illegal value in your A matrix.
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer :: n,kl,ku,klu,nrhs
integer :: ldc,lda,ldb,info,i,p
complex(kind=(kind(1.0d0))), dimension(lda,*) :: A
complex(kind=(kind(1.0d0))), dimension(ldb,*) :: B
complex(kind=(kind(1.0d0))),  dimension(:,:),pointer :: C
integer,  dimension(:),pointer :: ajmin
integer, dimension(64) :: spm
character :: trans
complex(kind=(kind(1.0d0))),  dimension(:),pointer :: work
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
double precision :: K,klud,nrhsd
integer :: infoloc
integer(8),parameter :: fout=6
logical :: test

spm(31) = 1 !Indicate that factorize and solve are done in one call.
klu = max(kl,ku)

trans='n'
klud = 1.0d0*klu 
nrhsd = 1.0d0*nrhs
if((spm(2) .eq. 2)) then
  K = 0.1d0*spm(7)
  spm(4) = int(10*(klud/(klud+K*nrhsd)*.5d0 + (      klud +       nrhsd)/(klud/K + nrhsd)))
  spm(5) = int(10*(klud/(klud+K*nrhsd)      + (1.5d0*klud + 2.0d0*nrhsd)/(klud/K + nrhsd)))
endif
ldc=kl+ku+1

!Give the system a hint that A and  C should be 'owned' by the various processors by copying them in a parallel section
call wallocate_2z(C,ldc,n,infoloc)
call wallocate_1i(Ajmin,spm(20)+1,infoloc)
call wallocate_1z(work,klu*klu*spm(10),infoloc)

call spikerl_calc_size_partitions(spm,Ajmin,n)
!$OMP PARALLEL DO default(shared) private(p)
do i=1,spm(20)
  p=i
  call ZLACPY('A',ldc,Ajmin(p+1)-Ajmin(p),A(1,Ajmin(p)),lda,C(1,Ajmin(p)),ldc)
enddo

call zspike_gbtrf(spm, n,kl,ku, C,lda,work, info)
call zspike_gbtrsi(spm,trans,n,kl,ku,nrhs,A,lda,C,lda,work,B,ldb)

call wdeallocate_2z(C)
call wdeallocate_1i(Ajmin)
call wdeallocate_1z(work)

end subroutine zspike_gbsv 





subroutine zspike_tune(spm)
! Purpose
! -------
! Auto-tuning- Calculate the partition size ratios, for COMPLE DOUBLE PRECISION values
! These ratios are described by a tuning constant, which is dependent on the Big-O behavior of the system banded factorization and solve.
! Discovering this constant requires a large, single threaded fatorize and solve.
! As a result, this function is fairly slow. 
! -------
!
! Arguments  
! ---------
! spm (out)
! An array of 64 integers, which control various features of the spike factorization and solve.
! These parameters are described in the description of spikeinit
! For this function we only modify the following values:
!
! spm(4)  = out. This describes the ratio between the size of the large first and last partitions, and the medium 2-thread spike partitions. It is given as 10X the ratio
! spm(5)  = out. This describes the ratio between the size of the large first and last partitions, and the small single-thread LU partitions. It is given as 10X the ratio
! spm(7)  = out. This describes the tuning constant from which these ratios are derived. It is given as 10X the tuning constant. 
! The derived relationship between these values should be:
! spm(4) = spm(7) + 5
! spm(5) = spm(7)*1.5 + 10 
! ---------
!=====================================================================
!  Braegan Spring -  2015
!=====================================================================
implicit none 
include 'f90_noruntime_interface.fi'
integer n,klu,lda,info
double precision :: t1,t2,t0,tim
integer, dimension(64) :: spm
complex(kind=(kind(1.0d0))), dimension(:,:),pointer :: A
complex(kind=(kind(1.0d0))), dimension(:,:),pointer :: B
real(kind=(kind(1.0d0))) :: nzero, norma
double precision :: K
complex(kind=((kind(1.0d0)))), parameter :: z_one=(1.0d0,0.0d0)
integer :: infoloc
integer(8),parameter :: fout=6
double precision, external :: OMP_GET_WTIME

n=10000
klu=200
nzero=1.0d-13

norma=4*klu*1.0d0
lda = 2*klu+1 

call wallocate_2z(A,lda,n,infoloc)
call wallocate_2z(B,n,klu,infoloc)

B=z_one
A=z_one
A(klu+1,:) = 3*klu*z_one

t0 = OMP_GET_WTIME()
call ZGBALU(n,klu,klu,A,lda,nzero,norma,info)
t1 = OMP_GET_WTIME()
call ZTBSM('L','N','U',n,klu,klu,A,lda,B,n)
t2 = OMP_GET_WTIME()

K = 2.0d0*(t2-t1)/(t1-t0)
spm(7) = int(K*10)
spm(4) = spm(7) +5
spm(5) = spm(7) +spm(7)/2 +10

if(spm(1) .eq. 1) then
  call wwrite_s(fout,'Large-bandwith partition ratios found:')
  call wwrite_d(fout, spm(4)*.1d0)
  call wwrite_t(fout)
  call wwrite_d(fout, spm(5)*.1d0)
  call wwrite_n(fout)
endif

call wdeallocate_2z(A)
call wdeallocate_2z(B)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE: Most of the functions below here aren't really intended to be called by an end user.  !
!                          They are not documented                                            !
!                                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









subroutine zspike_gbtrf2(spm,n,kl,ku,A,lda,rV,rW,red,info)
! Purpose
! -------
! This is the subroutine that ultimately does the work of the
! SPIKE factorization. However, there is a more user friendly 
! interface. The subroutine zspike_gbtrf should be used instead. 
! 
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! A (in/out) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! rV, rW, red, (out) 
! Work space, to contain various parts of the reduced system
! Must be communicated to the solve subroutine
!
! info (out)
! Information parameter
! info = 0 -> The banded LU factorization on each partition did not require boosting
! info = 1 -> The banded LU factorization on some partition required boosting. The answer will be approximate. Iterative refinement may be necessary.
! info = 2 -> There was an illegal value in your A matrix.
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm 
logical ::invred,test
integer :: t1,t2,t3,t4,t5,tim
integer :: n,kl,ku,lda,nbpart,nblevel,nthread
integer :: i,j,k,l,klu,p,p_adjust,myid,info,max_info
integer :: naux,m,nj,sizeA1,counter
integer,external :: omp_get_thread_num
integer, dimension(:), pointer :: Ajmin
integer, dimension(:,:), pointer :: keys
complex(kind=(kind(1.0d0))), dimension(lda,*):: A
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),2*max(kl,ku),*) :: red
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),max(kl,ku),spm(21),*)  :: rV,rW
!complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),spm(23)*2*max(kl,ku)) :: spike2_red
!complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),spm(23)*(kl+ku)) :: spike2_grj
complex(kind=(kind(1.0d0))), pointer, dimension(:,:) :: spike2_red
complex(kind=(kind(1.0d0))), pointer, dimension(:,:) :: spike2_grj
real(kind=(kind(1.0d0))) :: nzero, norma
complex(kind=(kind(1.0d0))),dimension(:,:),pointer::vw
complex(kind=(kind(1.0d0))),dimension(:,:),pointer::vw2
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
double precision, dimension(:),pointer :: timing_array
character(len=100) :: timing_print_format
integer :: infoloc
integer(8),parameter :: fout=6
double precision, external :: OMP_GET_WTIME
double precision :: start_time,end_time
double precision :: redsys_start_time, redsys_end_time
double precision :: parallel_start_time, parallel_end_time
character(len=1) :: norm
real(kind=(kind(1.0d0))), dimension(:), pointer :: normwork
real(kind=(kind(1.0d0))), external :: zlangb

call wallocate_2z(spike2_red, 2*max(kl,ku),spm(23)*2*max(kl,ku),infoloc)
call wallocate_2z(spike2_grj, 2*max(kl,ku),spm(23)*(kl+ku),infoloc)

start_time = OMP_GET_WTIME()

if(spm(14) .eq. 0) then
  !Norm infinity
  norm='i'
endif

if(spm(14) .eq. 1) then
  !Norm one
  norm='o'
endif

if(spm(14) .eq. 2) then
  !Norm Frobenius
  norm='f'
endif

call wallocate_1d(timing_array,spm(22),infoloc)

nzero=1.0d-13
!If we've only been given one thread, do the banded primitive and get out of here.
if(spm(20) == 1) then

  parallel_start_time = OMP_GET_WTIME()
  timing_array(1) = OMP_GET_WTIME()

  ! Lapack norm compute requires a work array for infinorm
  if(spm(14) .eq. 0) then
    call wallocate_1d(normwork,n,infoloc)
  else
    !not actually used but passing an unallocated array seems to be a problem 
    call wallocate_1d(normwork,1,infoloc)
  endif
  norma = zlangb( norm, n, kl, ku, A, lda, normwork )
  call wdeallocate_1d(normwork)  
  call ZGBALU(n,kl,ku,A,lda,nzero,norma,info)

  timing_array(1) = OMP_GET_WTIME() - timing_array(1)
  parallel_end_time = OMP_GET_WTIME()
  redsys_start_time = 0.0d0
  redsys_end_time   = 0.0d0
else

call wallocate_1i(Ajmin,spm(20)+1,infoloc)
call spikerl_calc_size_partitions(spm,Ajmin,n)

!Ajmin holds the positions of the sub-arrays. 
nthread = spm(22)
nbpart = spm(20)
nblevel = spm(21)
invred=.false.
klu=max(kl,ku)

call wallocate_2i(keys,4,spm(23),infoloc)
keys(1,1:spm(23)) = 0
!keys(2:4,1:spm(23)) = -1

! the keys is a shared section of memory that the simplified two-threaded spike partitions will use to communicate 
! when they have finished various steps of their work. 
do i=1,spm(23)
  keys(2,i) = 2*(i)-1
  keys(3,i) = 2*(i)
enddo

if(spm(23) > 0) then
  call wallocate_2z(vw2,Ajmin(3)    -Ajmin(2)       ,(ku+kl)*spm(23),infoloc) 
endif 
if(spm(20)-spm(23) > 0) then
 ! call wallocate_2z(vw,Ajmin(nbpart)-Ajmin(nbpart-1),(ku+kl)*spm(20),infoloc) 
  call wallocate_2z(vw,Ajmin(nbpart)-Ajmin(nbpart-1),(ku+kl)*(spm(24)-2),infoloc) 
endif

info = 0
max_info=0
parallel_start_time=OMP_GET_WTIME()
!$omp parallel do default(none) private(i,j,p,p_adjust,tim,counter,info,norma,infoloc,normwork) shared(nblevel,vw,vw2,spike2_red,spike2_grj,keys,A,rV,rW,spm,ajmin,nthread,max_info,timing_array) firstprivate(nbpart,lda,nzero,klu,kl,ku,norm)
do i=1,nthread,1

  timing_array(i) = OMP_GET_WTIME()
  call threadmap(i,p,nbpart,spm(23))
 
! <<<<<<<<<<< LU Factorize A(1) - A(nbpart-1) and Solving Ai.Vi=bi, A2..nbpart-1.W2..nbpart-1=c2..nbpart-1 >>>>>>>>>>> !

  ! Lapack norm compute requires a work array for infinorm
  if(spm(14) .eq. 0) then
    call wallocate_1d(normwork,Ajmin(p+1)-Ajmin(p),infoloc)
  else
    !not actually used but passing an unallocated array seems to be a problem 
    call wallocate_1d(normwork,1,infoloc)
  endif
  norma = zlangb( norm, Ajmin(p+1)-Ajmin(p), kl, ku, A(1,Ajmin(p)), lda, normwork )
  call wdeallocate_1d(normwork)  

  if(p == 1) then
    rV(:,:,1,p) = z_zero
    call ZLACPY( 'L',ku , ku, A(1,Ajmin(p+1)), lda-1,rV(2*klu-ku+1,1,1,p), 2*klu)

    call ZGBALU(Ajmin(p+1)-Ajmin(p),kl,ku,A(1,Ajmin(p)),lda,nzero,norma,info)
    call ZTBSM('L','N','U', klu, ku, kl, A(ku+1,Ajmin(p+1)-klu),lda, rV(klu+1,1,1,p), 2*klu)
! I need the L^-1b for the transpose case (so save it in an unused bit of rV) 
! This section is unused because there is no top of the V spike for the first partiton
    call ZLACPY( 'L',klu , klu, rV(2*klu-ku+1,1,1,p), 2*klu, rV(1,1,1,p), 2*klu)
    call ZTBSM('U','N','N', klu, ku, ku, A(1,   Ajmin(p+1)-klu),lda, rV(klu+1,1,1,p), 2*klu)
  endif  

  if(p==nbpart) then
    rW(:,:,1,p) = z_zero
    call ZLACPY( 'U', kl, kl, A(kl+ku+1,Ajmin(p)-kl), lda-1, rW(1,klu-kl+1,1,p), 2*klu)

    call ZGBAUL(Ajmin(p+1)-Ajmin(p),kl,ku,A(1,Ajmin(p)),lda,nzero,norma,info)
    call ZTBSM('U','N','U', klu, kl, ku,A(1,Ajmin(p)),lda,rW(1,klu-kl+1,1,p),2*klu)
!I need the U^-1 for the transpose case. 
    call ZLACPY( 'U', klu, klu, rW(1,klu-kl+1,1,p), 2*klu, rW(klu+1,klu-kl+1,1,p), 2*klu)
    call ZTBSM('L','N','N', klu, kl, kl,A(ku+1,Ajmin(p)),lda,rW(1,klu-kl+1,1,p),2*klu)
  endif

! Two thread partitions 
! Here we will call the simplified two-thread spike to work on the partition.
  if(p > 1 .and. p <= 1 + spm(23)) then 
!copy b and c in to vw

!Copying in to the area of vw worked on by the first thread.
    if(keys(2,p-1) == omp_get_thread_num()) then
      vw2(1:(Ajmin(p+1)-Ajmin(p))/2,(p-2)*(kl+ku)+1:(p-1)*(kl+ku))=z_zero
      call ZLACPY('U', kl, kl, A(kl+ku+1,Ajmin(p)-kl), lda-1, vw2(1,(p-2)*(kl+ku)+ku+1), Ajmin(p+1)-Ajmin(p))
    endif

!Copying in to the area of vw worked on by the second thread.
    if(keys(3,p-1) == omp_get_thread_num()) then
      vw2((Ajmin(p+1)-Ajmin(p))/2+1:Ajmin(p+1)-Ajmin(p),(p-2)*(kl+ku)+1:(p-1)*(kl+ku))=z_zero
      call ZLACPY('L', ku, ku, A(1,Ajmin(p+1))       , kl+ku, vw2(Ajmin(p+1)-Ajmin(p)-ku+1,(p-2)*(kl+ku)+1), Ajmin(p+1)-Ajmin(p))
    endif

    call zspike_GBTRFk2(Ajmin(p+1)-Ajmin(p),kl,ku,A(1,Ajmin(p)),lda,info,keys(1,p-1))
    call zspike_GBTRSk2('N',Ajmin(p+1)-Ajmin(p),kl+ku,kl,ku,A(1,Ajmin(p)),lda,vw2(1,(p-2)*(kl+ku)+1),Ajmin(p+1)-Ajmin(p),keys(1,p-1), spike2_red(1,(p-2)*2*klu+1), spike2_grj(1,(p-2)*(kl+ku)+1),.true.)
    
!copy vw out to rV
    if(keys(2,p-1) == omp_get_thread_num()) then
      rV(1:klu,:,1,p) = z_zero
      rW(1:klu,:,1,p) = z_zero
!Initializing the space into which the v tips are copied to zero. 

      !Copying in to the area of vw worked on by the first thread... 
      call ZLACPY('F', klu, ku, vw2(1,(p-2)*(kl+ku)+1)                           , Ajmin(p+1)-Ajmin(p), rV(1,1,1,p)            ,2*klu)
      call ZLACPY('F', klu, kl, vw2(1,(p-2)*(kl+ku)+ku+1)                        , Ajmin(p+1)-Ajmin(p), rW(1,klu-kl+1,1,p)     ,2*klu)
    endif 
    if(keys(3,p-1) == omp_get_thread_num()) then
      rV(klu+1:2*klu,:,1,p) = z_zero
      rW(klu+1:2*klu,:,1,p) = z_zero
      !Copying in to the area of vw worked on by the second thread... 
      call ZLACPY('F', klu, ku, vw2(Ajmin(p+1)-Ajmin(p)-klu+1,(p-2)*(kl+ku)+1)   , Ajmin(p+1)-Ajmin(p), rV(klu+1,1,1,p)        ,2*klu)
      call ZLACPY('F', klu, kl, vw2(Ajmin(p+1)-Ajmin(p)-klu+1,(p-2)*(kl+ku)+ku+1), Ajmin(p+1)-Ajmin(p), rW(klu+1,klu-kl+1,1,p) ,2*klu)
    endif

  endif

  if(p > 1 + spm(23) .and. p < nbpart) then 
    p_adjust = p - (spm(23))-2
    vw(:, p_adjust*(kl+ku)+1:(p_adjust+1)*(kl+ku)) = z_zero

    call ZLACPY( 'L', ku, ku, A(1,Ajmin(p+1)),        lda-1, vw(Ajmin(p+1)-Ajmin(p)-ku+1,p_adjust*(kl+ku)+1), Ajmin(p+1)-Ajmin(p))
    call ZLACPY( 'U', kl, kl, A(kl+ku+1,Ajmin(p)-kl), lda-1, vw(1,p_adjust*(kl+ku)+ku+1),                     Ajmin(p+1)-Ajmin(p))

    call ZGBALU(Ajmin(p+1)-Ajmin(p),kl,ku,A(1,Ajmin(p)),lda,nzero,norma,info)

    call ZTBSM('L','N','U', klu, ku, kl,A(ku+1,Ajmin(p+1)-klu),lda, vw(Ajmin(p+1)-Ajmin(p)-klu+1,p_adjust*(kl+ku)+1),Ajmin(p+1)-Ajmin(p))
    call ZTBSM('L','N','U',Ajmin(p+1)-Ajmin(p), kl, kl,A(ku+1,Ajmin(p)),lda,vw(1,p_adjust*(kl+ku)+ku+1),Ajmin(p+1)-Ajmin(p))
    call ZTBSM('U','N','N',Ajmin(p+1)-Ajmin(p), ku+kl, ku,A(1,Ajmin(p)),lda,vw(1,p_adjust*(kl+ku)+1),Ajmin(p+1)-Ajmin(p))


    rV(:,:,1,p) = z_zero
    rW(:,:,1,p) = z_zero

    call ZLACPY('F', klu, ku, vw(1,p_adjust*(kl+ku)+1)                           , Ajmin(p+1)-Ajmin(p), rV(1,1,1,p),2*klu)
    call ZLACPY('F', klu, kl, vw(1,p_adjust*(kl+ku)+ku+1)                        , Ajmin(p+1)-Ajmin(p), rW(1,klu-kl+1,1,p),2*klu)
    call ZLACPY('F', klu, ku, vw(Ajmin(p+1)-Ajmin(p)-klu+1,p_adjust*(kl+ku)+1)   , Ajmin(p+1)-Ajmin(p), rV(klu+1,1,1,p),2*klu)
    call ZLACPY('F', klu, kl, vw(Ajmin(p+1)-Ajmin(p)-klu+1,p_adjust*(kl+ku)+ku+1), Ajmin(p+1)-Ajmin(p), rW(klu+1,klu-kl+1,1,p),2*klu)
  endif

! Handle info
! We want errors to take presedence (as they are most critical -- unrecoverable)

  if(info .eq. 0) then
    !DGBALU was happy -- don't have to do anything
  else
    if(info .gt. 0) then
    !DGBALU had to boost 
      info = 1
    else
    !DGBALU found an error/illegal value 
      info = 2
    endif
  endif

  !$omp atomic
  max_info=max(info,max_info)
  timing_array(i) = OMP_GET_WTIME() - timing_array(i) 
enddo 
parallel_end_time=OMP_GET_WTIME()

info=max_info

if(spm(23) > 0) then
  call wdeallocate_2z(vw2) 
endif 
if(spm(20)-spm(23) > 0) then
  call wdeallocate_2z(vw)
endif
!<<<<<<<<<<<<<< Reduced systems initialization >>>>>>>>>>>>>>>>>>!

redsys_start_time = OMP_GET_WTIME()

! Reduced systems initialization - first level !
!$omp parallel do default(shared) private(p)
do p=1,nbpart-1,1
  call zspike_invred(invred,klu,rV(1,1,1,p),rW(1,1,1,p+1),red(1,1,p))
enddo

call zspike_prep_recn(invred,kl,ku,rV(1,1,1,1),rW(1,1,1,1),red(1,1,1),nbpart,nblevel)

redsys_end_time = OMP_GET_WTIME()

call wdeallocate_2i(keys)
call wdeallocate_1i(Ajmin)
call wdeallocate_2z(spike2_red)
call wdeallocate_2z(spike2_grj)

endif

end_time = OMP_GET_WTIME()

if(spm(1) .eq. 1) then
  call wwrite_s(fout,'---------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Factorization Time|')
  call wwrite_n(fout)
  call wwrite_s(fout,'------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|   LU or SPIKE on blocks    |')
  call wwrite_n(fout)
  call wwrite_s(fout,'--------------------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'| Thread |   Partition   |          Time         |')
  call wwrite_n(fout)
  call wwrite_s(fout,'--------------------------------------------------')
  call wwrite_n(fout)
  do i=1,spm(22)
    call threadmap(i,p,nbpart,spm(23))
    call wwrite_s(fout,'|    ')
    call wwrite_i(fout,i)
    call wwrite_s(fout,'   |       ')
    call wwrite_i(fout,p)
    call wwrite_s(fout,'       | ')
    call wwrite_d(fout,timing_array(i))
    call wwrite_s(fout,' |')
    call wwrite_n(fout)
  enddo 
  call wwrite_s(fout,'--------------------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Blocks Factorize |       ')
  call wwrite_d(fout,parallel_end_time-parallel_start_time)
  call wwrite_s(fout,'  |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Reduced System   |       ')
  call wwrite_d(fout,redsys_end_time-redsys_start_time)
  call wwrite_s(fout,'  |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Overall Factorize|       ')
  call wwrite_d(fout,end_time-start_time)
  call wwrite_s(fout,'  |')
  call wwrite_n(fout)
  call wwrite_s(fout,'--------------------------------------------------')
  call wwrite_n(fout)
  endif


call wdeallocate_1d(timing_array)


end subroutine zspike_gbtrf2






subroutine zspike_gbtrs2(spm,n,kl,ku,lda,nrhs,A,rV,rW,red,f,ldf)
! Purpose
! -------
! This is the subroutine that ultimately does the work of the
! SPIKE non-transpose solve. However, there is a more user friendly 
! interface. The subroutine zspike_gbtrs should be used instead. 
! 
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! A (in/out) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! rV, rW, red, (out) 
! Work space, to contain various parts of the reduced system
! Came from the factorization subroutine
! 
! f (in/out)
! The collection of vectors B.
!
! ldf (in)
! The leading dimension of f.
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: n,kl,ku,lda,nrhs, nbpart,nblevel,nthread,ldf
integer :: nj,sizeA1,naux,m,mysize,mystart
integer, dimension(:), pointer :: Ajmin
complex(kind=(kind(1.0d0))), dimension(lda,*):: A
complex(kind=(kind(1.0d0))), dimension(ldf,*):: f
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),2*max(kl,ku),*):: red
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),max(kl,ku),spm(21),*):: rV,rW
integer, dimension(:,:), pointer :: keys
integer ::i,j,klu,p,myid,info,k
integer,external :: omp_get_thread_num
!complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),spm(23)*2*max(kl,ku)) :: spike2_red
!complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),spm(23)*nrhs) :: spike2_grj
complex(kind=(kind(1.0d0))), pointer, dimension(:,:) :: spike2_red
complex(kind=(kind(1.0d0))), pointer, dimension(:,:) :: spike2_grj
complex(kind=(kind(1.0d0))), dimension(:,:), pointer:: g
complex(kind=(kind(1.0d0))), dimension(:,:), pointer :: g2
complex(kind=(kind(1.0d0))), dimension(:,:,:), pointer :: grj
logical :: invred
integer :: t1,t2,t3,t4,t5,tim
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
double precision, dimension(:,:),pointer :: timing_array
character(len=100) :: timing_print_format
integer :: infoloc
integer(8),parameter :: fout=6
double precision, external :: OMP_GET_WTIME
double precision :: start_time,end_time
double precision :: redsys_start_time, redsys_end_time
double precision :: parallel_start_time_1, parallel_end_time_1
double precision :: parallel_start_time_2, parallel_end_time_2

call wallocate_2z(spike2_red,2*max(kl,ku),spm(23)*2*max(kl,ku),infoloc)
call wallocate_2z(spike2_grj,2*max(kl,ku),spm(23)*nrhs,infoloc)

start_time = OMP_GET_WTIME()
call wallocate_2d(timing_array,spm(22),2,infoloc)

if(spm(20) == 1) then
  parallel_start_time_1 = OMP_GET_WTIME()
  timing_array(1,1) = OMP_GET_WTIME()
  call ZTBSM('L','N','U', n, nrhs, kl, A(ku+1,1),lda,f,ldf)
  timing_array(1,1) = OMP_GET_WTIME() - timing_array(1,1)
  parallel_end_time_1 = OMP_GET_WTIME()

  parallel_start_time_2 = OMP_GET_WTIME()
  timing_array(1,2) = OMP_GET_WTIME()
  call ZTBSM('U','N','N', n, nrhs, ku, A        ,lda,f,ldf)
  timing_array(1,2) = OMP_GET_WTIME() - timing_array(1,2)
  parallel_end_time_2 = OMP_GET_WTIME()
  redsys_start_time = 0.0d0
  redsys_end_time = 0.0d0
else

call wallocate_1i(Ajmin,spm(20)+1,infoloc)
call spikerl_calc_size_partitions(spm,Ajmin,n)
nthread = spm(22)
nbpart = spm(20) 
nblevel = spm(21)
invred=.false.

klu=max(kl,ku)
call wallocate_3z(grj,2*klu,nrhs,nbpart,infoloc) !! modified rhs
grj=z_zero

call wallocate_2i(keys,4,nbpart,infoloc)
keys(1,1:nbpart) = 0
do i=1,spm(23)
  keys(2,i) = 2*(i)-1
  keys(3,i) = 2*(i)
enddo

if(nbpart-(spm(23)+2) > 0) then
  call wallocate_2z(g,Ajmin(nbpart)-Ajmin(nbpart-1),nrhs*(nbpart-(spm(23)+2)),infoloc)
endif

if(spm(23) > 0) then
  call wallocate_2z(g2,Ajmin(3)-Ajmin(2),nrhs*spm(23),infoloc)
endif


parallel_start_time_1 = OMP_GET_WTIME()

!$omp parallel do firstprivate(t1,t2,i,p,tim,mystart,mysize,nrhs,kl,ku,lda,ldf,klu,nbpart,nthread), default(none), shared(timing_array,spm,Ajmin,A,f,g2,g,spike2_red,spike2_grj,grj,keys)
do i=1,nthread,1
  timing_array(i,1) = OMP_GET_WTIME()
  call threadmap(i,p,nbpart,spm(23))
!! Solving Ajgj=fj
  if (p<=nbpart-1) then
    if (p==1) then
      call ZTBSM('L','N','U', Ajmin(p+1)-Ajmin(p), nrhs,kl,A(ku+1,Ajmin(p)),lda,f(Ajmin(p),1),ldf)
      call ZLACPY('F', klu, nrhs, f(Ajmin(p+1)-klu,1),ldf, grj(klu+1,1,p), 2*klu )
      call ZTBSM('U','N','N', klu, nrhs, ku,A(1,Ajmin(p+1)-klu),lda,grj(klu+1,1,p),2*klu)
      grj(1:klu,1:nrhs,p)=z_zero
    else
      if(p <= spm(23)+1) then
    ! Spike version case
        mysize = Ajmin(3) - Ajmin(2)
        mystart = nrhs*(p-2)+1

        if(keys(2,p-1) == omp_get_thread_num()) then
          call ZLACPY('F',mysize/2,nrhs,f(Ajmin(p),1),ldf,g2(1,mystart),mysize)
        endif

        if(keys(3,p-1) == omp_get_thread_num()) then
          call ZLACPY('F',mysize/2,nrhs,f(Ajmin(p)+mysize/2,1),ldf,g2(1+mysize/2,mystart),mysize)
        endif
        call zspike_GBTRSk2('N',mysize,nrhs,kl,ku,A(1,Ajmin(p)),lda,g2(1,mystart),mysize,keys(1,p-1), spike2_red(1,(p-2)*2*klu+1),spike2_grj(1,(p-2)*(nrhs)+1),.false.)

        if(keys(2,p-1) == omp_get_thread_num()) then
          call ZLACPY('F',klu,nrhs,g2(1,mystart),mysize,grj(1,1,p),2*klu)
        endif
        if(keys(3,p-1) == omp_get_thread_num()) then
          call ZLACPY('F',klu,nrhs,g2(mysize-klu+1,mystart),mysize,grj(klu+1,1,p),2*klu)
        endif
      else 
        mysize = Ajmin(nbpart)-Ajmin(nbpart-1)
        mystart = nrhs*(p-(2+spm(23)))+1
        call ZLACPY('F',mysize,nrhs,f(Ajmin(p),1),ldf,g(1, mystart),mysize)
        call ZTBSM('L','N','U', mysize, nrhs, kl,A(ku+1,Ajmin(p)),lda,g(1, mystart),mysize)
        call ZTBSM('U','N','N', mysize, nrhs, ku,A(1,Ajmin(p)),lda,   g(1, mystart),mysize)

        call ZLACPY('F',klu,nrhs,g(1, mystart), mysize,grj(1,1,p),2*klu)
        call ZLACPY('F',klu,nrhs,g(mysize-klu+1, mystart),mysize,grj(klu+1,1,p),2*klu)

      endif
    endif
  endif
!! solving Anbpart gnbpart =f nbpart
  if (p==nbpart) then
    call ZTBSM('U','N','U', Ajmin(p+1)-Ajmin(p), nrhs, ku,A(1,Ajmin(p)),lda,f(Ajmin(p),1),ldf)
    call ZLACPY('F', klu, nrhs, f(Ajmin(p),1),ldf, grj(1,1,p), 2*klu )
    call ZTBSM('L','N','N', klu, nrhs, kl,A(ku+1,Ajmin(p)),lda,grj(1,1,p),2*klu)
    grj(klu+1:,1:nrhs,p)=z_zero
  endif
  timing_array(i,1) = OMP_GET_WTIME() - timing_array(i,1)
enddo
parallel_end_time_1 = OMP_GET_WTIME()

if(spm(23) > 0) then
  call wdeallocate_2z(g2)
endif
if(nbpart-(2+spm(23)) > 0) then
  call wdeallocate_2z(g)
endif

redsys_start_time = OMP_GET_WTIME()
call zspike_solve_recn(invred,kl,ku,nrhs,rV(1,1,1,1),rW(1,1,1,1),red(1,1,1),nbpart,nblevel,grj(1,1,1))
!<<<<<<<<<<<<<<<<<<<<<<<< Retrieval  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

!$omp parallel do
do p=1,nbpart,1
!!! Aixi=fi-Bi*x(i+1)b-Ci*x(i-1)t
  if(p==1) then
    call ZTRMM('L','L','N','N',ku,nrhs,-z_one,A(1,Ajmin(p+1)),lda-1,grj(1,1,p+1),2*klu)
    call ZTBSM('L','N','U',ku, nrhs,kl,A(ku+1,Ajmin(p+1)-ku),lda,grj(1,1,p+1),2*klu) ! finish U
  endif

  if(p==nbpart) then
    call ZTRMM('L','U','N','N',kl,nrhs,-z_one,A(kl+ku+1,Ajmin(p)-kl),lda-1,grj(2*klu-kl+1,1,p-1),2*klu)
    call ZTBSM('U','N','U', kl, nrhs, ku,A(1,Ajmin(p)),lda,grj(2*klu-kl+1,1,p-1),2*klu) ! finish L
  endif

  if (p<=nbpart-1 .and. p>=2) then
    call ZTRMM('L','L','N','N',ku,nrhs,-z_one,A(1,Ajmin(p+1)),lda-1,grj(1,1,p+1),2*klu)
    call ZTRMM('L','U','N','N',kl,nrhs,-z_one,A(kl+ku+1,Ajmin(p)-kl),lda-1,grj(2*klu-kl+1,1,p-1),2*klu)
  end if
enddo

!$omp parallel do
do p=1,nbpart,1
  if (p>=2) then
    do i=1,nrhs
      call ZAXPY(kl,z_one,grj(2*klu-kl+1,i,p-1),1,f(Ajmin(p),i),1) 
    enddo
  endif

  if (p<=nbpart-1) then
    do i=1,nrhs
      call ZAXPY(ku,z_one,grj(1,i,p+1),1,f(Ajmin(p+1)-ku,i),1) 
    enddo
  endif
enddo

redsys_end_time = OMP_GET_WTIME() 

!allocate(keys(4,nbpart))
keys(1,1:nbpart) = 0
do i=1,spm(23)
  keys(2,i) = 2*(i)-1
  keys(3,i) = 2*(i)
enddo

!! Solving Ajgj=fj
parallel_start_time_2 = OMP_GET_WTIME()
!$omp parallel do  private(t1,t2,tim,p), shared(spike2_red,spike2_grj,keys)
do i=1,nthread,1
  call threadmap(i,p,nbpart,spm(23))
  timing_array(i,2) = OMP_GET_WTIME()
  if(p==1) then
    call ZTBSM('U','N','N', Ajmin(p+1)-Ajmin(p), nrhs, ku,A(1,Ajmin(p)),lda,f(Ajmin(p),1),ldf)
  endif

  if(p==nbpart) then
    call ZTBSM('L','N','N', Ajmin(p+1)-Ajmin(p), nrhs, kl,A(ku+1,Ajmin(p)),lda,f(Ajmin(p),1),ldf)
  endif

  if(p > 1 .and. p <= 1 + spm(23)) then 
    call zspike_GBTRSk2('N',Ajmin(p+1)-Ajmin(p), nrhs, kl, ku, A(1,Ajmin(p)), lda, f(Ajmin(p),1), ldf, keys(1,p-1), spike2_red(1,(p-2)*2*klu+1),spike2_grj(1,(p-2)*nrhs+1),.false.)
  endif

  if(p > 1 + spm(23) .and. p < nbpart) then 
    call ZTBSM('L','N','U', Ajmin(p+1)-Ajmin(p), nrhs, kl, A(ku+1,Ajmin(p)), lda, f(Ajmin(p),1), ldf)
    call ZTBSM('U','N','N', Ajmin(p+1)-Ajmin(p) ,nrhs, ku, A(1,Ajmin(p)), lda, f(Ajmin(p),1), ldf)
  endif

!! solving Anbpart gnbpart =f nbpart
  timing_array(i,2) = OMP_GET_WTIME() - timing_array(i,2)
enddo
parallel_end_time_2 = OMP_GET_WTIME()
call wdeallocate_3z(grj)
call wdeallocate_2i(keys)
call wdeallocate_1i(Ajmin)
call wdeallocate_2z(spike2_red)
call wdeallocate_2z(spike2_grj)


endif

end_time = OMP_GET_WTIME()

if(spm(1) .eq. 1) then
  call wwrite_s(fout,'------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Solve Time|')
  call wwrite_n(fout)
  call wwrite_s(fout,'------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|  Solve Sweeps on blocks    |----------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'| Thread | Partition |         Time 1        |          Time 2       |')
  call wwrite_n(fout)
  do i=1,spm(22)
    call threadmap(i,p,nbpart,spm(23))
    call wwrite_s(fout,'|    ')
    call wwrite_i(fout,i)
    call wwrite_s(fout,'   |     ')
    call wwrite_i(fout,p)
    call wwrite_s(fout,'     | ')
    call wwrite_d(fout,timing_array(i,1))
    call wwrite_s(fout,' | ')
    call wwrite_d(fout,timing_array(i,2))
    call wwrite_s(fout,' |')
    call wwrite_n(fout)
  enddo 
  call wwrite_s(fout,'----------------------------------------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Blocks Solve        | ')
  call wwrite_d(fout,(parallel_end_time_1-parallel_start_time_1))
  call wwrite_s(fout,' | ')
  call wwrite_d(fout,(parallel_end_time_2-parallel_start_time_2))
  call wwrite_s(fout,' |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Reduced System Solve| ')
  call wwrite_d(fout,(redsys_end_time-redsys_start_time))
  call wwrite_s(fout,'                         |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Overall Solve       | ')
  call wwrite_d(fout,(end_time-start_time))
  call wwrite_s(fout,'                         |')
  call wwrite_n(fout)
  call wwrite_s(fout,'----------------------------------------------------------------------')
  call wwrite_n(fout)
endif

call wdeallocate_2d(timing_array)


end subroutine zspike_gbtrs2





subroutine zspike_gbtrst2(spm,trans,n,kl,ku,lda,nrhs,A,rV,rW,red,f,ldf)
! Purpose
! -------
! This is the subroutine that ultimately does the work of the
! SPIKE transpose solve. However, there is a more user friendly 
! interface. The subroutine zspike_gbtrs should be used instead. 
! 
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! A (in/out) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! rV, rW, red, (out) 
! Work space, to contain various parts of the reduced system
! Came from the factorization subroutine
! 
! f (in/out)
! The collection of vectors B.
!
! ldf (in)
! The leading dimension of f.
!=====================================================================
!  Braegan Spring - 2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: n,kl,ku,lda,nrhs, nbpart,nblevel,nthread
integer :: nj,sizeA1,naux,m,mysize,mystart
integer, dimension(:), pointer :: Ajmin
complex(kind=(kind(1.0d0))), dimension(lda,*):: A
complex(kind=(kind(1.0d0))), dimension(ldf,*):: f
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),2*max(kl,ku),*):: red
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),max(kl,ku),spm(21),*):: rV,rW
integer, dimension(:,:), pointer :: keys
integer ::i,j,klu,p,myid,info,k,ldf
double precision :: ratio
integer,external :: omp_get_thread_num
!complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),spm(23)*2*max(kl,ku)) :: spike2_red
!complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),spm(23)*nrhs) :: spike2_grj
complex(kind=(kind(1.0d0))), pointer, dimension(:,:) :: spike2_red
complex(kind=(kind(1.0d0))), pointer, dimension(:,:) :: spike2_grj
complex(kind=(kind(1.0d0))), dimension(:,:), pointer :: g
complex(kind=(kind(1.0d0))), dimension(:,:), pointer:: g2
complex(kind=(kind(1.0d0))), dimension(:,:,:), pointer :: grj
complex(kind=(kind(1.0d0))), dimension(:,:), pointer :: gaux
integer, dimension(2*max(kl,ku)*spm(20)) :: ipiv_red 
logical :: invred
integer :: t1,t2,t3,t4,t5,tim
double precision, dimension(:,:),pointer :: timing_array
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
integer :: infoloc
integer(8),parameter :: fout=6
double precision, external :: OMP_GET_WTIME
double precision :: start_time,end_time
double precision :: redsys_start_time, redsys_end_time
double precision :: parallel_start_time_1, parallel_end_time_1
double precision :: parallel_start_time_2, parallel_end_time_2
character :: trans

start_time = OMP_GET_WTIME()

call wallocate_2z(spike2_red,2*max(kl,ku),spm(23)*2*max(kl,ku),infoloc)
call wallocate_2z(spike2_grj,2*max(kl,ku),spm(23)*nrhs,infoloc)
call wallocate_2d(timing_array,spm(22),2,infoloc)

if(spm(20) == 1) then
  parallel_start_time_1 = OMP_GET_WTIME()
  timing_array(1,1) = OMP_GET_WTIME()
  call ZTBSM('U',trans,'N', n, nrhs, ku, A,lda,f,ldf)
  timing_array(1,1) = OMP_GET_WTIME() - timing_array(1,1)
  parallel_end_time_1 = OMP_GET_WTIME()

  parallel_start_time_2 = OMP_GET_WTIME()
  timing_array(1,2) = OMP_GET_WTIME()
  call ZTBSM('L',trans,'U', n, nrhs, kl, A(ku+1,1),lda,f,ldf)
  timing_array(1,2) = OMP_GET_WTIME() - timing_array(1,2)
  parallel_end_time_2 = OMP_GET_WTIME()
  redsys_start_time = 0.0d0
  redsys_end_time = 0.0d0
else

klu = max(kl,ku)

nbpart = spm(20)
nthread = spm(22)
nblevel = spm(21)
invred=.false.

call wallocate_1i(Ajmin,spm(20)+1,infoloc)
call wallocate_3z(grj,2*klu,nrhs,nbpart,infoloc) !! modified rhs
call wallocate_2z(gaux,2*klu,nrhs*nbpart,infoloc) !! modified rhs

call spikerl_calc_size_partitions(spm,Ajmin,n)
grj=z_zero
gaux=z_zero

!manually allocate 'keys' for the partitions on which we will use spike
call wallocate_2i(keys,4,nbpart,infoloc)
keys(1,1:nbpart) = 0
do i=1,spm(23)
  keys(2,i) = 2*(i)-1
  keys(3,i) = 2*(i)
enddo

parallel_start_time_1 = OMP_GET_WTIME()
!Set up the zero augmented right hand sides for the reduced system
!$omp parallel do default(shared) private(i,p,mysize,tim,t2,t1)
do i=1,nthread,1
  call threadmap(i,p,nbpart,spm(23))
  timing_array(i,1) = OMP_GET_WTIME()
  mysize = Ajmin(p+1)-Ajmin(p)
  if(p==1) then
    ! A1 is LU factorized
    !Augment the bottom of f1 with zero

    call ZLACPY('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, grj(klu+1,1,p), 2*klu)
    f(Ajmin(p+1)-klu:Ajmin(p+1)-1,1:nrhs) = z_zero

    ! Sweeps over the right hand side
    call ZTBSM('U',trans,'N',mysize,nrhs,ku,A(1,Ajmin(p)),   lda,f(Ajmin(p),1),ldf)

    ! Perform multiplication 
    ! We will need to unmodified values for f at a later point, so we pull out the ones that will need to work with
    call ZLACPY('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, gaux(klu+1,nrhs*(p-1)+1), 2*klu)
    call ZTRMM('L', 'L', trans,'N',ku,nrhs,z_one,rV(1,1,1,1),2*klu,gaux(klu+1+klu-ku,nrhs*(p-1)+1),2*klu)

  endif

  if(p==nbpart) then
   !Augment the top of the final part of f with zeroes

    call ZLACPY('F', klu, nrhs, f(Ajmin(p),1), ldf, grj(1,1,p), 2*klu)
    f(Ajmin(p):Ajmin(p)+klu-1,1:nrhs) = z_zero

    call ZTBSM('L',trans,'N',mysize,nrhs,kl,A(ku+1,Ajmin(p)),lda,f(Ajmin(p),1),ldf)

    !w's
    call ZLACPY('F', klu, nrhs, f(Ajmin(p),1), ldf, gaux(1,nrhs*(p-1)+1), 2*klu)
    call ZTRMM('L','U',trans,'N',kl,nrhs,z_one, rW(klu+1,klu-kl+1,1,p),2*klu,gaux(1,nrhs*(p-1)+1),2*klu)
  endif

  if(p > 1 + spm(23) .and. p < nbpart) then 

    !augment with 0
    call ZLACPY('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, grj(klu+1,1,p), 2*klu)
    f(Ajmin(p+1)-klu:Ajmin(p+1)-1,1:nrhs) = z_zero


    call ZLACPY('F', klu, nrhs, f(Ajmin(p),1), ldf, grj(1,1,p), 2*klu)
    f(Ajmin(p):Ajmin(p)+klu-1,1:nrhs) = z_zero

    ! Sweeps over the right hand side
    call ZTBSM('U',trans,'N',mysize,nrhs,ku,A(1,Ajmin(p)),   lda,f(Ajmin(p),1),ldf)
    call ZTBSM('L',trans,'U',mysize,nrhs,kl,A(ku+1,Ajmin(p)),lda,f(Ajmin(p),1),ldf)

    !v's
    call ZLACPY('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, gaux(klu+1,nrhs*(p-1)+1), 2*klu)
    call ZTRMM('L', 'L', trans,'N',ku,nrhs,z_one,A(1,Ajmin(p+1)),lda-1,gaux(klu+1+klu-ku,nrhs*(p-1)+1),2*klu)

    !w's    
    call ZLACPY('F', klu, nrhs, f(Ajmin(p),1), ldf, gaux(1,nrhs*(p-1)+1), 2*klu)
    call ZTRMM('L','U',trans,'N',kl,nrhs,z_one,A(lda,Ajmin(p)-kl),lda-1,gaux(1,nrhs*(p-1)+1),2*klu)

  endif

   
!Two thread partitions 
  if(p > 1 .and. p <= 1 + spm(23)) then 

    !augment with 0
    if(keys(3,p-1) == omp_get_thread_num()) then
      call ZLACPY('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, grj(klu+1,1,p), 2*klu)
      f(Ajmin(p+1)-klu:Ajmin(p+1)-1,1:nrhs) = z_zero
    endif
    if(keys(2,p-1) == omp_get_thread_num()) then
      call ZLACPY('F', klu, nrhs, f(Ajmin(p),1), ldf, grj(1,1,p), 2*klu)
      f(Ajmin(p):Ajmin(p)+klu-1,1:nrhs) = z_zero
    endif

    call zspike_GBTRSk2(trans,Ajmin(p+1)-Ajmin(p), nrhs, kl, ku, A(1,Ajmin(p)), lda, f(Ajmin(p),1), ldf, keys(1,p-1), spike2_red(1,(p-2)*2*klu+1),spike2_grj(1,(p-2)*nrhs+1))

    !v's
    if(keys(3,p-1) == omp_get_thread_num()) then
      call ZLACPY('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, gaux(klu+1,nrhs*(p-1)+1), 2*klu)
      call ZTRMM('L', 'L', trans,'N',ku,nrhs,z_one,A(1,Ajmin(p+1)),lda-1,gaux(klu+1+klu-ku,nrhs*(p-1)+1),2*klu)
    endif

    !w's    
    if(keys(2,p-1) == omp_get_thread_num()) then
      call ZLACPY('F', klu, nrhs, f(Ajmin(p),1), ldf, gaux(1,nrhs*(p-1)+1), 2*klu)
      call ZTRMM('L','U',trans,'N',kl,nrhs,z_one,A(lda,Ajmin(p)-kl),lda-1,gaux(1,nrhs*(p-1)+1),2*klu)
    endif

  endif


  timing_array(i,1) = OMP_GET_WTIME() - timing_array(i,1)
enddo
parallel_end_time_1 = OMP_GET_WTIME()

redsys_start_time = OMP_GET_WTIME()
!$omp parallel do default(shared) private(p,mysize,t1,t2,tim)
do p=1,nbpart,1
!  call threadmap(i,p,nbpart,spm(23))
! I know p-1+1 is silly, but I want to make it clear that this is inter-thread communication -- 
! the p-1 follows with what was used before, then the +1 stands for the next partition over
! So this works with the w spike from the next partition over
! so for grj, on the first index, 
! (p-1)*2*klu represents that we're working in partition p
! + klu indicates that we're writing to the bottom half
! + klu - ku indicates that we have implied 0's in gaux, because the dgemm that filled gaux
! began filling at index 1, rather than index ku
  if(p == 1) then
    do j=1,nrhs,1
      call ZAXPY(kl,-z_one,gaux(1,nrhs*(p-1+1)+j),1,grj(2*klu-kl+1,j,p),1)
    enddo 
  endif
 
  if(p == nbpart) then
! This works with the v spike from the previous partition
    do j=1,nrhs,1
      call ZAXPY(ku,-z_one,gaux(2*klu-ku+1,nrhs*(p-1-1)+j),1,grj(1,j,p),1)
    enddo
  endif 

  if(p > 1 .and. p < nbpart) then
    do j=1,nrhs,1
      call ZAXPY(kl,-z_one,gaux(1,nrhs*(p-1+1)+j),1,grj(2*klu-kl+1,j,p),1)
    enddo 
    do j=1,nrhs,1
      call ZAXPY(ku,-z_one,gaux(2*klu-ku+1,nrhs*(p-1-1)+j),1,grj(1,j,p),1)
    enddo
  endif
enddo

call wdeallocate_2z(gaux) !! modified rhs

!Reduced system
call zspike_solve_recn_transpose(trans,kl,ku,nrhs,rV(1,1,1,1),rW(1,1,1,1),red(1,1,1),nbpart,nblevel,grj(1,1,1)) 


redsys_end_time = OMP_GET_WTIME()

 
keys(1,1:nbpart) = 0
do i=1,spm(23)
  keys(2,i) = 2*(i)-1
  keys(3,i) = 2*(i)
enddo

! We did the sweeps A over f[middle] in the S stage, so we need only do the A sweeps on [yb;0;yt]
! That is, the vector we are solving for has the top and bottom tips equal to the things we got
! from the reduced system, and the middle all zeroes. 
if(nbpart-(spm(23)+2) > 0) then
  call wallocate_2z(g,Ajmin(nbpart)-Ajmin(nbpart-1),nrhs*(nbpart-(spm(23)+2)),infoloc)
endif

!g2 is the temp vector for the 2spike threads. It must be shared, because it is used to communicate between them.
if(spm(23) > 0) then
  call wallocate_2z(g2,Ajmin(3)-Ajmin(2),nrhs*spm(23),infoloc)
endif

parallel_start_time_2 = OMP_GET_WTIME()
!$omp parallel do default(shared) shared(g2,g) private(j,mystart,p,mysize,t1,t2,tim)
do i=1,nthread,1
  timing_array(i,2) = OMP_GET_WTIME()
  call threadmap(i,p,nbpart,spm(23))
  mysize = Ajmin(p+1)-Ajmin(p)
  if(p == 1) then
    call ZTBSM('U',trans,'N', klu, nrhs, ku, A(1,Ajmin(p+1)-klu), lda, grj(klu+1,1,p), 2*klu)

    do j=1,nrhs,1
      call ZAXPY(klu,z_one,grj(klu+1,j,p),1,f(Ajmin(p)+mysize-klu,j),1)
    enddo

    call ZTBSM('L',trans,'U', mysize, nrhs, kl, A(ku+1,Ajmin(p)), lda, f(1,1), ldf)
  endif

  if(p > 1 + spm(23) .and. p < nbpart) then 
    mystart = nrhs*(p-(2+spm(23)))+1
    g(:,mystart:mystart+nrhs-1) = z_zero
    call ZLACPY('F', klu, nrhs, grj(1,1,p), 2*klu, g(1,mystart), mysize)
    call ZLACPY('F', klu, nrhs, grj(klu+1,1,p), 2*klu, g(mysize-klu+1,mystart), mysize)
    call ZTBSM('U',trans,'N', mysize, nrhs, ku, A(1,Ajmin(p)),    lda, g(1,mystart), mysize)

!DTBSMPY combines the solve and sum, to save us trips to memory.
    call ZTBSMPY('L',trans,'U', mysize, nrhs, kl, A(ku+1,Ajmin(p)), lda, g(1,mystart), mysize, f(Ajmin(p),1), ldf)
! Keeping this around just in case...
!    call ZTBSM('L',trans,'U', mysize, nrhs, kl, A(ku+1,Ajmin(p)), lda, g(1,mystart), mysize)
!    do j=1,nrhs,1
!      call ZAXPY(mysize,1.0d0,g(1,mystart+j-1),1,f(Ajmin(p),j),1)
!    enddo
!    deallocate(g)
  endif

  if(p==nbpart) then 
    call ZTBSM('L',trans,'N',klu,nrhs,kl,A(ku+1,Ajmin(p)),lda,grj(1,1,p),2*klu)
    do j=1,nrhs,1
      call ZAXPY(klu,z_one,grj(1,j,p),1,f(Ajmin(p),j),1)
    enddo
    call ZTBSM('U',trans,'U',mysize,nrhs,ku,A(1,Ajmin(p)),   lda,f(Ajmin(p),1),ldf)
  endif

  if(p > 1 .and. p <= 1 + spm(23)) then 
    mystart = nrhs*(p-2)+1

    if(keys(2,p-1) == omp_get_thread_num()) then
      g2(1:mysize/2,mystart:mystart+nrhs-1) = z_zero
      call ZLACPY('F', klu, nrhs, grj(1,1,p), 2*klu, g2(1,mystart), mysize)
    endif

    if(keys(3,p-1) == omp_get_thread_num()) then
      g2(mysize/2+1:mysize,mystart:mystart+nrhs-1) = z_zero
      call ZLACPY('F', klu, nrhs, grj(klu+1,1,p), 2*klu, g2(mysize-klu+1,mystart), mysize)
    endif

    call zspike_GBTRSk2(trans,Ajmin(p+1)-Ajmin(p), nrhs, kl, ku, A(1,Ajmin(p)), lda, g2(1,mystart), mysize, keys(1,p-1), spike2_red(1,(p-2)*2*klu+1),spike2_grj(1,(p-2)*nrhs+1))

    if(keys(2,p-1) == omp_get_thread_num()) then
      do j=1,nrhs,1
        call ZAXPY(mysize/2,z_one,g2(1,mystart-1+j),1,f(Ajmin(p),j),1)
      enddo
    endif
    if(keys(3,p-1) == omp_get_thread_num()) then
      do j=1,nrhs,1
        call ZAXPY(mysize/2,z_one,g2(mysize/2+1,mystart-1+j),1,f(Ajmin(p)+mysize/2,j),1)
      enddo
    endif

  endif

  timing_array(i,2) = OMP_GET_WTIME() - timing_array(i,2)
enddo
parallel_end_time_2 = OMP_GET_WTIME()


if(spm(23) > 0) then
  call wdeallocate_2z(g2)
endif

if(nbpart-(spm(23)+2) > 0) then
  call wdeallocate_2z(g)
endif

call wdeallocate_2z(spike2_red)
call wdeallocate_2z(spike2_grj)
call wdeallocate_3z(grj)
call wdeallocate_2i(keys)
call wdeallocate_1i(Ajmin)
endif

end_time = OMP_GET_WTIME()

if(spm(1) .eq. 1) then
  call wwrite_s(fout,'------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Solve Time (Transpose)|')
  call wwrite_n(fout)
  call wwrite_s(fout,'------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|  Solve Sweeps on blocks    |----------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'| Thread | Partition |         Time 1        |          Time 2       |')
  call wwrite_n(fout)
  do i=1,spm(22)
    call threadmap(i,p,nbpart,spm(23))
    call wwrite_s(fout,'|    ')
    call wwrite_i(fout,i)
    call wwrite_s(fout,'   |     ')
    call wwrite_i(fout,p)
    call wwrite_s(fout,'     | ')
    call wwrite_d(fout,timing_array(i,1))
    call wwrite_s(fout,' | ')
    call wwrite_d(fout,timing_array(i,2))
    call wwrite_s(fout,' |')
    call wwrite_n(fout)
  enddo 
  call wwrite_s(fout,'----------------------------------------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Blocks Solve        | ')
  call wwrite_d(fout,(parallel_end_time_1-parallel_start_time_1))
  call wwrite_s(fout,' | ')
  call wwrite_d(fout,(parallel_end_time_2-parallel_start_time_2))
  call wwrite_s(fout,' |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Reduced System Solve| ')
  call wwrite_d(fout,(redsys_end_time-redsys_start_time))
  call wwrite_s(fout,'                         |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Overall Solve       | ')
  call wwrite_d(fout,(end_time-start_time))
  call wwrite_s(fout,'                         |')
  call wwrite_n(fout)
  call wwrite_s(fout,'----------------------------------------------------------------------')
  call wwrite_n(fout)
endif

call wdeallocate_2d(timing_array)

end subroutine zspike_gbtrst2



!Iterative refinment
subroutine zspike_itrefinement(spm,trans,n,kl,ku,nrhs,A,work,lda,C,ldc,B,oB,ldb)
!=====================================================================
!  Braegan Spring  2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: maxiter,n,kl,ku,lda,klu,ldb,nrhs,info,ldc
complex(kind=(kind(1.0d0))), dimension(lda,*) :: A
complex(kind=(kind(1.0d0))), dimension(ldb,*) :: B
complex(kind=(kind(1.0d0))), dimension(ldb,*) :: oB
complex(kind=(kind(1.0d0))), dimension(ldc,*) :: C
character :: trans
real(kind=(kind(1.0d0))),dimension(:),  pointer :: normB
real(kind=(kind(1.0d0))),dimension(:),  pointer :: normRes
integer,  dimension(:),pointer :: Ajmin
!double precision, dimension(ldb,nrhs) :: res
complex(kind=(kind(1.0d0))), dimension(:,:), pointer :: res
integer :: i,j,p
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
real(kind=(kind(1.0d0))), parameter :: d_zero_res=0.0d0
integer :: work_complete
complex(kind=(kind(1.0d0))), dimension(*) :: work
real(kind=(kind(1.0d0))) :: res_max
real(kind=(kind(1.0d0))), dimension(:),pointer  :: res_max_temp
integer :: infoloc
integer(8),parameter :: fout=6
integer :: norm
double precision, external :: OMP_GET_WTIME
double precision :: t1,t2


norm=spm(14)


j=0
work_complete = 0
maxiter = spm(11)

call wallocate_1i(Ajmin,spm(20)+1,infoloc)
call wallocate_2z(res,ldb,nrhs,infoloc)
call wallocate_1d(normB      ,nrhs          ,infoloc)
call wallocate_1d(normRes    ,nrhs          ,infoloc)

call spikerl_calc_size_partitions(spm,Ajmin,n)

!Each right hand side should have an individual norm. 

call zspike_vector_norm(spm,norm,oB,n,nrhs,normB)


call wallocate_1d(res_max_temp,spm(20),infoloc)

do while ((j<maxiter) .and. (work_complete .eq. 0)) 
  j=j+1

  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
    do i=1,nrhs
      res(Ajmin(p):Ajmin(p+1)-1,1:nrhs)=oB(Ajmin(p):Ajmin(p+1)-1,1:nrhs)
    enddo
  enddo
  res_max_temp = z_zero

!How far are we from the desired value? 
! Find RelativeResidual = norm(oB-AB)/norm(oB) 
! (since B contains the result) 

! Set res = oB-AB
t1 = OMP_GET_WTIME()
  call zspike_matmul(spm,Ajmin,trans,n,kl,ku,nrhs,C,lda,B,ldb,res)
t2 = OMP_GET_WTIME()
! RRES = norm(res)/norm(oB)

  call zspike_vector_norm(spm,norm,res,n,nrhs,normRes)
  res_max = d_zero_res
  do i=1,nrhs
    res_max = max(res_max,normRes(i)/normB(i))
  enddo

  if(spm(1) .eq. 1) then
    call wwrite_s(fout,'Num Refinements so far: ')
    call wwrite_i(fout,j-1)
    call wwrite_n(fout)
    call wwrite_s(fout,'Current relative residual norm: ')
    call wwrite_d(fout,res_max)
    call wwrite_n(fout)
  endif

  if(-log10(abs(res_max)) >= spm(12)) then 
   work_complete=1
    if(spm(1) .eq. 1) then
      call wwrite_s(fout,'Desired Residual Reached')
      call wwrite_n(fout)
    endif
  endif

 
  if(work_complete .eq. 0) then
    call zspike_gbtrs(spm,trans,n,kl,ku,nrhs,A,lda,work,res,ldb)
  
    !$OMP PARALLEL DO default(shared) private(i)
    do p=1,spm(20)
      do i=1,nrhs
        call ZAXPY(Ajmin(p+1)-Ajmin(p), z_one, res(Ajmin(p),i), 1, B(Ajmin(p),i), 1)
      enddo
    enddo
  endif

enddo
  call wdeallocate_1d(normB)
  call wdeallocate_1d(normRes) 
  call wdeallocate_1d(res_max_temp)
  call wdeallocate_2z(res)
  call wdeallocate_1i(Ajmin)

end subroutine zspike_itrefinement




subroutine zspike_matmul(spm,Ajmin,trans,n,kl,ku,nrhs,C,ldc,B,ldb,res)
!=====================================================================
!  Braegan Spring - 2015
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: n,kl,ku,lda,klu,ldb,nrhs,info,ldc
complex(kind=(kind(1.0d0))), dimension(ldb,*) :: B
complex(kind=(kind(1.0d0))), dimension(ldc,*) :: C
character :: trans
integer, dimension(spm(20)+1) :: Ajmin
complex(kind=(kind(1.0d0))), dimension(ldb,*) :: res
complex(kind=(kind(1.0d0))),  dimension(:,:),pointer :: temp_upper
complex(kind=(kind(1.0d0))),  dimension(:,:),pointer :: temp_lower
integer :: i,j,p
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
integer :: infoloc

if(trans .eq. 'n' .or. trans .eq. 'N') then 

  call wallocate_2z(temp_upper,ku,spm(20)*nrhs,infoloc)
  call wallocate_2z(temp_lower,kl,spm(20)*nrhs,infoloc)
  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
    if(p<spm(20)) then
      call ZLACPY('F',ku,nrhs,B(Ajmin(p+1), 1),ldb,temp_upper(1,(p-1)*nrhs+1),ku)
    endif
    if(p>1) then
      call ZLACPY('F',kl,nrhs,B(Ajmin(p)-kl,1),ldb,temp_lower(1,(p-1)*nrhs+1),kl)
    endif
  enddo

  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
    !Use the fancy structurally symmetric multiply if we can, otherwise use the general multiply. 
    if(kl == ku) then
      call ZSBMM('F', Ajmin(p+1)-Ajmin(p), nrhs, KL, -z_one, C(1,Ajmin(p)), ldc, B(Ajmin(p),1),LDB,z_one, res(Ajmin(p),1), LDB)
    else
      call ZGBMM(trans,'N', Ajmin(p+1)-Ajmin(p), nrhs, KL, KU, -z_one, C(1,Ajmin(p)), ldc, B(Ajmin(p),1),LDB,z_one, res(Ajmin(p),1), LDB)
    endif
!The triangular bits
    if(p<spm(20)) then
      call ZTRMM('L','L',trans,'N', ku, nrhs, z_one, C(1,Ajmin(p+1)), ldc-1, temp_upper(1,(p-1)*nrhs+1),ku)
      do i=1,nrhs 
        call ZAXPY(ku, -z_one, temp_upper(1,(p-1)*nrhs+i), 1, res(Ajmin(p+1)-ku,i), 1)
      enddo
    endif
    if(p>1) then
      call ZTRMM('L','U',trans,'N', kl, nrhs, z_one, C(ku+kl+1,Ajmin(p)-kl), ldc-1, temp_lower(1,(p-1)*nrhs+1),kl)
      do i=1,nrhs 
        call ZAXPY(kl, -z_one, temp_lower(1,(p-1)*nrhs+i), 1, res(Ajmin(p),i), 1)
      enddo
    endif
  enddo

else
  call wallocate_2z(temp_upper,kl,spm(20)*nrhs,infoloc)
  call wallocate_2z(temp_lower,ku,spm(20)*nrhs,infoloc)

  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
    if(p<spm(20)) then
      call ZLACPY('F',kl,nrhs,B(Ajmin(p+1), 1),ldb,temp_upper(1,(p-1)*nrhs+1),kl)
    endif
    if(p>1) then
      call ZLACPY('F',ku,nrhs,B(Ajmin(p)-ku,1),ldb,temp_lower(1,(p-1)*nrhs+1),ku)
    endif
  enddo

  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
    call ZGBMM(trans,'N', Ajmin(p+1)-Ajmin(p), nrhs, KU, KL, -z_one, C(1,Ajmin(p)), ldc, B(Ajmin(p),1),LDB,z_one, res(Ajmin(p),1), LDB)
  enddo

  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
!The triangular bits
    if(p>1) then
    !Multiply by the upper band-- so this should be the stuff that was pulled from the current partition (because transpose)
      call ZTRMM('L','L',trans,'N', ku, nrhs, z_one, C(1,Ajmin(p)), ldc-1, temp_lower(1,(p-1)*nrhs+1),ku)
      do i=1,nrhs 
        call ZAXPY(ku, -z_one, temp_lower(1,(p-1)*nrhs+i), 1, res(Ajmin(p),i), 1)
      enddo
    endif

    if(p<spm(20)) then
      call ZTRMM('L','U',trans,'N', kl, nrhs, z_one, C(ku+kl+1,Ajmin(p+1)-kl), ldc-1, temp_upper(1,(p-1)*nrhs+1),kl)
      do i=1,nrhs 
        call ZAXPY(kl, -z_one, temp_upper(1,(p-1)*nrhs+i), 1, res(Ajmin(p+1)-kl,i), 1)
      enddo
    endif
  enddo
endif

call wdeallocate_2z(temp_upper)
call wdeallocate_2z(temp_lower)

end subroutine zspike_matmul





subroutine zspike_multi(invred,klu,nrhs,A,xb,ldb,xt,ldt)
! Purpose:
! Solve of the reduced system A (different format possible) by x (xb,xt), result in x
! This is an internally used function, unlikely to be useful for the user.
! --------
! Arguments: 
!       A (input) complex(kind=(kind(1.0d0))) array, input matrix 
!       xb (input/output) complex(kind=(kind(1.0d0))) array, the multiplier and bottom half
!       of the result 
!       xt (input/output) complex(kind=(kind(1.0d0))) array, the multiplier and the top
!       half of the result  
!       ldb,ldt leading dimensions
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer:: klu,nrhs,ldb,ldt
complex(kind=(kind(1.0d0))), dimension(2*klu,*):: A
complex(kind=(kind(1.0d0))), dimension(ldb,*) :: xb
complex(kind=(kind(1.0d0))), dimension(ldt,*) :: xt
logical :: invred
integer:: Info_lap
complex(kind=(kind(1.0d0))), dimension(:,:),pointer :: A_aux
integer, dimension(:),pointer :: IPIV
character(len=1) :: TRANSA,TRANSB,SIDE
integer :: t1,t2,tim,t3,t4,tim2
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
integer :: infoloc


TRANSA='N'
TRANSB = 'N'
SIDE='L'

  call wallocate_1i(IPIV,klu,infoloc)
  call wallocate_2z(A_aux,klu,klu,infoloc)


!!! step (1) solve (I-WV)xt=xt-W*xb
  call ZGEMM (TRANSA, TRANSB, klu, nrhs, klu, -z_one, A(klu+1,1), 2*klu, xb, ldb,z_one, xt, ldt) 
!!! step (1) solve (I-WV)xt=xt-W*xb
  call ZLACPY('A',klu,klu,A,2*klu,A_aux,klu)

  call ZGESV( klu, nrhs, A_aux, klu, IPIV, xt, ldt, INFO_lap) 
!!! step (2) compute xb=xb-Vxt
  call ZGEMM (TRANSA, TRANSB, klu, nrhs, klu, -z_one, A(1,klu+1), 2*klu, xt, ldt,z_one, xb, ldb) 


  call wdeallocate_1i(IPIV)
  call wdeallocate_2z(A_aux)
end subroutine  zspike_multi





 subroutine  zspike_invred(invred,klu,V,W,redA)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Purpose: Invert  or store the truncated reduced system by blocs
! matrix to invert as the following form
!
!                |I V|^(-1) or |I-WV       V|
!                |W I|         |W          I|
! Arguments: 
!       V (input) complex(kind=(kind(1.0d0))) array, 1x2 block of input matrix 
!       W (input) complex(kind=(kind(1.0d0))) array, 2x1 block of the input matrix 
!       readA (output) complex(kind=(kind(1.0d0))) array, result matrix  
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer :: klu,i,info
complex(kind=(kind(1.0d0))), dimension(2*klu,*):: redA
complex(kind=(kind(1.0d0))), dimension(2*klu,*):: V
complex(kind=(kind(1.0d0))), dimension(2*klu,*):: W
logical:: invred
complex(kind=(kind(1.0d0))), dimension(:),pointer :: work
integer, dimension(:),pointer :: ipiv
character(len=1) :: TRANSA,TRANSB
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
integer :: infoloc



redA(1:klu,1:klu)=z_zero
        

redA(klu+1:2*klu,1:klu)=W(1:klu,1:klu)
redA(1:klu,klu+1:2*klu)=V(klu+1:2*klu,1:klu)
redA(klu+1:2*klu,klu+1:2*klu) = z_zero
do i=1,klu
  redA(i,i)=z_one
  redA(klu+i,klu+i)=z_one
end do


!! Step I-WV
TRANSA='N'
TRANSB = 'N'
call ZGEMM ( TRANSA, TRANSB, klu, klu, klu, -z_one, redA(klu+1,1), 2*klu, redA(1,klu+1), 2*klu,z_one, redA(1,1), 2*klu ) 

end subroutine zspike_invred




subroutine zspike_prep_recn(invred,kl,ku,rV,rW,redA,nbpart,nblevel) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepares the reduced system. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer :: nbpart,nblevel 
logical :: invred
integer,external :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads
integer :: j,sub,nbsub,l,i,kl,ku,rank,klu,s,k
integer:: p, myid
logical::test
integer::lck,a
complex(kind=(kind(1.0d0))),dimension(2*max(kl,ku),max(kl,ku),nblevel,*) :: rV,rW
complex(kind=(kind(1.0d0))),dimension(2*max(kl,ku),2*max(kl,ku),*) :: redA
complex(kind=(kind(1.0d0))), dimension(:,:,:), pointer :: aux1,aux2
character(len=1) :: TRANSA,TRANSB
integer :: t1,t2,tim,t3,t4,tim2
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
integer :: infoloc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
klu=max(kl,ku)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
call wallocate_3z(aux1,klu,klu,nbpart,infoloc)
call wallocate_3z(aux2,klu,klu,nbpart,infoloc)
!!!!!
!!! We associate to each CPU 1.0d0 particular reduce system (which is required for the recursion)
!!! Example(nbpart=8) for level 1 we get the reduce systems for CPU=1,3,5,7  
    !!                    for level 2 we get the reduce systems for CPU=2,6    
    !!                    for level 3 we get the reduce systems for CPU=4 
!!!!!

TRANSA='N'
TRANSB ='N'


do j=1,nblevel-1
!!! compute the spikes level j
  nbsub=2**(nblevel-j) !! Number of subsystem to solve for the level j
  sub=2**(j-1) !! Number of spikes by subsystem

!$omp parallel do private(l,p,t1,t2,t3,t4,tim)
  do i=1,nbpart,1
    p=i
    if ((mod(p-sub,nbpart/nbsub) == 0) .and. p < nbpart-nbpart/nbsub) then
      aux1(:,:,p)=z_zero
      call ZLACPY('A',klu,klu,rV(1,1,j,p+1),2*klu,aux2(1,1,p),klu)
      call zspike_multi(invred,klu,klu,redA(1,1,p),aux1(1,1,p),klu,aux2(1,1,p),klu) !solution of the  truncated reduce system 

      if (j>1) then
        do l=p-sub+1,p-1  !! CPU(/=i) involved in the calculation of the new V spikes (Vtop)
            !send i -> l  aux2 
            !recv l <- i  aux2
          call ZGEMM ( TRANSA, TRANSB, 2*klu, klu, klu, -z_one, rV(1,1,j,l), 2*klu, aux2(1,1,p), klu,z_zero, rV(1,1,j+1,l), 2*klu )
        enddo
      endif
      call ZGEMM ( TRANSA, TRANSB, klu, klu, klu, -z_one, rV(1,1,j,p), 2*klu, aux2(1,1,p), klu,z_zero, rV(1,1,j+1,p), 2*klu )
      call ZLACPY('A',klu,klu,aux1(1,1,p),klu,rV(klu+1,1,j+1,p),2*klu)

      do l=p+1,p+sub !! CPU(/=i) involved in the calculation of the new V spikes (Vbottom)
        !send i -> l  aux1
        !recv l <- i  aux1
        if (l/=p+1) then
          call ZLACPY('A',2*klu,klu,rV(1,1,j,l),2*klu,rV(1,1,j+1,l),2*klu)
          call ZGEMM ( TRANSA, TRANSB, 2*klu, klu, klu, -z_one, rW(1,1,j,l), 2*klu, aux1(1,1,p), klu,z_one,rV(1,1,j+1,l), 2*klu )
        else
          !send i -> i+1  aux2 
          !recv i+1 <- i  aux2
          call ZLACPY('A',klu,klu,aux2(1,1,p),klu,rV(1,1,j+1,l),2*klu)
          call ZLACPY('A',klu,klu,rV(klu+1,1,j,l),2*klu,rV(klu+1,1,j+1,l),2*klu)
          call ZGEMM ( TRANSA, TRANSB, klu, klu, klu, -z_one, rW(klu+1,1,j,l), 2*klu, aux1(1,1,p), klu,z_one, rV(klu+1,1,j+1,l), 2*klu )
        endif
      enddo
    endif

    p=i-1
    if ((p>= sub+nbpart/nbsub) .and. (mod(p-sub,nbpart/nbsub)==0)) then
      aux2(:,:,p-1)=z_zero
      call ZLACPY('A',klu,klu,rW(klu+1,1,j,p),2*klu,aux1(1,1,p-1),klu)
      call zspike_multi(invred,klu,klu,redA(1,1,p),aux1(1,1,p-1),klu,aux2(1,1,p-1),klu) !! solution of the  truncated reduce system 

      if(j>1) then
        do l=p-sub+1,p-1  !! CPU(/=i) involved in the calculation of the new W spikes (Wtop)
          !send i -> l  aux4 
          !recv l <- i  aux4
          call ZLACPY('A',2*klu,klu,rW(1,1,j,l),2*klu,rW(1,1,j+1,l),2*klu)
          call ZGEMM ( TRANSA, TRANSB, 2*klu, klu, klu, -z_one, rV(1,1,j,l), 2*klu, aux2(1,1,p-1), klu,z_one, rW(1,1,j+1,l), 2*klu )
        enddo
      endif

      call ZLACPY('A',klu,klu,rW(1,1,j,p),2*klu,rW(1,1,j+1,p),2*klu)
      call ZGEMM ( TRANSA, TRANSB, klu, klu, klu, -z_one, rV(1,1,j,p), 2*klu, aux2(1,1,p-1), klu,z_one, rW(1,1,j+1,p), 2*klu )
      call ZLACPY('A',klu,klu,aux1(1,1,p-1),klu,rW(klu+1,1,j+1,p),2*klu)

      do l=p+1,p+sub !! CPU(/=i) involved in the calculation of the new W spikes (Wbottom)
        !send i -> l  aux3
        !recv l <- i  aux3
        if (l/=p+1) then
          call ZGEMM ( TRANSA, TRANSB, 2*klu, klu, klu, -z_one, rW(1,1,j,l), 2*klu, aux1(1,1,p-1), klu,z_zero, rW(1,1,j+1,l), 2*klu ) 
        else
          !send i -> i+1  aux4 
          !recv i+1 <- i  aux4
          call ZLACPY('A',klu,klu,aux2(1,1,p-1),klu,rW(1,1,j+1,l),2*klu)
          call ZGEMM ( TRANSA, TRANSB, klu, klu, klu, -z_one, rW(klu+1,1,j,l), 2*klu, aux1(1,1,p-1), klu,z_zero, rW(klu+1,1,j+1,l), 2*klu )
        endif
      enddo
    endif
  enddo

!!!<<< all threads synchronize before computing reduce system for level j+1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! New reduce systems for level j+1 !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nbsub=2**(nblevel-j-1) !! Number of subsystem to solve for the level j
  sub=2**j !! Number of spikes by subsystem

!$omp parallel do default(shared) private(p)
  do i=1,nbpart
    p=i
    if ((mod(p-sub,nbpart/nbsub) == 0)) then
!!!!! Inverse of the truncated reduced system  
      call zspike_invred(invred,klu,rV(1,1,j+1,p),rW(1,1,j+1,p+1),redA(1,1,p)) !! give the inverse of the reduce system in the new bloc numbering
    endif
  enddo
enddo

! Remove Fortran runtime dependency
call wdeallocate_3z(aux1)
call wdeallocate_3z(aux2)
! End of removal

end subroutine zspike_prep_recn





subroutine zspike_solve_recn(invred,kl,ku,nrhs,rV,rW,redA,nbpart,nblevel,grj)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Solves the reduced system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer :: nbpart,nblevel
logical :: invred 
integer,external :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads
integer :: j,sub,nbsub,l,i,kl,ku,rank,klu,s,p,nrhs
character(len=1) :: TRANSA,TRANSB
complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: aux1,aux2
complex(kind=(kind(1.0d0))),dimension(2*max(kl,ku),max(kl,ku),nblevel,*) :: rV,rW
complex(kind=(kind(1.0d0))),dimension(2*max(kl,ku),2*max(kl,ku),*) :: redA
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),nrhs,*) :: grj
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
integer :: infoloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
klu=max(kl,ku)
!s=size(grj,2)
s=nrhs

! Remove Fortran runtime dependency
call wallocate_3z(aux1,klu,s,nbpart,infoloc)
call wallocate_3z(aux2,klu,s,nbpart,infoloc)
        
! End of removal

TRANSA='N'
TRANSB='N'

do j=1,nblevel

!!! compute the spikes level j
  nbsub=2**(nblevel-j) !! Number of subsystem to solve for the level j
  sub=2**(j-1) !! Number of spikes by subsystem
  aux1=z_zero
  aux2=z_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! CALCULATION INVOLVING  V and W SPIKES at level j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$omp parallel do
  do i=sub,nbpart,nbpart/nbsub
    call ZLACPY('A',klu,nrhs,grj(1,1,i+1),  2*klu,aux2(1,1,i),klu)
    call zspike_multi(invred,klu,nrhs,redA(1,1,i),grj(klu+1,1,i),2*klu,aux2(1,1,i),klu) !! solution of the  truncated reduce system 
    call ZLACPY('A',klu,nrhs,grj(klu+1,1,i),2*klu,aux1(1,1,i),klu)
  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$omp parallel do
  do i=sub,nbpart,nbpart/nbsub
    do l=i-sub+1,i-1  !! CPU(/=i) involved in the calculation with the V spikes 
             !send i -> l  aux2 
             !recv l <- i  aux2
      call ZGEMM ( TRANSA, TRANSB, 2*klu, s, klu, -z_one, rV(1,1,j,l), 2*klu, aux2(1,1,i), klu,z_one, grj(1,1,l), 2*klu )
    end do
    call ZGEMM ( TRANSA, TRANSB, klu, s, klu, -z_one, rV(1,1,j,i), 2*klu, aux2(1,1,i), klu,z_one, grj(1,1,i), 2*klu )
    do l=i+1,i+sub  !! CPU(/=i) involved in the calculation of the W spikes 
             !send i -> l  aux1
             !recv l <- i  aux1
      if (l/=i+1) then
                   !         grj(:,:)=grj(:,:)-matmul(sW(:,:,j),aux1) 
        call ZGEMM ( TRANSA, TRANSB, 2*klu, s, klu, -z_one, rW(1,1,j,l), 2*klu, aux1(1,1,i), klu,z_one, grj(1,1,l), 2*klu )
      else
                   !send i -> i+1  aux2 
                   !recv i+1 <- i  aux2
        call ZLACPY('A',klu,nrhs,aux2(1,1,i),klu,grj(1,1,l),2*klu)
        call ZGEMM ( TRANSA, TRANSB, klu, s, klu, -z_one, rW(klu+1,1,j,l), 2*klu, aux1(1,1,i), klu,z_one, grj(klu+1,1,l), 2*klu )
      end if
    end do
  end do
enddo

! Remove Fortran runtime dependency 
call wdeallocate_3z(aux1)
call wdeallocate_3z(aux2)
! End of removal
end subroutine zspike_solve_recn



subroutine zspike_solve_recn_transpose(trans,kl,ku,nrhs,rV,rW,redA,nbpart,nblevel,grj) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Transpose solve for the reduced system 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================================================================
!  Braegan Spring - 2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
character :: trans
integer   :: klu,kl,ku,nrhs
integer   :: nblevel, nbpart, nbsub, nspiketips, p
integer   :: i,j,k,l
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),max(kl,ku),nblevel,*) :: rW,rV
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),2*max(kl,ku),*)       :: redA
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),nrhs,*)               :: grj
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)

klu=max(kl,ku)
! grj holds the right hand side
! It is broken up by the number of partitons, because each partition works on a section at a time.
! redA holds the reduced systems 
! These are all the same size, and there is one per partition. 

! rV and rW hold the spike tips.
! The storage here is a bit weird
! The spikes for each level are part of the D matrix for the next level.
! They double in size as we move up the levels, starting at a size of 2xklu
! But these tips are stored in a weird way, because the dimensions of rV and rW are 2*klu x klu...
! Because this is the transpose case, we're going 'backwards' across the levels of the recursive systems. 

do i=nblevel,1,-1
! nsub is the number of subsystems that we'll be working on for this level. 
! it starts at 1, because the first matrix we work on is s_3 (because everything has been reversed in the transpose case), which has only one reduced system (in the middle)

! p = (2*j-1)*2**(i-1) maps j and i to the partiton (subsystem), p 
! The idea:
! Take the number of partitions total, and divide that by the number of partitions we'll be working on for this level to get the number of
! partitions we step over to get from one "partition we are working on in this level" to the next. 
! The starting partition is gotten by dividing the step size by two. 
! But since everything is a power of two, it can be done as below. 
! ... it makes sense if you draw it on graph paper!
  nbsub = 2**(nblevel-i)

! Build the modified right hand side 
  !$OMP PARALLEL DO default(none) private(p,k) shared(trans,i,klu,nrhs,rV,grj,nbpart,nbsub,rW)
  do j=1,nbsub

    p = (2*j-1)*2**(i-1)
    ! Now we have to build the modified right hand side.
    call ZGEMM(trans,'N',klu,nrhs,  klu,-z_one,rV(1,1,i,p),2*klu,grj(1,1,p),2*klu,z_one,grj(1,1,p+1),2*klu)
    do k=1,nbpart/(2*nbsub)-1,1
      call ZGEMM(trans,'N',klu,nrhs,2*klu,-z_one,rV(1,1,i,p-k),2*klu,grj(1,1,p-k),2*klu,z_one,grj(1,1,p+1),2*klu)
    enddo
  enddo

  !$OMP PARALLEL DO default(none) private(p,k) shared(trans,i,klu,nrhs,rV,grj,nbpart,nbsub,rW)
   do j=1,nbsub
    p = (2*j-1)*2**(i-1)
    call ZGEMM(trans,'N',klu,nrhs,  klu,-z_one,rW(klu+1,1,i,p+1),2*klu,grj(klu+1,1,p+1),2*klu,z_one,grj(klu+1,1,p),2*klu)
    do k=1,nbpart/(2*nbsub)-1,1
      call ZGEMM(trans,'N',klu,nrhs,2*klu,-z_one,rW(1    ,1,i,p+k+1),2*klu,grj(1,1  ,p+k+1),2*klu,z_one,grj(klu+1,1,p),2*klu)
    enddo
  enddo
  !$OMP PARALLEL DO default(shared) private(p)
  do j=1,nbsub
    p = (2*j-1)*2**(i-1)
    call zspike_multi_transpose(trans,klu,nrhs,redA(1,1,p),grj(1,1,p+1),2*klu,grj(klu+1,1,p),2*klu)
  enddo
enddo

end subroutine zspike_solve_recn_transpose




subroutine zspike_multi_transpose(trans,klu,nrhs,red,xb,ldxb,xt,ldxt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recursive spike calls SPIKE on the reduced system for the main system.
! This solves the reduced system for the SPIKE calls on the main reduced system.
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================================================================
!  Braegan Spring - 2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer   :: klu, nrhs,i,ldxb,ldxt
integer   :: info
character :: trans
integer, dimension(klu)                   :: IPIV
complex(kind=(kind(1.0d0))), dimension(2*klu,*)    :: xb,xt
complex(kind=(kind(1.0d0))), dimension(2*klu,*) :: red
!complex(kind=(kind(1.0d0))), dimension(klu,*)     :: red_aux 
complex(kind=(kind(1.0d0))), pointer, dimension(:,:)     :: red_aux 
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
integer :: infoloc

call wallocate_2z(red_aux,klu,klu,infoloc)
! The mathematical format of the reduced system after block LU is 
!(where w and v are actually the tips of their respective spikes):
! __   __    __       __
! |I  0 |    | I    V  |
! |W  I |    | 0  I-WV |
! --   --    --       --
! Which transposes to 
! __         __    __     __
! |I  0       |    | I  WT | = LU
! |VT (I-WV)T |    | 0  I  |
! --         --    --     --
! But it is actually stored with (I-WV) in the top-left corner
! Because of the identities in the top left and bottom right corners, 
! most of this solve can actually be done with DGEMM's

call ZLACPY('A',klu,klu,red,2*klu,red_aux,klu)
! Set up modified right hand side (to take advantage of the fact that we already know the top
! part of the vector(s) that we are solving the reduced system for)
call ZGEMM(trans,'N',klu,nrhs,klu,-z_one,red(1,klu+1),2*klu,xt,ldxt,z_one,xb,ldxb)

call ZGETRF(klu,klu,red_aux,klu,ipiv,info)
call ZGETRS(trans,klu,nrhs,red_aux,klu,ipiv,xb,ldxb,info)
! Now for U we already know the bottom part of the solution, so we just need to do a multiply and 
! subtract for the top
call ZGEMM(trans,'N',klu,nrhs,klu,-z_one,red(klu+1,1),2*klu,xb,ldxb,z_one,xt,ldxt)
call wdeallocate_2z(red_aux)
end subroutine zspike_multi_transpose




!The simple two partition case.
subroutine zspike_GBTRFk2(n,kl,ku,A,lda,info,keys)
!=====================================================================
!  Braegan Spring  - 2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer :: i,kl,ku,klu,lda,n,info,n1,n2,k,j
complex(kind=(kind(1.0d0))), dimension(lda,*) :: A
!complex(kind=(kind(1.0d0))), dimension(max(kl,ku),2*max(kl,ku)) :: vw
integer,external :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads,omp_in_parallel,omp_get_ancestor_thread_num
!complex(kind=(kind(1.0d0))), dimension(kl,kl) :: w
!complex(kind=(kind(1.0d0))), dimension(ku,ku) :: v
complex(kind=(kind(1.0d0))), pointer, dimension(:,:) :: w
complex(kind=(kind(1.0d0))), pointer, dimension(:,:) :: v
real(kind=(kind(1.0d0))) :: norma,nzero
double precision :: t1,t2,t3,t4,t5,t6
integer, dimension(4) :: keys
complex(kind=(kind(1.0d0))),parameter :: z_zero=(0.0d0,0.0d0)
real(kind=(kind(1.0d0))), parameter :: d_zero_res=0.0d0
integer :: infoloc

call wallocate_2z(w,kl,kl,infoloc)
call wallocate_2z(v,ku,ku,infoloc)

klu=max(kl,ku)
norma = d_zero_res
nzero = d_zero_res

if(mod(n,2)==0) then
  n1 = n/2
  n2 = n/2
else
  n1 = n/2
  n2 = n/2+1
endif
! Conceptually the thing we're pulling out of the top is a lower triangular
! and the thing we're pulling out of the bottom is upper triangular.
! But they look opposite because banded is weird.

! i=omp_get_thread_num()+1
! Top partition
! Solving A[n/2-ku:n/2,n/2-(kl+ku):n/2]*V[n/2-ku:n/2]=B=A[n/2-ku:n/2,n/2:n/2+ku]
! But in banded.

!lu factorize

! Here is the logic that will make us happy even if this is only called with one thread:
! If there is no block, block other threads, and claim the key for the first section. 
! Then, unblock other threads and go to factorize the first partition. 
! If the key has been claimed for the first partition by some other thread, jump over that code. 
! Place a block after the end of the first partition factorization section. 
! Claim the key for the second partition factorization if it is unclaimed. 
! And so on. 

! Using this scheme we will be able to avoid talking to all other threads/using critical sections if it is unnecessary...

! Test the key. If the key is not claimed, enter a critical region. Test again (because the key may have been claimed between the time you checked and the time you entered the critical region)

if(keys(2) == -1) then
  !$omp critical (keyclaim_region_1) 
  if(keys(2) == -1) then
    keys(2) = omp_get_thread_num()
    !$omp flush
  endif
  !$omp end critical (keyclaim_region_1)
endif

if(keys(2) == omp_get_thread_num()) then
!  call wallocate_2z(v,ku,ku,infoloc)
  v=z_zero
  call ZGBALU(n1,kl,ku,A(1,1),lda,nzero,norma,info)
  call ZLACPY('L',ku,ku,A(1,n1+1),kl+ku,v(1,1),ku) 
  call ZTBSM('L','N','U', ku, ku, kl, A(ku+1,n1+1-ku), lda, v(1,1), ku)
  call ZLACPY('L',ku,ku,v(1,1),ku,A(1,n1+1),lda-1)
endif

if(keys(3) == -1) then
!$omp critical (keyclaim_region_2) 
  !$omp flush
  if(keys(3) == -1) then
    keys(3) = omp_get_thread_num()
    endif
!$omp end critical (keyclaim_region_2)
endif

if(keys(3) == omp_get_thread_num()) then
!  call wallocate_2z(w,kl,kl,infoloc)
  w=z_zero
  call ZGBAUL(n2,kl,ku,A(1,n1+1),lda,nzero,norma,info)
  call ZLACPY('U',kl,kl,A(kl+ku+1, n1+1-kl),kl+ku,w(1,1),kl)
  call ZTBSM('U','N','U',kl, kl, ku, A(1, n1+1), lda, w(1,1), kl)
  call ZLACPY('U',kl,kl,w(1,1),kl,A(kl+ku+1, n1+1-kl),kl+ku)
endif
call wdeallocate_2z(w)
call wdeallocate_2z(v)
end subroutine zspike_GBTRFk2








subroutine zspike_GBTRSk2(trans,n,nrhs,kl,ku,A,lda,B,ldb,keys,red,grj,zeroskip)
!=====================================================================
!  Braegan Spring  2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer                              :: n,nrhs,kl,ku,lda,ldb,info
complex(kind=(kind(1.0d0))), dimension(lda,*)   :: A
complex(kind=(kind(1.0d0))), dimension(ldb,*):: B
logical :: zeroskip
character                            :: trans
integer, dimension(4)                :: keys
!complex(kind=(kind(1.0d0))), target, dimension(2*max(kl,ku),2*max(kl,ku)+nrhs) :: work
!complex(kind=(kind(1.0d0))), dimension(:,:), pointer :: red
!complex(kind=(kind(1.0d0))), dimension(:,:), pointer :: grj
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),*) :: red
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),*) :: grj 
!integer,external :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads,omp_in_parallel

if(trans == 'N' .or. trans=='n') then
!  red => work(1:2*max(kl,ku),1:2*max(kl,ku))
!  grj => work(1:2*max(kl,ku),2*max(kl,ku)+1:2*max(kl,ku)+nrhs)
  call zspike_simple_solve(n,nrhs,kl,ku,A,lda,B,ldb,info,keys,red,grj,zeroskip)
else
!if(omp_get_thread_num() == 1) then
!endif
!  red => work(1:2*max(kl,ku),1:2*max(kl,ku))
!  grj => work(1:2*max(kl,ku),2*max(kl,ku)+1:2*max(kl,ku)+nrhs)
  call zspike_simple_solve_transpose(n,trans,nrhs,kl,ku,A,lda,B,ldb,info,keys,red,grj)
endif

end subroutine zspike_GBTRSk2




subroutine zspike_simple_solve(n,nrhs,kl,ku,A,lda,f,ldf,info,keys,red,grj,zeroskip)
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer :: j,i,kl,ku,klu,p,nrhs,local_start,local_nrhs,n,lda,info,t1,t2,t3,tim,n1,n2,ldf
logical :: invred
logical :: zeroskip
complex(kind=(kind(1.0d0))), dimension(lda,*):: A
integer,external :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads,omp_in_parallel
complex(kind=(kind(1.0d0))), dimension(ldf,*):: f
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),*) :: red
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),*) :: grj
!complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),2*max(kl,ku)+nrhs) :: work
!complex(kind=(kind(1.0d0))), pointer, dimension(:,:):: red
!complex(kind=(kind(1.0d0))), pointer, dimension(:,:) :: grj
!complex(kind=(kind(1.0d0))) :: nzero, norma
real(kind=(kind(1.0d0))) :: nzero, norma
!integer, dimension(2*max(kl,ku)) :: ipiv
integer, dimension(2*max(kl,ku)) :: ipiv
integer , dimension(4) :: keys
complex(kind=(kind(1.0d0))),parameter :: z_one=(1.0d0,0.0d0),z_zero=(0.0d0,0.0d0)
integer :: infoloc

invred=.false.
klu=max(kl,ku)

!Note : In general, keys(2) belongs to the one that works one the top part of the f vector, and keys(3) belongs to the one that works on the bottom part of the f vector

!call wallocate_1i(ipiv,2*max(kl,ku),infoloc)
!allocate(ipiv(2*max(kl,ku)))

!allocate(red(2*klu,2*klu))
!allocate(grj(2*klu,nrhs))


if(mod(n,2)==0) then
  n1 = n/2
  n2 = n/2
else
  n1 = n/2
  n2 = n/2+1 
endif

! Automatic keyclaiming region -- if you haven't set the keys the first thread to get here will grab a key.

if(keys(2) == -1) then
  !$omp critical (keyclaim_region_3) 
  if(keys(2) == -1) then
    keys(2) = omp_get_thread_num()
    !$omp flush
  endif
  !$omp end critical (keyclaim_region_3)
endif

if(keys(2) == omp_get_thread_num()) then   
!! Solving Ajgj=fj
! Unpacking the semi-solved v from A and doing the U sweep to form the tip of v.
  if(zeroskip) then
    local_start = ku+1
    local_nrhs  = kl
  else
    local_start = 1
    local_nrhs  = nrhs
  endif

  red(1:klu,1:2*klu)=z_zero
  do i=1,klu
    red(i,i)=z_one
  enddo
!  grj(1:klu,1:nrhs) = 0.0d0
! Now begin the real work of the solve stage
  call ZTBSM('L','N','U', n1, local_nrhs,kl,A(ku+1,1),lda,f(1,local_start),ldf)
  call ZLACPY('F', klu, nrhs,f(n1+1-klu,1),ldf, grj(1,1),2*klu)
  call ZTBSM('U','N','N',klu,nrhs,ku,A(1,n1+1-klu),lda,grj(1,1),2*klu)

  call ZLACPY('L',ku,ku,A(1,n1+1),kl+ku,red(klu-ku+1,1+klu),2*klu)
  call ZTBSM('U','N','N', klu, ku, ku, A(1,n1+1-klu), lda, red(1,1+klu), 2*klu)

  !$omp atomic
  keys(1)= keys(1)+ 1
endif


! Automatic keyclaiming region -- if you haven't set the keys the first thread to get here will grab a key.
if(keys(3) == -1) then
  !$omp critical (keyclaim_region_4) 
  if(keys(3) == -1) then
    keys(3) = omp_get_thread_num()
    !$omp flush
  endif
  !$omp end critical (keyclaim_region_4)
endif

if(keys(3) == omp_get_thread_num()) then
!! solving Anbpart gnbpart =f nbpart
! Unpack the semi-solved w from A and do the L sweep.

  if(zeroskip) then
    local_nrhs  = ku
  else
    local_nrhs  = nrhs
  endif

  red(klu+1:2*klu,1:2*klu) = z_zero
  do i=1,klu
    red(klu+i,klu+i)=z_one
  enddo
!  grj(klu+1:2*klu,1:nrhs) = 0.0d0

  call ZTBSM('U','N','U', n2, local_nrhs, ku,A(1,n1+1),lda,f(n1+1,1),ldf)
  call ZLACPY('F', klu, nrhs, f(n1+1,1),ldf, grj(klu+1,1), 2*klu )
  call ZTBSM('L','N','N',klu,nrhs,kl,A(ku+1,n1+1),lda,grj(klu+1,1),2*klu)

  call ZLACPY('U',kl,kl,A(kl+ku+1,n1+1-kl),kl+ku,red(klu+1,klu-kl+1),2*klu)
  call ZTBSM('L','N','N',klu, kl, kl, A(ku+1,n1+1), lda, red(klu+1,klu-kl+1),2*klu)

  !$omp atomic
  keys(1)= keys(1)+ 1
endif 

!do j=1,2*klu
!enddo

do while(keys(1)< 2)
!$omp flush
enddo

if(keys(2) == omp_get_thread_num()) then
  !$omp flush
  call ZGESV( 2*klu, nrhs, red, 2*klu, ipiv, grj, 2*klu, info) 
  !$omp atomic
  keys(1)= keys(1)+ 1
endif


do while(keys(1)< 3)
!$omp flush
enddo

!<<<<<<<<<<<<<<<<<<<<<<<< Retrieval  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

if(keys(2) == omp_get_thread_num()) then
  call ZTRMM('L','L','N','N',ku,nrhs,-z_one,A(1,n1+1),kl+ku,grj(klu+1,1),2*klu)
  do i=1,nrhs
    call ZAXPY(ku,z_one,grj(klu+1,i),1,f(n/2+1-ku,i),1) 
  enddo
  call ZTBSM('U','N','N', n1, nrhs, ku,A(1,1),lda,f(1,1),ldf)
  !$omp atomic
  keys(1)= keys(1)+ 1
endif

if(keys(3) == omp_get_thread_num()) then
  call ZTRMM('L','U','N','N',kl,nrhs,-z_one,A(kl+ku+1,n1+1-kl),kl+ku,grj(klu-kl+1,1),2*klu)
  do i=1,nrhs
    call ZAXPY(kl,z_one,grj(klu-kl+1,i),1,f(n1+1,i),1) 
  enddo
  call ZTBSM('L','N','N', n2, nrhs, kl,A(ku+1,n1+1),lda,f(n1+1,1),ldf)
  !$omp atomic
  keys(1)= keys(1)+ 1
endif 

do while(keys(1)< 5)
!$omp flush
enddo
!deallocate(grj,red)

end subroutine zspike_simple_solve


subroutine zspike_simple_solve_transpose(n,trans,nrhs,kl,ku,A,lda,f,ldf,info,keys,red,grj)
!=====================================================================
!  Braegan Spring  - 2015
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer :: j,i,k,kl,ku,klu,p,nrhs,n,lda,info,n1,n2,info_red,ldf
complex(kind=(kind(1.0d0))), dimension(lda,*):: A
complex(kind=(kind(1.0d0))), dimension(ldf,*):: f
!complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),nrhs) :: aux
complex(kind=(kind(1.0d0))), pointer, dimension(:,:) :: aux
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),2*max(kl,ku)) :: red
complex(kind=(kind(1.0d0))), dimension(2*max(kl,ku),nrhs) :: grj
integer, dimension(2*max(kl,ku)) :: ipiv_red
integer :: t1,t2,t3,t4,t5,t6,t7,t8,tim,t9,infoloc
character(len=30) :: format
integer, dimension(4) :: keys
integer, external :: omp_get_thread_num
complex(kind=(kind(1.0d0))),parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
character :: trans


call wallocate_2z(aux,2*max(kl,ku),nrhs,infoloc)

klu=max(kl,ku)
!allocate(red(2*klu,2*klu))
!allocate(grj(2*klu,nrhs))
!allocate(aux(2*klu,nrhs))

if(mod(n,2)==0) then
  n1 = n/2
  n2 = n/2
else
  n1 = n/2
  n2 = n/2+1 
endif

if(keys(3) == -1) then
  !$omp critical (keyclaim_region_3) 
  if(keys(3) == -1) then
    keys(3) = omp_get_thread_num()
    !$omp flush
  endif
  !$omp end critical (keyclaim_region_3)
endif

if(keys(3) == omp_get_thread_num()) then
  call ZLACPY('F', klu, nrhs, f(n1+1,1), ldf, grj(klu+1,1), 2*klu)
  red(1:2*klu,1:klu) = z_zero
!    red(1:klu,1:2*klu) = 0.0d0
  do j=1,klu
    red(j,j)=z_one
  enddo
  !$omp atomic
  keys(1)= keys(1)+ 1
endif

if(keys(2) == -1) then
  !$omp critical (keyclaim_region_3) 
  if(keys(2) == -1) then
    keys(2) = omp_get_thread_num()
    !$omp flush
  endif
  !$omp end critical (keyclaim_region_3)
endif

if(keys(2) == omp_get_thread_num()) then
  call ZLACPY('F', klu, nrhs, f(n1+1-klu,1), ldf, grj(1,1), 2*klu)
  red(1:2*klu,1+klu:2*klu) = z_zero
!    red(klu+1:2*klu,1:2*klu) = 0.0d0
  do j=1,klu
    red(klu+j,klu+j)=z_one
  enddo
  !$omp atomic
  keys(1)= keys(1)+ 1
endif 

do while(keys(1)< 2)
!$omp flush
enddo

!Set up f-vector - save tip and augment it with zeroes
if(keys(2) == omp_get_thread_num()) then

  call ZLACPY('L',ku,ku,A(1,n1+1),lda-1,red(klu-ku+1,1+klu),2*klu)
  call ZTBSM('U','N','N', klu, ku, ku, A(1,n1+1-klu), lda, red(1,1+klu), 2*klu)

  f(n1-klu+1:n1,1:nrhs) = z_zero
! Perform the up sweep on f, then multiply it by the already L solved v (the L solve was done in the factorization stage)
  call ZTBSM('U',trans,'N',n1, nrhs, ku, A(1,1), lda, f(1,1), ldf)
  call ZLACPY('F', klu, nrhs, f(n1+1-klu,1), ldf,aux(1,1), 2*klu)

! This multiples by the v~ 
!    call ZTRMM('L','L',trans,'N',ku, nrhs, z_one, A(1,n1+1), kl+ku, f(n1-ku+1,1), n) ! "upper" band 
  call ZTRMM('L','L',trans,'N',ku, nrhs, z_one, A(1,n1+1), lda-1, aux(klu-ku+1,1), 2*klu) ! "upper" band 
! Later we'll assume that vw contains the tips of the actual v and w spikes (instead of the semi-solved version from the factorizaiton stage), so do their sweep here.
  do j=1,nrhs
!      call ZAXPY(ku,-z_one,f(n1-ku+1,j),1, grj(klu+1,j),1)
    call ZAXPY(ku,-z_one,aux(klu-ku+1,j),1, grj(klu+1,j),1)
  enddo
!  call ZLACPY('F', klu, nrhs, aux(1,1),2*klu,f(n1+1-klu,1), 2*klu)
  !$omp atomic
  keys(1)= keys(1)+ 1
endif  

if(keys(3) == omp_get_thread_num()) then
  call ZLACPY('U',kl,kl,A(kl+ku+1,n1+1-kl),lda-1,red(klu+1,klu-kl+1),2*klu)
  call ZTBSM('L','N','N',klu, kl, kl, A(ku+1,n1+1), lda, red(klu+1,klu-kl+1),2*klu)

  f(n1+1:n1+klu,1:nrhs) = z_zero
  call ZTBSM('L',trans,'N',n2, nrhs, kl, A(ku+1,n1+1), lda, f(n1+1,1), ldf)
  call ZLACPY('F', klu, nrhs, f(n1+1,1), ldf, aux(klu+1,1), 2*klu)
!    call ZTRMM('L','U',trans,'N',kl, nrhs,  1.0d0, A(kl+ku+1,n1+1-kl), kl+ku, f(n1+1,1),n) ! "lower" band 
  call ZTRMM('L','U',trans,'N',kl, nrhs, z_one, A(kl+ku+1,n1+1-kl), lda-1, aux(klu+1,1),2*klu) ! "lower" band 
  do j=1,nrhs
!      call ZAXPY(kl,-1.0d0,f(n1+1,j),1, grj(klu-kl+1,j),1)
  call ZAXPY(kl,-z_one,aux(klu+1,j),1, grj(klu-kl+1,j),1)
  enddo
!    call ZLACPY('F', klu, nrhs, aux(klu+1,1), 2*klu, f(n1+1,1), n)
  !$omp atomic
  keys(1)= keys(1)+ 1
endif

do while(keys(1)< 4)
!$omp flush
enddo

!Communication -- only effects red, which is shared
if(keys(2) == omp_get_thread_num()) then
do j=1,2*klu 
enddo

  call ZGETRF(2*klu,2*klu,red,2*klu,ipiv_red,info_red)
  call ZGETRS(trans,2*klu,nrhs,red,2*klu,ipiv_red,grj,2*klu,info_red)

do j=1,2*klu 
enddo

  !$omp atomic
  keys(1)= keys(1)+ 1
endif

do while(keys(1)< 5)
!$omp flush
enddo

if(keys(2) == omp_get_thread_num()) then
  call ZTBSM('U',trans,'N',klu, nrhs, ku, A(1,n1-klu+1), lda, grj(1,1), 2*klu)
  do j=1,nrhs
    call ZAXPY(klu,z_one,grj(1,j),1, f(n1-klu+1,j),1)
  enddo
  call ZTBSM('L',trans,'U',n1, nrhs,kl,A(1+ku,1),lda,f(1,1),ldf)
endif


if(keys(3) == omp_get_thread_num()) then
  call ZTBSM('L',trans,'N',klu,  nrhs, kl, A(ku+1,n1+1),   lda, grj(klu+1,1),2*klu)
  do j=1,nrhs
    call ZAXPY(klu,z_one,grj(klu+1,j),1, f(n1+1,j),1)
  enddo
  call ZTBSM('U',trans,'U', n2, nrhs, ku, A(1,n1+1), lda, f(n1+1,1), ldf)
endif
 
!deallocate(grj,red)
call wdeallocate_2z(aux)

end subroutine zspike_simple_solve_transpose



subroutine zspike_vector_norm(spm,norm,x,n,nrhs,dnorm)
! Purpose
! -------
! Computes the norm of a given vector or set of vectors. 
! options are norm one, norm 2 and norm infinity.
! Retains the SPIKE partitioning scheme. 
! -------
!
! Arguments  
! ---------
! spm   (in) : The spike paramater array. Described in spikeinit, in the file spike_smp_utilities.f90. Contains information about the size of each partition.
!
! norm  (in) : They type of norm computed. 0 for infinorm, 1 for norm one, and 2 for infinorm.
! 
! x     (in) : The vector or set of vectors. Dimension (N by nrhs)
!
! n     (in) : The number of rows in x. This is also the leading dimensions of x. 
!
! nrhs  (in) : The number of vectors in x.
!
! dnorm (out): Array to hold the norms of the x vectors. Dimension nrhs. The norm for the i'th column of x is stores in dnorm(i)
!
!---------
!=====================================================================
!  Braegan Spring - Eric Polizzi 2015
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
complex(kind=(kind(1.0d0))), dimension(n,*)          :: x
integer, dimension(64)                               :: spm
integer                                              :: norm
integer                                              :: n,nrhs
real(kind=(kind(1.0d0))), dimension(nrhs)              :: dnorm
real(kind=(kind(1.0d0))), dimension(:,:), pointer     :: dnormTemp
integer                                              :: imax
integer                                              :: p,i,infoloc
real(kind=(kind(1.0d0))), external                    :: dzasum,dznrm2
real(kind=(kind(1.0d0))), external                    :: dasum,dnrm2
integer, external                                    :: izamax
integer, dimension(spm(20)+1)                        :: Ajmin

call spikerl_calc_size_partitions(spm,Ajmin,n)
call wallocate_2d(dnormTemp,spm(20),nrhs,infoloc)

!One norm
if(norm==1) then
! Take the norm of each partition, then take the norm of those norms to get the overall norm. 
  !$OMP PARALLEL DO default(shared) private(i)  
  do p=1,spm(20)  
    do i=1,nrhs 
      dnormTemp(p,i) = dzasum(Ajmin(p+1)-Ajmin(p) ,x(Ajmin(p),i),1)
    enddo
  enddo

  do i=1,nrhs 
    dnorm(i) = dasum(spm(20),dnormTemp(1,i),1)
  enddo
endif


!Infinorm
if(norm==0) then
  !$OMP PARALLEL DO default(shared) private(i,imax)  
  do p=1,spm(20)  
    do i=1,nrhs 
      imax = izamax(Ajmin(p+1)-Ajmin(p),x(Ajmin(p),i),1)
      dnormTemp(p,i) = x(Ajmin(p)+imax,i)
    enddo
  enddo
  do i=1,nrhs
    imax = izamax(spm(20),dnormTemp(1,i),1)
    dnorm(i) = dnormTemp(imax,i)
  enddo
endif


!2 norm (or euclidean)
if(norm==2) then
  !$OMP PARALLEL DO default(shared) private(i)  
  do p=1,spm(20)  
    do i=1,nrhs 
      dnormTemp(p,i) = dznrm2(Ajmin(p+1)-Ajmin(p) ,x(Ajmin(p),i),1)
    enddo
  enddo

  do i=1,nrhs 
    dnorm(i) = dnrm2(spm(20),dnormTemp(1,i),1)
  enddo
endif
call wdeallocate_2d(dnormTemp)

end subroutine zspike_vector_norm

