
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
!!!!!!!!!!!!!!!!! SPIKE-SMP -  common routines                     !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! List of routines:

!!!!!! Documented !!!!!!!!
! spikeinit
! spikeinit_nthread


!!!! Auxiliary !!!!!!!!!
!spikeinit_default
!spikerl_renumber_threads
!spikerl_calc_num_levels
!spikerl_calc_num_partitions
!spikerl_calc_size_partitions
!threadmap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spikeinit(spm,n,klu)
implicit none
integer, dimension(64) :: spm
integer :: i,n,klu
integer, external :: omp_get_max_threads
!=========================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================



! Purpose
! -------
! Initializes the Spike ParaMater array to sane defaults
! This function must be called before any factorization or solve calls
! -------
!
! Arguments  
! ---------
! spm (out)
! An array of 64 integers, which control various features of the spike factorization and solve.
! Defined as follows:
! 
! User paramaters
!1  = print flag 1 for timing information, 0 for no printing (except errors, naturally)
!2  = spikeinit optimize flag. 0: Use user defined ratios. 1: Set tuning parameters to large NRHS. 2: DSPIKE_GBSV will attempt to find tuning ratios (See dspike_gbsv)
!4  = 10 times ratio large. Defaults to 22 (So ratio is 2.2)
!5  = 10 times ratio small. Defaults to 35 (Ratio 3.5)
!7  = 10 times constant K -- tuning constant for system 
!10 = number of klu*klu blocks required for work array
!11 = max number of iterations for iterative solve
!12 = exponent for residual tolerance for double precision (and double complex)-- I.E, res tolerance is 10^-spm(12)
!13 = exponent for residual tolerance for single precision (and single complex)-- I.E, res tolerance is 10^-spm(13)
!14  = Norm type used. 0=infininorm 1, 1=norm1, 2=norm2 (frob or euc depending on the situation). Norm 1 is default.
!
! Internal parameters (Unlikely to be useful to users)
!20  = #partitions for level0
!21  = #levels
!22  = #threads used by SPIKE
!23  = #two-thread partitions
!24 =  #one-thread partitions
!25  = total #threads available
! ---------


call spikeinit_default(spm)


! By default assume that we want to use every thread for spike
spm(25) = OMP_GET_MAX_THREADS()
spm(22) = spm(25)
spm(20) = spm(22)


call spikerl_calc_num_partitions(spm)
call spikerl_calc_num_levels(spm)

! Reduce number of threads if appropriate
call spikerl_renumber_threads(spm,n,klu)

spm(10) = 4*(spm(21)*spm(20)+(spm(20)-1))
!spm(26)=max(1,(spm(25)-spm(22))/spm(22))
spm(24)=spm(22)-2*spm(23)


end subroutine spikeinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spikeinit_default(spm)
implicit none
integer, dimension(64) :: spm
!=========================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================



spm=0

spm(4) = 22
spm(5) = 35
spm(7) = 16
spm(11) = 3
spm(12) = 12
spm(13) = 5
spm(14) = 1

end subroutine spikeinit_default



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spikeinit_nthread(spm,n,klu,nthread)
implicit none
integer, dimension(64) :: spm
integer :: i,nthread,n,klu
integer, external :: omp_get_max_threads
! Purpose
! -------
! Initializes the Spike ParaMater array to sane defaults and manually sets number of threads
! This function must be called before any factorization or solve calls
! -------
!
! Arguments  
! ---------
! spm (out)
! An array of 64 integers, which control various features of the spike factorization and solve.
! Defined as follows:
! 
! User paramaters
!1  = print flag 1 for timing information, 2 for partition information, 0 for no printing (except errors, naturally)
!2  = spikeinit optimize flag. 0: Use user defined ratios. 1: Set tuning parameters to large NRHS. 2: DSPIKE_GBSV will attempt to find tuning ratios (See dspike_gbsv)
!4  = 10 times ratio large. Defaults to 22 (So ratio is 2.2)
!5  = 10 times ratio small. Defaults to 35 (Ratio 3.5)
!7  = 10 times constant K -- tuning constant for system 
!10 = number of klu*klu blocks required for work array
!11 = max number of iterations for iterative solve
!12 = exponent for residual tolerance for double precision (and double complex)-- I.E, res tolerance is 10^-spm(12)
!13 = exponent for residual tolerance for single precision (and single complex)-- I.E, res tolerance is 10^-spm(13)
!14  = Norm type used. 0=norm 1, 1=infinorm, 2=norm2 (frob or euc depending on the situation). Norm 1 is default.
!
! Internal paramaters (Unlikely to be useful to users)
!20  = #partitions for level0
!21  = #levels
!22  = #threads
!23  = #two-thread paritions
!24 =  #one-thread partitions
!25  = total #threads available
! ---------
!=========================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================


call spikeinit_default(spm)


! By default assume that we want to use every thread for spike
spm(25) = nthread
spm(22) = spm(25)
spm(20) = spm(22)


call spikerl_calc_num_partitions(spm)
call spikerl_calc_num_levels(spm)

! Reduce number of threads if appropriate
call spikerl_renumber_threads(spm,n,klu)

spm(10) = 4*(spm(21)*spm(20)+(spm(20)-1))
!spm(26)=max(1,(spm(25)-spm(22))/spm(22))
spm(24)=spm(22)-2*spm(23)

end subroutine spikeinit_nthread



subroutine spikerl_calc_num_levels(spm)
! Purpose
! -------
! Calculates the number of levels of recursion for recursive spike.
! This is an internal utility routine, unlikely to be useful to the user.
!=========================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================

implicit none
integer, dimension(64) :: spm
integer i
spm(21) = 0
i = spm(20)
do while(i/=1)
  i = i/2
  spm(21) = spm(21) + 1
enddo
end subroutine spikerl_calc_num_levels


subroutine spikerl_calc_num_partitions(spm)
implicit none
integer, dimension(64) :: spm
integer, dimension(spm(20)+1) :: Ajmin
integer :: x,y,N
! Purpose
! -------
! This calculates the number of partitions, and which type they are (single- threaded large partitions, two threaded spike partitions, or single threaded small partitions)
! It is not really intended for external use.
!
!x hold the number of 1 thread partitions , except the first and last ones, because they are special
!y holds the number of 2 thread (spike) partitions 
!N holds the number of partitions total. Need this to be 2^i (a restriction of recursive spike)
!
!The only weird case is where we have 2^n - 1 partitions initially
!This makes us end up with -1 = x. So we'll just cheat and tell it we
!have 1 less partition... 
!
!using iand to check if we're 2^n - 1 
!One thread, though, is a weird special case. We just immediately leave because we're going to just call the banded prim anyway.
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================


if(spm(22) == 1) return

if(iand(spm(22)+1,spm(22)) == 0) spm(22)=spm(22)-1
spm(20) = spm(22)
y = 0 
x = spm(20) - 2
N = spm(20)

! The code below removes two single-thread partitions and adds one two thread partition, until the total number of threads is a power of two.
! (It is likely that there is a more clever way to do this, but this method is bulletproof and simple)
do while(iand(N,N-1) .NE. 0)
  y=y+1
  x=x-2
  N=N-1
enddo

spm(20) = N
spm(23) = y

end subroutine spikerl_calc_num_partitions





subroutine spikerl_calc_size_partitions(spm,Ajmin,n)
implicit none
integer, dimension(64) :: spm
integer, dimension(*) :: Ajmin
integer :: n,x,y,num_partitions,i,size_large,size_small,size_huge
double precision :: R_large, R_small
integer(8), parameter :: fout=6
!=========================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================


!Number of bigger mid-partitions (2 threaded)
y = spm(23)
!Number of smaller mid-partitions (1 threaded)
x = spm(20)-(2+spm(23))

if(spm(2) .eq. 1) then
  spm(4) = 10
!  spm(5) = 20
endif

R_large = 0.1d0 * spm(4)
R_small = 0.2d0 * spm(4)
!R_small = 0.1d0 * spm(5)

Ajmin(1) = 1

size_large = n/(R_large*(2+x/R_small+y/R_large))
size_large= size_large + mod(size_large,2)
size_small = n/(R_small*(2+x/R_small+y/R_large))

! Dumping extra elements into the first and last matrix evenly
size_huge = (n-(x*size_small + y*size_large))/2

Ajmin(2) = size_huge+1

do i=3, y+2
  Ajmin(i)=Ajmin(i-1)+size_large 
enddo

do i=y+3, spm(20)
  Ajmin(i)=Ajmin(i-1)+size_small 
enddo

Ajmin(spm(20)+1)=n+1

end subroutine spikerl_calc_size_partitions






subroutine threadmap(i,p,nbpart,n2thread)
implicit none 
integer :: i,p,nbpart,n2thread
! Purpose
! -------
! Gives the partition on which a thread should be working
! -------
!
! Arguments
! ---------
! i (in) 
! This is an index number used to identify the thread.
! For example, in an OMP_PARALLEL_DO section, each thread
! is identified by a loop index value.
!
! p (out)
! This gives the partition on which the thread should be working.
! Generally, the code in the factorization and solve functions
! filters threads by matching with this value
!
! nbpart (in) 
! The total number of partitions into which SPIKE has broken the matrix
!
! n2thread (in)
! The number of partitions which are two-threaded (spike) partitions.
! After a call to spike_calc_num_partitions (or spikeinit), 
! spm(23) contains this value.

! This method of allocating threads attempts to keep the 
! large partitions on the 'first' core, which is likely 
! to have been the one to initialize the matrix.
! It benchmarked worse than the alternative version
! given below. 
!=========================================================================
!  Braegan Spring  2015
!=====================================================================



if(.false.) then

  if(i==1) p=1
  if(i==2) p=nbpart
  if(i > 2 + n2thread*2) then 
    p = i - (1+n2thread) 
  endif

  if(i > 2 .and. i <= 2 + n2thread*2) then 
    if(mod(i,2) == 1) then
     !Odd threads -- first one could be thread 3 which should work on partition 2
      p = (i+1)/2
    else
    !Even threads -- first one could be thread 4 which should also work on partition 2
      p = i/2
    endif
  endif 
endif

!This attempts to keep two-core spike threads on the same package (in multi package systems)
if(.true.) then 
if(i==1) then
  p=1
endif
if(i==nbpart+n2thread) then 
  p=nbpart 
endif

if(i > 1 + n2thread*2 .and. i < nbpart+n2thread) then 
  p = i - (n2thread) 
endif

if(i > 1 .and. i <= 1 + n2thread*2) then 
  if(mod(i,2) == 1) then
     !Odd threads -- first one could be thread 3 which should work on partition 2
    p = (i+1)/2
  else
    !Even threads -- first one could be thread 4 which should also work on partition 2
    p = 1+i/2
  endif
endif

endif

end subroutine






subroutine spikerl_renumber_threads(spm,n,klu)
implicit none 
include 'f90_noruntime_interface.fi'
integer :: n,klu,i
integer,  dimension(:),pointer :: ajmin
integer, dimension(64) :: spm
character :: trans
complex(kind=(kind(1.0d0))),  dimension(:),pointer :: work
complex(kind=(kind(1.0d0))), parameter :: z_one=(1.0d0,0.0d0), z_zero=(0.0d0,0.0d0)
integer :: infoloc
integer(8),parameter :: fout=6
logical :: test
! Purpose: test if the matrix is banded enough- Reduce the number of threads used
!          by SPIKE if needed
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2015
!=====================================================================


if(spm(22) > 1) then

  call wallocate_1i(Ajmin,spm(20)+1,infoloc)
  call spikerl_calc_size_partitions(spm,Ajmin,n) ! recalculate partition size

  test = .true.
  if (((Ajmin(spm(20))-Ajmin(spm(20)-1))>4*klu) .and. ((Ajmin(3)-Ajmin(2))>4*klu)) test=.false. 

  do while (test)

  spm(22)=spm(22)-1 ! reduce threads 
  call spikerl_calc_num_partitions(spm) ! recalculate partitions #
  if (spm(22)>1) then
    call spikerl_calc_size_partitions(spm,Ajmin,n) ! recalculate partition size
    call spikerl_calc_num_levels(spm)
    if (((Ajmin(spm(20))-Ajmin(spm(20)-1))>4*klu) .and. ((Ajmin(3)-Ajmin(2))>4*klu)) test=.false. 
  else
    test=.false.
    spm(20)=1
  endif
  enddo
  call wdeallocate_1i(Ajmin)
endif


end subroutine
