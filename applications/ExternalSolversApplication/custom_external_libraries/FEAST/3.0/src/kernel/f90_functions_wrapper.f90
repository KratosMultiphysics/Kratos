!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FILE created by Eric Polizzi 2009 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
         double precision function wdatan2(x,y)
         double precision :: x,y
         wdatan2=atan2(x,y)
         end function wdatan2

         double precision function wdsin(x)
         double precision :: x
         wdsin=sin(x)
         end function wdsin

         double precision function wdcos(x)
         double precision :: x
         wdcos=cos(x)
         end function wdcos

         real function wsatan2(x,y)
         real :: x,y
         wsatan2=atan2(x,y)
         end function wsatan2

         real function wssin(x)
         real :: x
         wssin=sin(x)
         end function wssin

         real function wscos(x)
         real :: x
         wscos=cos(x)
         end function wscos

        subroutine wallocate_3i(A, firstNo, secondNo, thirdNo, alloc_info)  
        integer, dimension(:,:,:),pointer :: A   
        integer :: firstNo   
        integer :: secondNo   
        integer :: thirdNo 
        integer :: alloc_info
        allocate(A(firstNo,secondNo,thirdNo))
        alloc_info=0
        end subroutine wallocate_3i  
   
        subroutine wallocate_2i(A, rowNo, columNo, alloc_info)  
        integer, dimension(:,:),pointer :: A  
        integer :: rowNo  
        integer :: columNo  
        integer :: alloc_info
        allocate(A(rowNo,columNo))
        alloc_info=0
        end subroutine wallocate_2i 
  
        subroutine wallocate_1i(A, element_No, alloc_info)   
        integer, dimension(:),pointer :: A   
        integer :: element_No   
        integer :: alloc_info
        allocate(A(element_No))
        alloc_info=0
        end subroutine wallocate_1i   
  
        subroutine wallocate_3d(A, firstNo, secondNo, thirdNo, alloc_info)   
        double precision, dimension(:,:,:),pointer :: A   
        integer :: firstNo   
        integer :: secondNo 
        integer :: thirdNo   
        integer :: alloc_info
        allocate(A(firstNo,secondNo,thirdNo))
        alloc_info=0
        end subroutine wallocate_3d   
 
        subroutine wallocate_2d(A, rowNo, columNo, alloc_info)  
        double precision, dimension(:,:),pointer :: A  
        integer :: rowNo  
        integer :: columNo  
        integer :: alloc_info
        allocate(A(rowNo,columNo))
        alloc_info=0
        end subroutine wallocate_2d  
  
        subroutine wallocate_1d(A, element_No, alloc_info)   
        double precision, dimension(:),pointer :: A   
        integer :: element_No   
        integer :: alloc_info
        allocate(A(element_No))
        alloc_info=0
        end subroutine wallocate_1d   


        subroutine wallocate_3z(A, firstNo, secondNo, thirdNo, alloc_info)   
        complex(kind=kind(1.0d0)),dimension(:,:,:),pointer :: A   
        integer :: firstNo   
        integer :: secondNo 
        integer :: thirdNo   
        integer :: alloc_info
        allocate(A(firstNo,secondNo,thirdNo))
        alloc_info=0
        end subroutine wallocate_3z   
 
        subroutine wallocate_2z(A, rowNo, columNo, alloc_info)  
        complex(kind=kind(1.0d0)),dimension(:,:),pointer :: A  
        integer :: rowNo  
        integer :: columNo  
        integer :: alloc_info
        allocate(A(rowNo,columNo))
        alloc_info=0
        end subroutine wallocate_2z  
  
        subroutine wallocate_1z(A, element_No, alloc_info)   
        complex(kind=kind(1.0d0)), dimension(:),pointer :: A   
        integer :: element_No   
        integer :: alloc_info
        allocate(A(element_No))
        alloc_info=0
        end subroutine wallocate_1z   

        subroutine wallocate_3s(A, firstNo, secondNo, thirdNo, alloc_info)   
        real, dimension(:,:,:),pointer :: A   
        integer :: firstNo   
        integer :: secondNo 
        integer :: thirdNo   
        integer :: alloc_info
        allocate(A(firstNo,secondNo,thirdNo))
        alloc_info=0
        end subroutine wallocate_3s
 
        subroutine wallocate_2s(A, rowNo, columNo, alloc_info)  
        real, dimension(:,:),pointer :: A  
        integer :: rowNo  
        integer :: columNo  
        integer :: alloc_info
        allocate(A(rowNo,columNo))
        alloc_info=0
        end subroutine wallocate_2s  
  
        subroutine wallocate_1s(A, element_No, alloc_info)   
        real, dimension(:),pointer :: A   
        integer :: element_No   
        integer :: alloc_info
        allocate(A(element_No))
        alloc_info=0
        end subroutine wallocate_1s   


        subroutine wallocate_3c(A, firstNo, secondNo, thirdNo, alloc_info)   
        complex,dimension(:,:,:),pointer :: A   
        integer :: firstNo   
        integer :: secondNo 
        integer :: thirdNo   
        integer :: alloc_info
        allocate(A(firstNo,secondNo,thirdNo))
        alloc_info=0
        end subroutine wallocate_3c   
 
        subroutine wallocate_2c(A, rowNo, columNo, alloc_info)  
        complex,dimension(:,:),pointer :: A  
        integer :: rowNo  
        integer :: columNo  
        integer :: alloc_info
        allocate(A(rowNo,columNo))
        alloc_info=0
        end subroutine wallocate_2c  
  
        subroutine wallocate_1c(A, element_No, alloc_info)   
        complex, dimension(:),pointer :: A   
        integer :: element_No   
        integer :: alloc_info
        allocate(A(element_No))
        alloc_info=0
        end subroutine wallocate_1c   


        subroutine wdeallocate_3i(A)    
        integer, dimension(:,:,:),pointer :: A
        deallocate(A)    
        end subroutine wdeallocate_3i  
 
        subroutine wdeallocate_2i(A)   
        integer, dimension(:,:),pointer :: A 
        deallocate(A)   
        end subroutine wdeallocate_2i   
   
        subroutine wdeallocate_1i(A)    
        integer, dimension(:),pointer :: A
        deallocate(A)      
        end subroutine wdeallocate_1i    
 
        subroutine wdeallocate_3d(A)   
        double precision, dimension(:,:,:),pointer :: A   
        deallocate(A)  
        end subroutine wdeallocate_3d   
 
        subroutine wdeallocate_2d(A)  
        double precision, dimension(:,:),pointer :: A  
        deallocate(A)  
        end subroutine wdeallocate_2d
  
        subroutine wdeallocate_1d(A)   
        double precision, dimension(:),pointer :: A   
        deallocate(A)  
        end subroutine wdeallocate_1d 

        subroutine wdeallocate_3z(A)   
        complex(kind=kind(1.0d0)), dimension(:,:,:),pointer :: A   
        deallocate(A)  
        end subroutine wdeallocate_3z   
 
        subroutine wdeallocate_2z(A)  
        complex(kind=kind(1.0d0)), dimension(:,:),pointer :: A  
        deallocate(A)  
        end subroutine wdeallocate_2z
  
        subroutine wdeallocate_1z(A)   
        complex(kind=kind(1.0d0)), dimension(:),pointer :: A   
        deallocate(A)  
        end subroutine wdeallocate_1z 
 
        subroutine wdeallocate_3s(A)   
        real, dimension(:,:,:),pointer :: A   
        deallocate(A)  
        end subroutine wdeallocate_3s   
 
        subroutine wdeallocate_2s(A)  
        real, dimension(:,:),pointer :: A  
        deallocate(A)  
        end subroutine wdeallocate_2s
  
        subroutine wdeallocate_1s(A)   
        real, dimension(:),pointer :: A   
        deallocate(A)  
        end subroutine wdeallocate_1s 

        subroutine wdeallocate_3c(A)   
        complex, dimension(:,:,:),pointer :: A   
        deallocate(A)  
        end subroutine wdeallocate_3c   
 
        subroutine wdeallocate_2c(A)  
        complex, dimension(:,:),pointer :: A  
        deallocate(A)  
        end subroutine wdeallocate_2c
  
        subroutine wdeallocate_1c(A)   
        complex, dimension(:),pointer :: A   
        deallocate(A)  
        end subroutine wdeallocate_1c 
!!$ 

subroutine wwrite_n(file_id)
 integer(8) :: file_id
write(file_id,*)
end subroutine wwrite_n

subroutine wwrite_t(file_id)
 integer(8) :: file_id
write(file_id,'(3X)',advance='no')
end subroutine wwrite_t



subroutine wwrite_s(file_id,buffer)
 integer(8) :: file_id
 character(len=*) :: buffer
write(file_id,'(A)',advance='no') buffer
end subroutine wwrite_s


subroutine wwrite_i(file_id,buffer)
 integer(8) :: file_id
 integer :: buffer
if (abs(buffer)<10000) then
write(file_id,'(I4)',advance='no') buffer
else
write(file_id,'(I8)',advance='no') buffer
endif

end subroutine wwrite_i


subroutine wwrite_f(file_id,buffer)
 integer(8) :: file_id
 real :: buffer
write(file_id,'(E15.7)',advance='no') buffer
end subroutine wwrite_f


subroutine wwrite_d(file_id,buffer)
 integer(8) :: file_id
 double precision :: buffer
write(file_id,'(E25.16)',advance='no') buffer
end subroutine wwrite_d



