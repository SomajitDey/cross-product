!About:
!This program introduces a 3x3 matrix inverse function that uses the 
!array expression for cross-product. See matrix_inv below
!
!Author: Somajit Dey <sdphys_rs@caluniv.ac.in> January 2021
!Copyright (C) 2021 Somajit Dey
!License: GNU LGPL v2.1-or-later

program matrix_inversion
implicit none
real(8),dimension(3,3)::matrix,inverted
real(8)::determinant
logical::is_invertible

write(*,*)'Input matrix:'
call input(matrix)
write(*,*)

inverted=matrix_inv(matrix,determinant,is_invertible)

write(*,*)'Determinant:'
write(*,*)determinant
write(*,*)

if(.NOT.is_invertible)STOP 'Not invertible'

write(*,*)'Inverse:'
call output(inverted)
write(*,*)

write(*,*)'Product:'
call output(matmul(matrix,inverted),tolerance=1E-8)

contains

function matrix_inv(matrix,determinant,is_invertible)
real(8),dimension(3,3),intent(in)::matrix
real(8),dimension(3,3)::matrix_inv
real(8),intent(out),optional::determinant
logical,intent(out),optional::is_invertible
real(8),dimension(3,3)::tmp1,tmp2,tmp
real(8)::dotproduct
real(8),parameter::tol=2*tiny(tol) !Tolerance
tmp1=cshift(matrix,1,2)
tmp2=cshift(matrix,-1,2)
tmp=cshift((cshift(tmp1,-1,1)*tmp2-tmp1*cshift(tmp2,-1,1)),-1,1)
dotproduct=dot_product(tmp(:,1),matrix(:,1)) !Could be 2 or 3 instead of 1
if(abs(dotproduct)>tol)then
    matrix_inv=transpose(tmp)/dotproduct
else
    matrix_inv=0.0 !To return a null-matrix in case of non-invertibility
endif

!Optional extras follow:
if(present(determinant))determinant=dotproduct
if(present(is_invertible))then
    if(abs(dotproduct)>tol)then
        is_invertible=.true.
    else
        is_invertible=.false.
    endif
endif
end function matrix_inv

subroutine input(matrix)
real(8),dimension(3,3),intent(out)::matrix
integer::row,col
do row=1,3
    read(*,*)(matrix(row,col),col=1,3)
enddo
end subroutine input

subroutine output(matrix,tolerance)
real(8),dimension(3,3)::matrix
integer::row,col
real,intent(in),optional::tolerance
if(present(tolerance))where(abs(matrix)<tolerance)matrix=0.0
do row=1,3
    write(*,*)(real(matrix(row,col)),col=1,3)
enddo
end subroutine output

end program matrix_inversion