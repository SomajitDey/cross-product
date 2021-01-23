!About:
!This program introduces a 3x3 matrix inverse function that uses the 
!array expression for cross-product. See matrix_inv below
!
!Author: Somajit Dey <sdphys_rs@caluniv.ac.in> January 2021
!Copyright (C) 2021 Somajit Dey
!License: GNU LGPL v2.1-or-later

program matrix_inversion
implicit none
real(8),dimension(3,3)::matrix

write(*,*)'Input matrix:'
call input(matrix)
write(*,*)
write(*,*)'Inverse:'
call output(matrix_inv(matrix))
write(*,*)
write(*,*)'Product'
call output(matmul(matrix,matrix_inv(matrix)),tolerance=1E-8)

contains

function matrix_inv(matrix)
real(8),dimension(3,3),intent(in)::matrix
real(8),dimension(3,3)::matrix_inv
real(8),dimension(3,3)::tmp1,tmp2,tmp
real(8)::dotproduct
real(8),parameter::tol=2*tiny(tol)
integer::i
logical::is_invertible=.true.

tmp1=cshift(matrix,-1,2)
tmp2=cshift(matrix,1,2)
tmp=cshift((cshift(tmp1,-1,1)*tmp2-tmp1*cshift(tmp2,-1,1)),-1,1)

do concurrent (i=1:3)
    dotproduct=dot_product(tmp(:,i),matrix(:,i))
    if(abs(dotproduct)>tol)then
        tmp(:,i)=tmp(:,i)/dotproduct
    else
        is_invertible=.false.
    endif
enddo
if(.NOT.is_invertible)STOP 'Not invertible'

matrix_inv=transpose(tmp)
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