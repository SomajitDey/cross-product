!About:
!This program describes how the cross-product of two vectors may be
!obtained using an array expression which optimizes for speed.
!
!Author: Somajit Dey <sdphys_rs@caluniv.ac.in> January 2021
!Copyright (C) 2021 Somajit Dey
!License: GNU LGPL v2.1-or-later

program crossproduct
implicit none
integer,dimension(3)::v1,v2,crossprod,expectation
do
    write(*,*)'Give 1st vector:'
    read(*,*)v1
    write(*,*)'Give 2nd vector:'
    read(*,*)v2

!My One-Liner:
    crossprod=cshift((cshift(v1,-1)*v2-v1*cshift(v2,-1)),-1)

    write(*,*)'The cross-product is:'
    write(*,'(3(I0,1X))')crossprod
    
!Conventional Determination of CrossProduct:
    expectation(1)=v1(2)*v2(3)-v1(3)*v2(2)
    expectation(2)=v1(3)*v2(1)-v1(1)*v2(3)
    expectation(3)=v1(1)*v2(2)-v1(2)*v2(1)

!Check if my one-liner gives the same result as expectation    
    if(any(crossprod/=expectation))STOP 'Oops'

    write(*,*)
    write(*,*)'Press ENTER to continue'
    read(*,*)
enddo    
end program crossproduct