!*******************************************************************************
!  SOMAR - Stratified Ocean Model with Adaptive Refinement
!  Developed by Ed Santilli & Alberto Scotti
!  Copyright (C) 2014 University of North Carolina at Chapel Hill
!
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
!  USA
!
!  For up-to-date contact information, please visit the repository homepage,
!  https://github.com/somarhub.
!*******************************************************************************


! ------------------------------------------------------------------------------
! Solves a symmetric, tridiagonal system.
! The matrix indices are the COL number
!   a - sub-diagonal (means it is the diagonal below the main diagonal)
!   b - the main diagonal
!   c - sup-diagonal (means it is the diagonal above the main diagonal)
!   d - right part
!   x - the answer
!   n - number of equations
! ------------------------------------------------------------------------------
      subroutine solve_tridiag(a,b,c,d,x,n)
        implicit none

        integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: b,d
        real(8),dimension(n-1),intent(in) :: a,c
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: cp,dp
        real(8) :: m
        integer i

! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n-1
           m = b(i)-cp(i-1)*a(i-1)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i-1))/m
         enddo

! Last iter a special case to avoid calling c(n)
         i = n
         m = b(i)-cp(i-1)*a(i-1)
         dp(i) = (d(i)-dp(i-1)*a(i-1))/m

! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do

      end subroutine solve_tridiag


! ------------------------------------------------------------------------------
! Solves a symmetric, tridiagonal system.
! The matrix indices are the COL number
!   a - sub-diagonal (means it is the diagonal below the main diagonal)
!   b - the main diagonal
!   d - right part
!   x - the answer
!   n - number of equations
! ------------------------------------------------------------------------------
      subroutine solve_symtridiag(a,b,d,x,n)
        implicit none

        integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: b,d
        real(8),dimension(n-1),intent(in) :: a
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: cp,dp
        real(8) :: m
        integer i

! initialize c-prime and d-prime
        cp(1) = a(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n-1
           m = b(i)-cp(i-1)*a(i-1)
           cp(i) = a(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i-1))/m
         enddo

! Last iter a special case to avoid calling a(n)
         i = n
         m = b(i)-cp(i-1)*a(i-1)
         dp(i) = (d(i)-dp(i-1)*a(i-1))/m

! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do

      end subroutine solve_symtridiag
