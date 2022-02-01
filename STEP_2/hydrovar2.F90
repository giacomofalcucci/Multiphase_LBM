!---------------------------------------------
        subroutine hydrovar2
!----------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

! Calculation of velocities and pseudopotential
        do j = 0, ny+1
           do i = 0, nx+1
              rhod2(i,j) = f(0,i,j) + &
                           f(1,i,j) + &
                           f(2,i,j) + &
                           f(3,i,j) + &
                           f(4,i,j) + &
                           f(5,i,j) + &
                           f(6,i,j) + &
                           f(7,i,j) + &
                           f(8,i,j)
              rho2 = 1.d0 / rhod2(i,j)
              u2(i,j) = ( f(1,i,j) - f(3,i,j) + f(5,i,j)   &
                        - f(6,i,j) - f(7,i,j) + f(8,i,j) ) * rho2
              v2(i,j) = ( f(5,i,j) + f(2,i,j) + f(6,i,j) &
                        - f(7,i,j) - f(4,i,j) - f(8,i,j) ) * rho2
           enddo
        enddo

        return
        end subroutine hydrovar2

