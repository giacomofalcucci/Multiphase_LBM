!---------------------------------------------
        subroutine hydrovar1
!----------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)

! Calculation of velocities and pseudopotential
        do j = 0, ny+1
           do i = 0, nx+1
              rhod1(i,j) = f(0,i,j) + &
                           f(1,i,j) + &
                           f(2,i,j) + &
                           f(3,i,j) + &
                           f(4,i,j) + &
                           f(5,i,j) + &
                           f(6,i,j) + &
                           f(7,i,j) + &
                           f(8,i,j)
              rho1 = 1.d0 / rhod1(i,j)        
              u1(i,j) = ( f(1,i,j) - f(3,i,j) + f(5,i,j)            &
                        - f(6,i,j) - f(7,i,j) + f(8,i,j) ) * rho1 
              v1(i,j) = ( f(5,i,j) + f(2,i,j) + f(6,i,j)            &
                        - f(7,i,j) - f(4,i,j) - f(8,i,j) ) * rho1
           enddo
        enddo

        return
        end subroutine hydrovar1
