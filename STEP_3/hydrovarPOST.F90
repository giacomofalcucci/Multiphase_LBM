!---------------------------------------------
        subroutine hydrovarPOST
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

! computing mean velocity (pre/post collision)
              u1(i,j)=(u1(i,j)+u2(i,j))*0.5d0
              v1(i,j)=(v1(i,j)+v2(i,j))*0.5d0

           enddo
        enddo
#ifdef DEBUG
        write(6,*) "Completed subroutine hydrovarPOST"
#endif

        end subroutine hydrovarPOST

