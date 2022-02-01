!---------------------------------------------
        subroutine hydrovarPRE
!----------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)

! Calculation of velocities and pseudopotential
        do j = 0, ny+1
           do i = 0, nx+1
              rhod1(i,j) = fp0(i,j) + &
                           fp1(i,j) + &
                           fp2(i,j) + &
                           fp3(i,j) + &
                           fp4(i,j) + &
                           fp5(i,j) + &
                           fp6(i,j) + &
                           fp7(i,j) + &
                           fp8(i,j)
              rho1 = 1.d0 / rhod1(i,j)        
              u1(i,j) = ( fp1(i,j) - fp3(i,j) + fp5(i,j)            &
                        - fp6(i,j) - fp7(i,j) + fp8(i,j) ) * rho1 
              v1(i,j) = ( fp5(i,j) + fp2(i,j) + fp6(i,j)            &
                        - fp7(i,j) - fp4(i,j) - fp8(i,j) ) * rho1
!                
              u2(i,j) = u1(i,j)
              v2(i,j) = v1(i,j)
!
! (from previous force subroutine)
              psi(i,j) =rhopsi * (1.d0 -exp(- rhod1(i,j) / rhopsi))

           enddo
        enddo

#ifdef DEBUG
        write(6,*) "Completed subroutine hydrovarPRE"
#endif

        end subroutine hydrovarPRE
