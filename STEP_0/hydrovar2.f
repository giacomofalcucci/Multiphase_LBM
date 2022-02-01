c---------------------------------------------
        subroutine hydrovar2

        implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------

        do j = 1, ny
           do i = 1, nx
              rhod2(i,j) = 0.d0
              do k = 0, npop-1
                 rhod2(i,j) = rhod2(i,j) + f(k,i,j)
              enddo
           enddo
        enddo

!       write(*,*) '!DEBUG2! pop', f(5,32,32)

c Calculation of velocities and pseudopotential
        do j = 1, ny
           do i = 1, nx
              rho2 = 1.d0 / rhod2(i,j)
              u2(i,j) = ( f(1,i,j) - f(3,i,j) + f(5,i,j)  
     &                  - f(6,i,j) - f(7,i,j) + f(8,i,j) ) * rho2
              v2(i,j) = ( f(5,i,j) + f(2,i,j) + f(6,i,j)
     &                  - f(7,i,j) - f(4,i,j) - f(8,i,j) ) * rho2
           enddo
        enddo

        return
        end subroutine hydrovar2

