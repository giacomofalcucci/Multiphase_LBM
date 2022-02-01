c----------------------------------------------------------
        subroutine collis

        implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------------------------
        do j = 1, ny
           do i = 1, nx
              do k = 0, npop-1
                 f(k,i,j) = f(k,i,j) * (1.0d0 - omega)
     &                   + omega * feq(k,i,j)
              enddo
           enddo
        enddo

        return 
        end subroutine collis

