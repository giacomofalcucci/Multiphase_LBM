c----------------------------------------------------------
        subroutine collis

        implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------------------------
        do j = 1, ny
           do i = 1, nx
              f(0,i,j) = f(0,i,j) * (1.0d0 - omega) + omega * feq(0,i,j)
              f(1,i,j) = f(1,i,j) * (1.0d0 - omega) + omega * feq(1,i,j)
              f(2,i,j) = f(2,i,j) * (1.0d0 - omega) + omega * feq(2,i,j)
              f(3,i,j) = f(3,i,j) * (1.0d0 - omega) + omega * feq(3,i,j)
              f(4,i,j) = f(4,i,j) * (1.0d0 - omega) + omega * feq(4,i,j)
              f(5,i,j) = f(5,i,j) * (1.0d0 - omega) + omega * feq(5,i,j)
              f(6,i,j) = f(6,i,j) * (1.0d0 - omega) + omega * feq(6,i,j)
              f(7,i,j) = f(7,i,j) * (1.0d0 - omega) + omega * feq(7,i,j)
              f(8,i,j) = f(8,i,j) * (1.0d0 - omega) + omega * feq(8,i,j)
           enddo
        enddo

        return 
        end subroutine collis

