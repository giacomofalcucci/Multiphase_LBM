c --------------------------------------------------
        subroutine move
c----------------------------------------------
        implicit double precision(a-h,o-z)
        include 'muphase.par'
c---------------------------------------------
        do j = 0,ny+1
           do i = 0, nx+1 
              fp(1,i,j)=f(1,i,j)
              fp(2,i,j)=f(2,i,j)
              fp(3,i,j)=f(3,i,j)
              fp(4,i,j)=f(4,i,j)
              fp(5,i,j)=f(5,i,j)
              fp(6,i,j)=f(6,i,j)
              fp(7,i,j)=f(7,i,j)
              fp(8,i,j)=f(8,i,j)
           enddo
        enddo

        do j = 1,ny
           do i = 1, nx
              f(1,i,j) = fp(1,i-1,j)
              f(2,i,j) = fp(2,i,j-1)
              f(3,i,j) = fp(3,i+1,j)
              f(4,i,j) = fp(4,i,j+1)
              f(5,i,j) = fp(5,i-1,j-1)
              f(6,i,j) = fp(6,i+1,j-1)
              f(7,i,j) = fp(7,i+1,j+1)
              f(8,i,j) = fp(8,i-1,j+1)
           enddo
        enddo

        return
        end subroutine move
