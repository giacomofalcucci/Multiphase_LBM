!-----------------------------------------------------------
        subroutine pbc
!-----------------------------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)

        do j = 1, ny
          do k = 0, npop-1
            f(k,0,j) = f(k,nx,j)
            f(k,nx+1,j) = f(k,1,j)
          enddo
        enddo

        do i = 1, nx
          do k = 0, npop-1
            f(k,i,0) = f(k,i,ny)
            f(k,i,ny+1) = f(k,i,1)
          enddo
        enddo

        do k = 0, npop-1
          f(k,0,0) = f(k,nx,0)
          f(k,nx+1,0) = f(k,nx+1,ny)
          f(k,nx+1,ny+1) = f(k,1,ny+1)
          f(k,0,ny+1) = f(k,0,1)
        enddo

        return
        end subroutine pbc
