!-------------------------------------------------------------
        subroutine mbc
!-------------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)
! EAST case
        do j = 1,ny
           f(1,0,j) = f(1,nx,j)
           f(5,0,j) = f(5,nx,j)
           f(8,0,j) = f(8,nx,j)
        enddo

! WEST case
        do j = 1,ny
           f(3,nx+1,j) = f(3,1,j)
           f(6,nx+1,j) = f(6,1,j)
           f(7,nx+1,j) = f(7,1,j)
        enddo

! NORTH case
        do i = 1,nx
           f(4,i,ny+1) = f(2,i,ny)
           f(8,i,ny+1) = f(6,i,ny)
           f(7,i,ny+1) = f(5,i,ny)
        enddo

! SOUTH case
        do i = 1,nx
           f(2,i,0) = f(4,i,1)
           f(6,i,0) = f(8,i,1)
           f(5,i,0) = f(7,i,1)
        enddo

        return
        end subroutine mbc
