! --------------------------------------------------
        subroutine move
!---------------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)

        do j = 1,ny
           do i = 1, nx
              fp(0,i,j) = f(0,i,j)
              fp(1,i,j) = f(1,i-1,j)
              fp(2,i,j) = f(2,i,j-1)
              fp(3,i,j) = f(3,i+1,j)
              fp(4,i,j) = f(4,i,j+1)
              fp(5,i,j) = f(5,i-1,j-1)
              fp(6,i,j) = f(6,i+1,j-1)
              fp(7,i,j) = f(7,i+1,j+1)
              fp(8,i,j) = f(8,i-1,j+1)
           enddo
        enddo

! border fix

        j = 0
        do i = 0, nx+1 
           fp(0,i,j)=f(0,i,j)
           fp(1,i,j)=f(1,i,j)
           fp(2,i,j)=f(2,i,j)
           fp(3,i,j)=f(3,i,j)
           fp(4,i,j)=f(4,i,j)
           fp(5,i,j)=f(5,i,j)
           fp(6,i,j)=f(6,i,j)
           fp(7,i,j)=f(7,i,j)
           fp(8,i,j)=f(8,i,j)
        enddo
!
        j = ny+1
        do i = 0, nx+1 
           fp(0,i,j)=f(0,i,j)
           fp(1,i,j)=f(1,i,j)
           fp(2,i,j)=f(2,i,j)
           fp(3,i,j)=f(3,i,j)
           fp(4,i,j)=f(4,i,j)
           fp(5,i,j)=f(5,i,j)
           fp(6,i,j)=f(6,i,j)
           fp(7,i,j)=f(7,i,j)
           fp(8,i,j)=f(8,i,j)
        enddo
!
        i = 0
        do j = 0, ny+1
           fp(0,i,j)=f(0,i,j)
           fp(1,i,j)=f(1,i,j)
           fp(2,i,j)=f(2,i,j)
           fp(3,i,j)=f(3,i,j)
           fp(4,i,j)=f(4,i,j)
           fp(5,i,j)=f(5,i,j)
           fp(6,i,j)=f(6,i,j)
           fp(7,i,j)=f(7,i,j)
           fp(8,i,j)=f(8,i,j)
        enddo
!
        i = nx+1
        do j = 0, ny+1
           fp(0,i,j)=f(0,i,j)
           fp(1,i,j)=f(1,i,j)
           fp(2,i,j)=f(2,i,j)
           fp(3,i,j)=f(3,i,j)
           fp(4,i,j)=f(4,i,j)
           fp(5,i,j)=f(5,i,j)
           fp(6,i,j)=f(6,i,j)
           fp(7,i,j)=f(7,i,j)
           fp(8,i,j)=f(8,i,j)
        enddo
!
#ifdef DEBUG
        write(6,*) "Completed subroutine move"
#endif
        end subroutine move
