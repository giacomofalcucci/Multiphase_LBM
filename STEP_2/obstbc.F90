!--------------------------------------------------------
        subroutine obstbc
!--------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        k = nx / 4

        do j = ny/2-nobst/2+1,ny/2+nobst/2
           f(1,k+1,j) = f(3,k+1,j)
           f(3,k  ,j) = f(1,k  ,j)
        enddo

        do j = ny/2-nobst/2,ny/2+nobst/2+1
           f(5,k+1,j) = f(7,k+1,j)
           f(8,k+1,j) = f(6,k+1,j)
           f(7,k,  j) = f(5,k,  j)
           f(6,k,  j) = f(8,k,  j)
        enddo

        return
        end subroutine obstbc
