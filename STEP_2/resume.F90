!----------------------------------------------------------
        subroutine resume 
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        write(6,*)'################################'
        write(6,*)'reading  popolations ....'
        write(6,*)'################################'

        do j=0,ny+1
           do i=0,nx+1
              do k=0,npop-1
                 read(111)f(k,i,j)
              enddo
           enddo
        enddo

        do j=1,ny
           do i=1,nx
              read(113)u1(i,j),v1(i,j)
           enddo
        enddo

        do j=-1,ny+2
           do i=-1,nx+2
              read(112)rhod1(i,j)
           enddo
        enddo

        return
        end subroutine resume

