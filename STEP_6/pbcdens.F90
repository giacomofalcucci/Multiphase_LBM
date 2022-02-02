!-----------------------------------------------------------
        subroutine pbcdens
!-----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

!GA: OpenMP seems not usefull here, too "small" loop probably
!!$OMP PARALLEL DEFAULT(NONE)                            &
!!$OMP PRIVATE(i,j)                                      &
!!$OMP SHARED(nx,ny)                                     &
!!$OMP SHARED(rhod1)                   
!!$OMP DO
        do j = 1, ny
          rhod1(0,j) = rhod1(nx,j)
          rhod1(nx+1,j) = rhod1(1,j)
        enddo

!!$OMP DO
        do i = 1, nx
          rhod1(i,0) = rhod1(i,ny)
          rhod1(i,ny+1) = rhod1(i,1)
        enddo
!!$OMP END PARALLEL

        rhod1(0,0) = rhod1(nx,0)
        rhod1(nx+1,0) = rhod1(nx+1,ny)
        rhod1(nx+1,ny+1) = rhod1(1,ny+1)
        rhod1(0,ny+1) = rhod1(0,1)

        return
        end
