!--------------------------------------------------
        subroutine inithydro1
!---------------------------------------------------
! ------- modules
        use storage
        implicit double precision(a-h,o-z)

#ifdef NVFORTRAN        
        real*8 rand 
#else
        real*4 rand ! hack for intel compiler
#endif        

! modified to take into account the possibility of localized peturbation
        do j = 1, ny
           do i = 1, nx
              u1(i,j) = u0
              v1(i,j) = 0.d0
              u2(i,j) = u0
              v2(i,j) = 0.d0
           enddo
        enddo

! Case 1 (single bubble)
        if (icond == 1) then 
           do j = 0,ny+1
              do i =0,nx+1
                 if(((i-nx/2)**2+(j-ny/2)**2).lt.400) then
                    rhod1(i,j)=  &
                    2.410d0*(1.d0 + 0.01d0 * (rand(0) - 0.5d0) * 2.d0) 
                 else
                    rhod1(i,j)=0.1250d0
                 endif
              enddo
           enddo
        endif

! Case 2 (many bubbles)
        if (icond == 2) then 
           do j = 1,ny
              do i =1,nx
                 rhod1(i,j)=rhoin*(1.d0+0.01d0*(rand(0)-0.5d0)*2.d0)
              enddo
           enddo
        endif

        if ((icond .gt. 2).or.(icond.lt.1)) then 
           write(6,*) "ERROR: option non supported"
           stop
        endif

        return
        end subroutine inithydro1

