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

! set everything to zero
        do j = 1, ny
           do i = 1, nx
              u1(i,j) = 0
              v1(i,j) = 0
              u2(i,j) = 0
              v2(i,j) = 0
           enddo
        enddo


! Case 1 (single bubble)
        if (icond == 1) then 
!                
! for a squared box
!           radius2 = (ny*0.16)**2
! for a rectangular box
           radius2 = (ny*0.25)**2
           do j = 0,ny+1
              do i =0,nx+1
!              
                 u1(i,j) = u0
                 v1(i,j) = v0
                 u2(i,j) = u0
                 v2(i,j) = v0
!            
                 if (((i-nx/2)**2+(j-ny/2)**2).lt.radius2) then
                    rhod1(i,j)= 1.939027909286642
!                    rhod1(i,j)=  &
!                    2.410d0*(1.d0 + 0.01d0* (rand(0) - 0.5d0) * 2.d0) 
                 else
                    rhod1(i,j)=0.1572145251524053
!                    rhod1(i,j)=0.1250d0
                 endif
              enddo
           enddo
        endif

! Case 2 (many bubbles)
        if (icond == 2) then 
!                
! modified to take into account the possibility of localized peturbation
           do j = 1, ny
              do i = 1, nx
                 u1(i,j) = u0
                 v1(i,j) = v0
                 u2(i,j) = u0
                 v2(i,j) = v0
              enddo
           enddo

           do j = 1,ny
              do i =1,nx
                 rhod1(i,j)=rhoin*(1.d0+0.01d0*(rand(0)-0.5d0)*2.d0)
              enddo
           enddo
        endif

! Case 3 (two bubbles colliding)
        if (icond == 3) then
           radius2 = (ny*0.15)**2
!
! left bubble                
           do j = 0,ny+1
              do i = 0,nx/2
                 if (((i-nx/4)**2+(j-ny/2-60)**2).lt.radius2) then
                    rhod1(i,j)= 1.939027909286642
                    u1(i,j) = u0
                    v1(i,j) = v0
                    u2(i,j) = u0
                    v2(i,j) = v0
!
                 else
                    rhod1(i,j)=0.1572145251524053
                 endif
              enddo
           enddo
!
! right bubble                
           do j = 0,ny+1
              do i = nx/2+1, nx+1
                 if (((i-3*nx/4)**2+(j-ny/2)**2).lt.radius2) then
                    rhod1(i,j)= 1.939027909286642
!
                    u1(i,j) = -u0
                    v1(i,j) = -v0
                    u2(i,j) = -u0
                    v2(i,j) = -v0
!
                 else
                    rhod1(i,j)=0.1572145251524053
                 endif
              enddo
           enddo
        endif

        if ((icond .gt. 3).or.(icond.lt.1)) then 
           write(6,*) "ERROR: option non supported"
           stop
        endif

        return
        end subroutine inithydro1

