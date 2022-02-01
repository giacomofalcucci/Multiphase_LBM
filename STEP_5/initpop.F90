! --------------------------------------------------
        subroutine initpop
!---------------------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)
           
        if (dump.eq.1)then
           call resume
        else
           do j = 1, ny
              do i = 1, nx
                 f0(i,j) = rhod1(i,j)*w(0)
                 f1(i,j) = rhod1(i,j)*w(1)
                 f2(i,j) = rhod1(i,j)*w(2)
                 f3(i,j) = rhod1(i,j)*w(3)
                 f4(i,j) = rhod1(i,j)*w(4)
                 f5(i,j) = rhod1(i,j)*w(5)
                 f6(i,j) = rhod1(i,j)*w(6)
                 f7(i,j) = rhod1(i,j)*w(7)
                 f8(i,j) = rhod1(i,j)*w(8)
              enddo
           enddo
        endif

        rhoaver = 0.d0
        do j = 1, ny
           do i = 1, nx
               rhoaver = rhoaver + f0(i,j) +  & 
     &                   f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)+ &
     &                   f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j)    
           enddo
        enddo

        rhoaver = rhoaver / dfloat(nx*ny)
        dinvrho = 1.d0 / rhoaver
        print*,'average density at t=0:', rhoaver
        write(*,*) rhod1(25,25), rhod1(32,32), rhod1(60,60)
        write(*,*) 'rhoaver =', rhoaver
        write(*,*) 'pop -->', f5(32,32)
        return

#ifdef DEBUG
        write(6,*) "Completed subroutine initpop"
#endif


        end subroutine initpop
