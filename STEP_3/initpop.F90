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
                 do k = 0, npop-1
                    f(k,i,j) = rhod1(i,j)*w(k)
                 enddo
              enddo
           enddo
        endif

        rhoaver = 0.d0
        do j = 1, ny
           do i = 1, nx
              do k = 0, npop-1
                 rhoaver = rhoaver + f(k,i,j)
              enddo
           enddo
        enddo

        rhoaver = rhoaver / dfloat(nx*ny)
        dinvrho = 1.d0 / rhoaver
        print*,'average density at t=0:', rhoaver
        write(*,*) rhod1(25,25), rhod1(32,32), rhod1(60,60)
        write(*,*) 'rhoaver =', rhoaver
        write(*,*) 'pop -->', f(5,32,32)
        return

        end subroutine initpop
