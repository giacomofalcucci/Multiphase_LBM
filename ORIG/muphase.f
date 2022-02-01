c ================================
	program bgk2D
c ================================
	implicit double precision(a-h,o-z)
	include 'muphase.par'

c ------- input parameters

        omega=1.d0 

	call input
	
c -------initialisation

	call inithydro1
        call initpop

c ------- MAIN LOOP

    	do 10 istep = 1,nsteps
           call pbcdens
           call pbc
           call move
           call hydrovar1
           call hydrovar2
           if (iforce) then
              call force
           endif 
           call equili
           call collis
           call hydrovar1
           call media

	   if (mod(istep,ndiag).eq.0) then
	      call out0d
            call energy
            call calcsurf
	   endif

	   if (mod(istep,nout).eq.0) then
	      call out1d(frce)
              call out2d(frce)
	   endif

           if(mod(istep,100).eq.0) then
            write(6,*) 'completed time step ',istep
           endif

           if(mod(istep,20000).eq.0) then
             call savepop
           endif

!
!           if(mod(istep,nsteps).eq.0) then
!             call laplace    !  it needs only for Case 1
!           endif

c-------- end of main loop
            
10	continue


	end
c--------------------------------------------------
	subroutine input

	implicit double precision(a-h,o-z)
	include 'muphase.par'
c---------------------------------------------------
        open(5,file='muphase.inp')

	print*,' Number of steps'
        read(5,*)nsteps

	print*,' Number of steps between printing profile'
        read(5,*)nout

	print*,' Number of steps between performing diagnostics'
	read(5,*)ndiag

        print*,'dt'
        read(5,*)dt

        print*,'viscosity'
        read(5,*)visc

        print*,' Coupling gnn'
        read(5,*)gnn

	print*,'Insert rhopsi'
        read(5,*)rhopsi

        print*,'Applied force  (.TRUE. or .FALSE.) ?'
	read(5,*)iforce 

	print*,' Initial density'
	read(5,*)rhoin
     
        print*,' Initial X velocity component'
        read(5,*)u0

        print*,' Initial Y velocity component'
        read(5,*)v0

	print*,' Final velocity for the Poise force'
	read(5,*) uf

	print*,' Linear obstacle ?'
	read(5,*) iobst

	if (iobst) then
	    print*,' Length of the obstacle (multiple of 2)'
	    read(5,*) nobst
	endif

        print*,' File for output: 5 chars'
        read(5,'(A)')fileout

        print*,' read populations dump (0 or 1)'
        read(5,*)dump

        close(5)

        open(11,file=fileout//'.rho_u_v_1D')
        open(51,file=fileout//'.ruv2d')
        open(111,file='dump_pop',status='unknown',
     &        form='unformatted')

        open(112,file='dump_rhod1',status='unknown',
     &        form='unformatted')

        open(113,file='dump_u1_v1',status='unknown',
     &        form='unformatted')

        open(55,file=fileout//'.forces_1D')

        print*,'*******************************************************'
        print*,'         Lattice BGK model, 2D with 9 velocities'
        print*,'             multiphse code with Shan-Chen EoS'
        print*,'*******************************************************'
        print*,'  developed and released by Prof. GIACOMO FALCUCCI, PhD'
        print*,'               CFD School @ CINECA'
        print*,'                         2021'
        call sleep(1)
	print*,'Number of cells :',nx,'*',ny
	print*,'Nsteps :',nsteps
	print*,'Relaxation frequency :',omega
	print*,'Coupling gnn :',gnn
	print*,'Coupling gnnn :',gnnn
	print*,'Applied force :',iforce
	print*,'Initial velocity for this Poiseuille force :',u0
	if (iobst) then
	    print*,' Linear Obstacle with length :',nobst
	endif
	write(6,*)'Output file :',fileout

c constants

	cs2  = 1.0d0 / 3.0d0
	cs22 = 2.0d0 * cs2
	cssq = 2.0d0 / 9.0d0

! input weights and storage in w(o:npop-1) array

!	w0 = 4.0d0 / 9.0d0
!	w1 = 1.0d0 / 9.0d0
!	w2 = 1.0d0 / 36.0d0
!	w(0) = w0
!	do i = 1, 4
!           w(i) = w1
!           w(i+4) = w2
!        end do

! ATTENZIONE: i pesi di Mauro sono a meno di cs2 !!!!!!!!!!!!!!!!!!!

        w1 =4.d0/21.d0/3.d0
        w2 =4.d0/45.d0/3.d0

        w4 =1.d0/60.d0/3.d0
        w5 =2.d0/315.d0/3.d0
        w8 =1.d0/5040.d0/3.d0

        w0=1.d0-(4.d0*w1+4.d0*w2+4.d0*w4+8.d0*w5+4.d0*w8)

! calcolo delle distanze reticolari al quadrato per la subroutine ENERGY        

        c1_2=1.d0
        c2_2=2.d0

        c4_2=4.d0
        c5_2=5.d0
        c8_2=8.d0 
        


        w(0) = 4.d0/9.d0
        do i = 1, 4
           w(i) = 1.d0/9.d0
           w(i+4) = 1.d0/36.d0
        end do

!       w1=1.d0/9.d0
!       w2=1.d0/36.d0
!       w0=4.d0/9.d0 
        
        



c reduced density

	den = rhoin/float(npop) 


c scaling

        dx = dt

c calculation of omega

        omega = 1.d0/(3.*visc*(dt*dt)/(dx*dx) + 0.5*dt)

!	visc = (1.0d0 / omega - 0.5d0) * cs2
	print*,' Viscosity :',visc,omega,w(0)

c calculation of the constant applied force

	fpois = 8.0d0 * visc * uf / dfloat(ny) / dfloat(ny)
        fpois = rhoin*fpois/6.  ! # of biased populations
	print*,' Intensity of the applied force ',fpois
	
	return
	end
c--------------------------------------------------
	subroutine inithydro1
	
	implicit double precision(a-h,o-z)
	include 'muphase.par'
        real*4 rand  ! hack for intel compiler 
c---------------------------------------------------
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

!      do j = 0,ny+1
!       do i =0,nx+1
!
!        if(((i-nx/2)**2+(j-ny/2)**2).lt.400) then
!          rhod1(i,j)=2.410d0*(1.d0 + 0.01d0 * (rand(0) - 0.5d0) * 2.d0) 
!         else
!          rhod1(i,j)=0.1250d0
!        endif
! 
!       enddo
!      enddo

! Case 2 (many bubbles)

        do j = 1,ny
           do i =1,nx
              rhod1(i,j)=rhoin*(1.d0 + 0.01d0*(rand(0)-0.5d0)*2.d0)
           enddo
        enddo
        write(6,*)rand(0)
!        stop

	return
	end

c --------------------------------------------------
	subroutine initpop
	
	implicit double precision(a-h,o-z)
	include 'muphase.par'
c---------------------------------------------------
           
        if(dump.eq.1)then
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
        write(*,*) 'popolazione', f(5,32,32)
        return
        end
c----------------------------------------------
	subroutine move
c----------------------------------------------
	implicit double precision(a-h,o-z)
	include 'muphase.par'
c---------------------------------------------
        do j = 0,ny+1
           do i = 0, nx+1 
              do k=0,npop-1
                 fp(k,i,j)=f(k,i,j)
              enddo
           enddo
        enddo


	do j = 1,ny
	   do i = 1, nx

              f(2,i,j) = fp(2,i,j-1)
              f(6,i,j) = fp(6,i+1,j-1)
	      f(1,i,j) = fp(1,i-1,j)
	      f(5,i,j) = fp(5,i-1,j-1)
   	      f(4,i,j) = fp(4,i,j+1)
	      f(8,i,j) = fp(8,i-1,j+1)
	      f(3,i,j) = fp(3,i+1,j)
	      f(7,i,j) = fp(7,i+1,j+1)
	      f(0,i,j) = fp(0,i,j)
           enddo
	enddo

	return
	end
c---------------------------------------------
	subroutine hydrovar1

	implicit double precision(a-h,o-z)
	include 'muphase.par'
c----------------------------------------

	do j = 1, ny
	  do i = 1, nx
            rhod1(i,j) = 0.d0
            do k = 0, npop-1
              rhod1(i,j) = rhod1(i,j) + f(k,i,j)
            enddo
	  enddo
	enddo

!        write(*,*) '!DEBUG! popolazione', f(5,32,32)

c Calculation of velocities and pseudopotential

	do j = 1, ny
	  do i = 1, nx

           rho1 = 1.d0 / rhod1(i,j)        

	   u1(i,j) = ( f(1,i,j) - f(3,i,j) + f(5,i,j) - 
     &              f(6,i,j) - f(7,i,j) + f(8,i,j) ) * rho1 

	   v1(i,j) = ( f(5,i,j) + f(2,i,j) + f(6,i,j)
     &              - f(7,i,j) - f(4,i,j) - f(8,i,j) ) * rho1

     
          enddo
	enddo
		
	return
	end

c---------------------------------------------
        subroutine hydrovar2

        implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------

        do j = 1, ny
          do i = 1, nx
            rhod2(i,j) = 0.d0
            do k = 0, npop-1
              rhod2(i,j) = rhod2(i,j) + f(k,i,j)
            enddo
          enddo
        enddo

!        write(*,*) '!DEBUG2! popolazione', f(5,32,32)


c Calculation of velocities and pseudopotential

        do j = 1, ny
          do i = 1, nx

           rho2 = 1.d0 / rhod2(i,j)

           u2(i,j) = ( f(1,i,j) - f(3,i,j) + f(5,i,j) -
     &              f(6,i,j) - f(7,i,j) + f(8,i,j) ) * rho2

           v2(i,j) = ( f(5,i,j) + f(2,i,j) + f(6,i,j)
     &              - f(7,i,j) - f(4,i,j) - f(8,i,j) ) * rho2

          enddo
        enddo

        return
        end

c-------------------------------------------------
	subroutine equili

	implicit double precision(a-h,o-z)
	include 'muphase.par'

c-------------------------------------------------
	do j = 1, ny
	   do i = 1, nx
	     usq = u1(i,j) * u1(i,j) 
	     vsq = v1(i,j) * v1(i,j)
	     sumsq = (usq + vsq) / cs22
	     sumsq2 = sumsq * (1.0d0 - cs2) / cs2
	     u22 = usq / cssq 
             v22 = vsq / cssq
	     ui = u1(i,j) / cs2
	     vi = v1(i,j) / cs2
	     uv = ui * vi
             rhoij = rhod1(i,j)

	     

	     feq(0,i,j) = (4.d0/9.d0)*rhoij*(1.0d0 - sumsq)

	     feq(1,i,j) = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + u22 + ui)
	     feq(2,i,j) = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + v22 + vi)
	     feq(3,i,j) = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + u22 - ui)
	     feq(4,i,j) = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + v22 - vi)

	     feq(5,i,j) = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 + ui + vi 
     &		 + uv)
	     feq(6,i,j) = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 - ui + vi 
     &	 - uv)
	     feq(7,i,j) = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 - ui - vi 
     &	 + uv)
	     feq(8,i,j) = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 + ui - vi
     &   - uv)
	
	   enddo
	enddo

	return
	end
c----------------------------------------------------------
	subroutine collis

	implicit double precision(a-h,o-z)
	include 'muphase.par'

c----------------------------------------------------------
!        omt=omega * dt  
	 do j = 1, ny
	   do i = 1, nx
	     do k = 0, npop-1
	       f(k,i,j) = f(k,i,j) * (1.0d0 - omega)
     &                   + omega * feq(k,i,j)
	     end do
	   end do
	 end do

	return 
	end

c===========================================
	subroutine media
c ==========================================
	implicit double precision(a-h,o-z)
	include 'muphase.par'
c--------------------------------------------------------

          do i=1,nx
            do j=1,ny
 
              u1(i,j)=(u1(i,j)+u2(i,j))*0.5d0
              v1(i,j)=(v1(i,j)+v2(i,j))*0.5d0
   
            end do
          end do
	return
	end
c ==========================================
	subroutine force
c ==========================================
	implicit double precision(a-h,o-z)
	include 'muphase.par'
c--------------------------------------------------------

        frce = fpois   ! constant external force  
        do i=-1,nx+2
         do j=-1,ny+2    ! ho levato le BC per psi!!!
     
            psi(i,j) =rhopsi * ( 1.d0 - exp(- rhod1(i,j) / rhopsi))
             
         enddo
        enddo      


        wnn=w1/cs2
        wnnn=w2/cs2

        ww_n  = ww1/cs2
        ww_nn = ww2/cs2


 	do j = 1, ny
          do i = 1, nx
c non local Lennard-Jones (may put a density functional, not density itself)

         fnnx  = psi(i,j)*(psi(i+1,j)-psi(i-1,j))
         fnny  = psi(i,j)*(psi(i,j+1)-psi(i,j-1))
         fnnnx = psi(i,j)*((psi(i+1,j-1)+psi(i+1,j+1)
     &                  -psi(i-1,j-1) - psi(i-1,j+1)))
         fnnny = psi(i,j)*((psi(i+1,j+1) + psi(i-1,j+1)
     &                  - psi(i-1,j-1) - psi(i+1,j-1)))


! interactions

         f2x =-(gnnn * fnnx * w1 + gnnn * fnnnx * w2)-  ! prima era solo -(gnn * fnnx * w1 + gnn * fnnnx * w2)
     &         (gnn * fnnx*1./9.+ gnn * fnnnx *1./36.)  ! aggiunte per completare la simmetria del reticolo (Chib)
         f2y =-(gnnn * fnny * w1 + gnnn * fnnny * w2)-  ! prima era solo -(gnn * fnnx * w1 + gnn * fnnnx * w2)
     &         (gnn * fnny*1./9.+ gnn * fnnny *1./36.)  ! aggiunte per completare la simmetria del reticolo (Chib) 


        force_x(i,j)=f2x  !+f2xn
        force_y(i,j)=f2y  !+f2yn

        f_x(i,j) =f2x
        f_2x(i,j)=f2xn

        f_y(i,j) =f2y
        f_2y(i,j)=f2yn



       end do
       end do


! SHIFT EQUILIBRIO

       do i=1,nx
          do j=1,ny

           u1(i,j)=u1(i,j)+force_x(i,j)/(omega*rhod1(i,j))
       
           v1(i,j)=v1(i,j)+force_y(i,j)/(omega*rhod1(i,j)) 

          end do
       end do 


        return
        end
c =========================
	subroutine pbc
c =========================
	implicit double precision(a-h,o-z)
	include 'muphase.par'
c-----------------------------------------------------------

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
	end
c =========================
	subroutine pbcdens
c =========================
	implicit double precision(a-h,o-z)
	include 'muphase.par'
c-----------------------------------------------------------

        do j = 1, ny
          rhod1(0,j) = rhod1(nx,j)
          rhod1(nx+1,j) = rhod1(1,j)
        enddo

        do i = 1, nx
          rhod1(i,0) = rhod1(i,ny)
          rhod1(i,ny+1) = rhod1(i,1)
        enddo

        rhod1(0,0) = rhod1(nx,0)
        rhod1(nx+1,0) = rhod1(nx+1,ny)
        rhod1(nx+1,ny+1) = rhod1(1,ny+1)
        rhod1(0,ny+1) = rhod1(0,1)

	return
	end
c-------------------------------------------------------------
	subroutine mbc
	
	implicit double precision(a-h,o-z)
	include 'muphase.par'
c-------------------------------------------------------------
c EAST case

	 do j = 1,ny
           f(1,0,j) = f(1,nx,j)
           f(5,0,j) = f(5,nx,j)
           f(8,0,j) = f(8,nx,j)
        enddo

c WEST case

	do j = 1,ny
           f(3,nx+1,j) = f(3,1,j)
           f(6,nx+1,j) = f(6,1,j)
           f(7,nx+1,j) = f(7,1,j)
        enddo

c NORTH case

	do i = 1,nx
	   f(4,i,ny+1) = f(2,i,ny)
	   f(8,i,ny+1) = f(6,i,ny)
	   f(7,i,ny+1) = f(5,i,ny)
	enddo

c SOUTH case

	do i = 1,nx
           f(2,i,0) = f(4,i,1)
           f(6,i,0) = f(8,i,1)
           f(5,i,0) = f(7,i,1)
        enddo

	return
	end
c ==========================
	subroutine obstbc
c ==========================
	implicit double precision(a-h,o-z)
	include 'muphase.par'
c--------------------------------------------------------
	k = nx / 4	
	
        do j = ny/2-nobst/2+1,ny/2+nobst/2
           f(1,k+1,j) = f(3,k+1,j)
	   f(3,k  ,j) = f(1,k,j)
        enddo

        do j = ny/2-nobst/2,ny/2+nobst/2+1
	   f(5,k+1,j) = f(7,k+1,j)
	   f(8,k+1,j) = f(6,k+1,j)
	   f(7,k,  j) = f(5,k,  j)
	   f(6,k,  j) = f(8,k,  j)
        enddo

	return
	end
c ===================================
	subroutine out1d(frce)
c ===================================
	implicit double precision(a-h,o-z)
	include 'muphase.par'
c----------------------------------------------------------

c transverse profiles
	do j = 1,ny
	  write(11,*) j,rhod1(nx/2,j),u1(nx/2,j),v1(nx/2,j)
          write(55,*) j,force_x(nx/2,j),force_y(nx/2,j)
	enddo
	write(11,'(bn)')
        write(55,'(bn)')

	return
	end

c =================================
	subroutine out2d(frce)
c =================================
	implicit double precision(a-h,o-z)
	include 'muphase.par'
c----------------------------------------------------------
	do j = 1, ny
	 do i = 1, nx
          write(51,99) i,j,rhod1(i,j),u1(i,j),v1(i,j)
         enddo
          write(51,'(bn)')
	enddo

          write(51,'(bn)')
          write(51,'(bn)')

 99     format(2I6,3(1x,e13.6))

        
        return
        end

c ============================
        subroutine energy
c ============================
        implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------------------------

! Energy

        g_eff=gnn+gnnn

        en_pot=0.d0
        en_cin=0.d0

        do j = 1, ny
         do i = 1, nx
            en_pot=en_pot+g_eff*psi(i,j)*psi(i,j)*cs2	!/2.d0
            en_cin=en_cin+rhod1(i,j)*cs2
         enddo
        enddo

        en_pot_1=0.d0
        en_pot_2=0.d0

!        write(6,*)'i valori dei c_2 per la subroutine energy sono'
!        write(6,*) c1_2,c2_2,c4_2,c5_2,c8_2


         do j = 1, ny
          do i = 1, nx

! computing forces to find total potential energy

           en_pot_1=en_pot_1+0.5*gnn*psi(i,j)*
     & (w(1)*c1_2*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+
     & w(5)*c2_2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2=en_pot_2+0.5d0*gnnn*psi(i,j)*
     &  (w1*c1_2*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+
     &   w2*c2_2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1))+
     &   w4*c4_2*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+
     &   w5*c5_2*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+
     &  w8*c8_2*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))


         enddo
        enddo

         write(86,98) istep,en_pot_1,en_pot_2,en_pot_1+en_pot_2,en_cin,
     &                en_pot,(en_pot_1+en_pot_2)/en_cin

98      format(1I6,1x,6(1x,e13.6)) 

! potential energy without "C_i^2" (Prof. Succi 11 - 03 - 2007)

        en_pot_1_bis=0.d0
        en_pot_2_bis=0.d0

         do j = 1, ny
          do i = 1, nx


           en_pot_1_bis=en_pot_1_bis+0.5*gnn*psi(i,j)*
     &  (4./9.)*psi(i,j)+
     & (w(1)*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+
     &  w(5)*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2_bis=en_pot_2_bis+0.5*gnnn*psi(i,j)*
     &  (w0*psi(i,j)+
     &   w1*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+
     &   w2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1))+
     &   w4*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+
     &   w5*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+
     &   w8*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))


         enddo
        enddo

         write(85,96) istep,en_pot_1_bis,en_pot_2_bis,
     &                en_pot_1_bis+en_pot_2_bis,en_pot

96      format(1I6,1x,4(1x,e13.6))


!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
! Diagnostica a meno di G1 e G2
!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]


        en_pot_ter=0.d0

!        do j = 1, ny
!         do i = 1, nx
!            en_pot_ter=en_pot_ter+psi(i,j)*psi(i,j)*cs2   !/2.d0
!         enddo
!        enddo



        en_pot_1_ter=0.d0
        en_pot_2_ter=0.d0

        en_pot_1_bulk=0.d0
        en_pot_2_bulk=0.d0


!         write(*,*)'i w valgono,',w(0),w(1),w(5)


         do j = 1, ny
          do i = 1, nx

           en_pot_1_ter=en_pot_1_ter+0.5d0*psi(i,j)*
     & ((4.d0/9.d0)*psi(i,j)+
     & w(1)*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+
     & w(5)*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2_ter=en_pot_2_ter+0.5d0*psi(i,j)*
     &  (w0*psi(i,j)+
     &   w1*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+
     &   w2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1))+
     &   w4*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+
     &   w5*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+
     &   w8*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))



           en_pot_1_bulk=en_pot_1_bulk+0.5d0*psi(i,j)*
     & ((4.d0/9.d0)*psi(i,j)+
     &  w(1)*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+
     &  w(5)*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2_bulk=en_pot_2_bulk+0.5d0*psi(i,j)*
     &  (w0*psi(i,j)+
     &   w1*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+
     &   w2*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+
     &   w4*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+
     &   w5*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)+
     &            psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+
     &   w8*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))


         enddo
        enddo

97       format(1I6,1x,5(1x,e13.6))

         en_2_1=0.d0
         en_2_2=0.d0
         en_2_1_bulk=0.d0
         en_2_2_bulk=0.d0

         do j = 1, ny
          do i = 1, nx

          en_2_1=en_2_1+0.5d0*psi(i,j)*
     &  (w0*psi(i,j)+
     &   w1*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+
     &   w2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))


          en_2_2=en_2_2+0.5d0*psi(i,j)*
     &   (w4*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+
     &   w5*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+
     &   w8*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))



          en_2_1_bulk=en_2_1_bulk+0.5d0*psi(i,j)*
     &  (w0*psi(i,j)+
     &   w1*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+
     &   w2*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))

          en_2_2_bulk=en_2_2_bulk+0.5d0*psi(i,j)*
     &  (w4*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+
     &   w5*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)+
     &            psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+
     &   w8*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))


         enddo
        enddo

         write(81,91) istep,en_pot_1_ter,en_pot_2_bulk,
     &                en_2_1,en_2_2,
     &                en_2_1_bulk,en_2_2_bulk,en_cin

         write(89,*) istep,gnn*en_pot_1_ter+gnnn*(en_2_1+en_2_2)


91         format(1I6,1x,7(1x,e13.6))

        

!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

        nliq=0

        do i=1,nx
          do j=1,ny

            if(rhod1(i,j).gt.(0.7d0)) then
              nliq=nliq+1
            endif

          enddo
        enddo 


!  number of bubbles...

 
        write(84,*)istep,nliq

        return
        end
c ============================
	subroutine out0d
c ============================
	implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------------------------

	densit = 0.0d0

	do k = 0, npop-1
	   do j= 1, ny
	      do i = 1, nx
	         densit = densit + f(k,i,j)
	      enddo
	   enddo
	enddo

	densit = densit / dfloat(nx*ny) 

	do i = 1, nx
	   do j= 1, ny
              param(i,j)=0.d0
	      do k = 0, npop-1
	         param(i,j) = param(i,j) + f(k,i,j)
	      enddo
	   enddo
	enddo
	do i = 1, nx
	   do j= 1, ny
	      param(i,j) = (param(i,j)-rhoaver)**2.d0
	   enddo
	enddo
        paramord=0.d0
	do i = 1, nx
	   do j= 1, ny
	      paramord = paramord+param(i,j)
	   enddo
	enddo
        paramord = paramord / dfloat(nx*ny)
        paramord = dsqrt(paramord) / rhoaver

	umoy = 0.0d0
	vmoy = 0.0d0

	do j = 1, ny
	   do i = 1, nx
	      umoy = umoy + u1(i,j)
	      vmoy = vmoy + v1(i,j)
	   enddo
	enddo
	
	umoy = umoy / dfloat(nx*ny)
	vmoy = vmoy / dfloat(nx*ny)

c	print*,'istep density umoy and vmoy'
        write(*,*)'Debug'  
        write(*,'(6e18.10)') dfloat(istep),densit,umoy,vmoy

        write(92,'(I6,3(1x,e13.6))')istep,densit,umoy,vmoy

c	print*,'paramord'
c        write(*,'(6e18.10)') paramord

        rhomax = -100000000.d0
        rhomin = 100000000.d0
        do i = 1, nx
          do j = 1, ny
            rhomax = max(rhomax,rhod1(i,j))
            rhomin = min(rhomin,rhod1(i,j)) 
        
         enddo
        enddo

        densit = densit / dfloat(nx*ny)  ! bulk 

        print *,' rhomax rhomin'
        write(*,'(6e18.10)')densit,rhomax,rhomin

       print*,'============================================='  


!       p=0.d0
!
!       do j = 1,ny
!          do i = 1,nx
!            pid  =rhod1(i,j)*cs2
!            pnid =0.5d0*(gnn+gnnn)*cs2*psi(i,j)*psi(i,j)     !0.5d0*0.7365d0/3.d0*gnn*psi(i,j)*psi(i,j)+
!!                                                          &  0.2635d0/3.d0*gnnn*psi(i,j)*psi(i,j) !*psi(i,j)*psi(i,j)
!            p(i,j)=pid+pnid
!             if(p(i,j).gt.pmax)then
!                pmax = p(i,j)
!             endif
!             if(p(i,j).lt.pmin)then
!                pmin = p(i,j)
!             endif
!           enddo
!       enddo
!

	return
	end
c ============================
        subroutine savepop
c ============================
        implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------------------------

       write(6,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
       write(6,*)'savin populations  ....'
       write(6,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'  


       rewind(111)
       rewind(113)
       rewind(112)

        do j=0,ny+1
         do i=0,nx+1
          do k=0,npop-1 
            write(111)f(k,i,j)      
          enddo
         enddo
        enddo 
     
        do j=1,ny
         do i=1,nx
            write(113)u1(i,j),v1(i,j)
         enddo
        enddo

        do j=-1,ny+2
         do i=-1,nx+2
            write(112)rhod1(i,j)
         enddo
        enddo


!        write(114)f


 
        write(6,*)'.....done!!! :^)'
   
       return
       end
       
c ============================
        subroutine resume 
c ============================
        implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------------------------

       write(6,*)'################################'
       write(6,*)'LETTURA popolazioni in corso....'
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
       end

c ============================
        subroutine calcsurf
c ============================
        implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------------------------
   

         do j=1,ny
           do i=1,nx
              iflag(i,j)=0
              if(rhod1(i,j).gt.0.865d0)then 
                iflag(i,j)=1
              endif
           enddo
         enddo

         isurf=0

         do j=1,ny
           do i=1,nx

              if(iflag(i,j).ne.iflag(i,j+1).or.
     &           iflag(i,j).ne.iflag(i,j-1).or.
     &           iflag(i,j).ne.iflag(i+1,j).or.
     &           iflag(i,j).ne.iflag(i-1,j)) then
                 isurf=isurf+1
              endif

           enddo
         enddo    

         isurf=isurf/2.d0

         write(88,*)istep,isurf


          return
          end

c ============================
        subroutine fourier
c ============================
        implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------------------------


           integer L,i,j,k,k1,q,kx,ky,interf,in,jn
           parameter (L=nx)!lattice size,must be 2^n,the same in x & y directions
   
           real*8 dx !lattice spacing
!           parameter (dx=1.d0)
           real*8 pi
           real*8 dq !spacing in Fourier space
           real*8 C,av,xc,yc,RRRR,z,s1,radius
           parameter (pi=3.14159265d0)
           real*8 G(-L/2:L/2,-L/2:L/2) !structure factor of density
           real*8 fourdata(1:2*L**2) !FFT density
           real*8 x1(1:L,1:L) !density field
           real*8 x(1:L,1:L) !density field with zero average
           integer nn(1:2)
           real*8 SC(0:L),SCC(0:L)


           nn(1)=L
           nn(2)=L
           dq=2.d0*pi/(L*dx) !spacing in Fourier space
   
   
!              open(23,file='BGK01.ruv2d',status='old')
!                  do j=1,L
!                   do i=1,L
!                    read(23,*)kx,ky,x1(i,j)
!                   enddo
!                 enddo
!                close(23)
!   
!             do j=1,L
!               do i=1,L
!                 x1(i,j)=dsin(2.d0*pi/L*8*(i))
!                 write(49,*)x1(i,j)
!               enddo
!             enddo


             do j=1,L
               do i=1,L

                  x1(i,j)=rhod1(i,j)

               enddo
             enddo
   
                av=0.d0
                do i=1,L
                  do j=1,L
   
                    av=av+x1(i,j)
   
                  enddo
                enddo
                av=av/L/L
   
                    print*,av
   
                do i=1,L
                  do j=1,L
   
                    x(i,j)=x1(i,j)-av
   
                  enddo
                enddo
   
   
                   do i=1,L
                      do j=1,L
   
                   fourdata(2*((i-1)*L+j-1)+1) = x(i,j)
                   fourdata(2*((i-1)*L+j-1)+2) = 0.d0
   
                      enddo
                   enddo
   
   
                   call fourn(fourdata,nn,2,1)
   
           do i=1,L
           do j=1,L
   
           C= ((fourdata(2*((i-1)*L+j-1)+1))**2 
     &       +(fourdata(2*((i-1)*L+j-1)+2))**2)!**(0.5d0)
           C=C/(L*L)
   
           if(i.le.L/2.and.j.le.L/2)then
           kx=i-1
           ky=j-1
           G(kx,ky)=C
           elseif(i.gt.L/2.and.j.le.L/2)then
           kx=i-L-1
           ky=j-1
           G(kx,ky)=C
           elseif(i.le.L/2.and.j.gt.L/2)then
           kx=i-1
           ky=j-L-1
           G(kx,ky)=C
           elseif(i.gt.L/2.and.j.gt.L/2)then
           kx=i-L-1
           ky=j-L-1
           G(kx,ky)=C
           endif
           enddo
           enddo
   
!           do ky=-L/2,L/2-1
!           do kx=-L/2,L/2-1
!             write(50,*)kx*dq,ky*dq,G(kx,ky)
!           enddo
!           enddo
   
           do k=0,L
           SC(k)=0.d0
           SCC(k)=1.d0
           enddo
           do kx=-L/2,L/2-1
           do ky=-L/2,L/2-1
           xc=float(kx)+0.5
           yc=float(ky)+0.5
           RRRR=(xc*xc+yc*yc)**(0.5)
           do k=1,L
           if(RRRR.lt.(float(k)+0.5).and.RRRR.gt.(float(k)-0.5))then
           SC(k)=SC(k)+G(kx,ky)  !(kx*kx+ky*ky)**(0.5d0)*G(kx,ky)
           SCC(k)=SCC(k)+1.
           endif
           enddo
           enddo
           enddo
           Z=0.d0
           S1=0.d0
           do k=1,L/2
           SC(k)=SC(k)/SCC(k)
   
            write(177,*)k,SC(k)
   
           Z=Z+dq*SC(k)
           q=k
           S1=S1+dq*abs(dq*q)*SC(k)
           enddo
            write(177,'(bn)')
            write(177,'(bn)')  

           radius=pi*Z/S1
   
                  interf=0
                  do i=1,L
                  do j=1,L
                  in=i+1
                  if(in.gt.L)in=1
                  if(x(i,j)*x(in,j).le.0.d0)then
                  interf=interf+1
                  endif
                  jn=j+1
                  if(jn.gt.L)jn=1
                  if(x(i,j)*x(i,jn).le.0.d0)then
                  interf=interf+1
                  endif
                  enddo
                  enddo
   
   
!                  print*,'radius =',radius,'1/interf =',1.d0/interf
   
                  write(117,*)istep,radius
!     &           (float(L)/2)*radius/6.283185307179586476925287d0
   
!c radius e 1/interf sono due misure diverse della taglia dei domini

!c radius e' ottenuta come l'inverso del primo momento del fattore di
!c    struttura mediato circolarmente
   
!c 1/interf e' l'inverso della lunghezza delle interfacce del sistema
   
!         do kx=-L/2,L/2
!            do ky=-L/2,L/2
!               write(130,*)kx,ky,G(kx,ky)
!               if(G(kx,ky).gt.100)then
!                  write(*,*) 'coordinatei picco', kx, ky
!               endif
!            enddo
!         enddo
   
           return
           end
   
   
   
   
           subroutine fourn(data,nn,ndim,isign)
           implicit real*8 (a-h), real*8 (o-z)
           dimension nn(ndim),data(*)
           ntot=1
           do idim=1,ndim
              ntot=ntot*nn(idim)
           enddo
           nprev=1
           do idim=1,ndim
              n=nn(ndim)
              nrem=ntot/(n*nprev)
              ip1=2*nprev
              ip2=ip1*n
              ip3=ip2*nrem
              i2rev=1
              do i2=1,ip2,ip1
                 if(i2.lt.i2rev) then
                    do i1=i2,i2+ip1-2,2
                       do i3=i1,ip3,ip2
                          i3rev=i2rev+i3-i2
                          tempr=data(i3)
                          tempi=data(i3+1)
                          data(i3)=data(i3rev)
                          data(i3+1)=data(i3rev+1)
                          data(i3rev)=tempr
                          data(i3rev+1)=tempi
                       enddo
                    enddo
                 endif
                 ibit=ip2/2
1              if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
                    i2rev=i2rev-ibit
                    ibit=ibit/2
                    goto 1
                 endif
                 i2rev=i2rev+ibit
              enddo
              ifp1=ip1
2            if(ifp1.lt.ip2) then
                 ifp2=2*ifp1
                 theta=isign*6.283185307179586476925287d0/(ifp2/ip1)
                 wpr=-2d0*sin(5d-1*theta)**2
                 wpi=sin(theta)
                 wr=1d0
                 wi=0d0
                 do i3=1,ifp1,ip1
                    do i1=i3,i3+ip1-2,2
                       do i2=i1,ip3,ifp2
                          k1=i2
                          k2=k1+ifp1
                          tempr=wr*data(k2)-wi*data(k2+1)
                          tempi=wr*data(k2+1)+wi*data(k2)
                          data(k2)=data(k1)-tempr
                          data(k2+1)=data(k1+1)-tempi
                          data(k1)=data(k1)+tempr
                          data(k1+1)=data(k1+1)+tempi
                       enddo
                    enddo
                    wtemp=wr
                    wr=wr*wpr-wi*wpi+wr
                    wi=wi*wpr+wtemp*wpi+wi
                 enddo
                ifp1=ifp2
                 goto 2
              endif
              nprev=n*nprev
           enddo
           return
           end
   
  

c ============================
        subroutine laplace
c ============================
        implicit double precision(a-h,o-z)
        include 'muphase.par'
c----------------------------------------------------------

        p=0.d0

        do j = 1,ny
          do i = 1,nx
            pid  =rhod1(i,j)*cs2
            pnid =0.5d0*(gnn+gnnn)*cs2*psi(i,j)*psi(i,j)
!0.5d0*0.7365d0/3.d0*gnn*psi(i,j)*psi(i,j)+
!                                                          &
!                                                          0.2635d0/3.d0*gnnn*psi(i,j)*psi(i,j)
!                                                          !*psi(i,j)*psi(i,j)
            p(i,j)=pid+pnid
             if(p(i,j).gt.pmax)then
                pmax = p(i,j)
             endif
             if(p(i,j).lt.pmin)then
                pmin = p(i,j)
             endif

           write(65,*)i,j,p(i,j)

           enddo

           write(65,'(bn)')

         enddo

         do i=1,nx
            write(66,*)i,p(i,ny/2)
         enddo


            rhoaver=0.d0

        rhomax = -100000000.d0
        rhomin = 100000000.d0
        do i = 1, nx
          do j = 1, ny
            rhomax = max(rhomax,rhod1(i,j))
            rhomin = min(rhomin,rhod1(i,j))

         enddo
        enddo

         write(*,*) 'rho_min=',rhomin,'rho_max=',rhomax
         write(*,*) '......................................'

         write(*,*) '..computing mean rho...'

         do j=1,ny
          do i=1,nx
              rhoaver = rhoaver + rhod1(i,j)
          enddo
         enddo

         write(*,*)'A_1 =',gnn+gnnn
         write(*,*)'A_2 =',gnn+1.5d0*gnnn


              rhoaver = rhoaver/(nx*ny)

         write(*,*) 'rhoaver =', rhoaver

         radius = ((rhoaver*(nx*ny)-rhomin*(nx*ny))/(3.1415926536d0*
     &             (rhomax-rhomin)))**(0.5d0)

         write(*,*) 'Bubble radius (Case 1)', radius

          surf_tens = (p(nx/2,ny/2)-p(nx/16,ny/16))*radius

          write(*,*) 'Laplace test: gamma = ', surf_tens

          return
          end
 
