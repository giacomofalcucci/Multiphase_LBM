c----------------------------------------------------------
        subroutine fourier

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
           end subroutine fourier
   
   
   
   
c----------------------------------------------------------
           subroutine fourn(data,nn,ndim,isign)

           implicit real*8 (a-h), real*8 (o-z)
           dimension nn(ndim),data(*)
c----------------------------------------------------------
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
           end subroutine fourn
   
  

