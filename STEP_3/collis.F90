!----------------------------------------------------------
        subroutine collis
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        do j = 1, ny
           do i = 1, nx

! forcing

! non local Lennard-Jones (may put a density functional, not density itself)
              fnnx  = psi(i,j)*(psi(i+1,j)-psi(i-1,j))
              fnny  = psi(i,j)*(psi(i,j+1)-psi(i,j-1))
              fnnnx = psi(i,j)*((psi(i+1,j-1) + psi(i+1,j+1)    &
                               - psi(i-1,j-1) - psi(i-1,j+1)))
              fnnny = psi(i,j)*((psi(i+1,j+1) + psi(i-1,j+1)    &
                               - psi(i-1,j-1) - psi(i+1,j-1)))

! interactions
              f2x =-(gnnn * fnnx * w1 + gnnn * fnnnx * w2)- &  
                    (gnn * fnnx*1./9.+ gnn * fnnnx *1./36.)  
              f2y =-(gnnn * fnny * w1 + gnnn * fnnny * w2)- & 
                    (gnn * fnny*1./9.+ gnn * fnnny *1./36.)

! SHIFT Equilibrium
              u1(i,j)=u1(i,j)+f2x/(omega*rhod1(i,j))
              v1(i,j)=v1(i,j)+f2y/(omega*rhod1(i,j))
           enddo
        enddo

        do j = 1, ny
           do i = 1, nx
! compute equilibrium              
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

              feq0 = (4.d0/9.d0)*rhoij*(1.0d0 - sumsq)

              feq1 = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + u22 + ui)
              feq2 = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + v22 + vi)
              feq3 = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + u22 - ui)
              feq4 = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + v22 - vi)

              feq5 = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 +ui+vi+uv)
              feq6 = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 -ui+vi-uv)
              feq7 = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 -ui-vi+uv)
              feq8 = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 +ui-vi-uv)

! compute correction              
              f(0,i,j) = fp(0,i,j) * (1.0d0 - omega) + omega * feq0
              f(1,i,j) = fp(1,i,j) * (1.0d0 - omega) + omega * feq1
              f(2,i,j) = fp(2,i,j) * (1.0d0 - omega) + omega * feq2
              f(3,i,j) = fp(3,i,j) * (1.0d0 - omega) + omega * feq3
              f(4,i,j) = fp(4,i,j) * (1.0d0 - omega) + omega * feq4
              f(5,i,j) = fp(5,i,j) * (1.0d0 - omega) + omega * feq5
              f(6,i,j) = fp(6,i,j) * (1.0d0 - omega) + omega * feq6
              f(7,i,j) = fp(7,i,j) * (1.0d0 - omega) + omega * feq7
              f(8,i,j) = fp(8,i,j) * (1.0d0 - omega) + omega * feq8

           enddo
        enddo

#ifdef DEBUG
        write(6,*) "Completed subroutine collis"
#endif
        end subroutine collis

