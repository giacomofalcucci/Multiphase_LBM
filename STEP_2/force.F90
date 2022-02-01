!----------------------------------------------------------
        subroutine force
!--------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        frce = fpois   ! constant external force  
        do i=-1,nx+2
           do j=-1,ny+2    ! ho levato le BC per psi!!!
              psi(i,j) =rhopsi * (1.d0 -exp(- rhod1(i,j) / rhopsi))
           enddo
        enddo      

        wnn=w1/cs2
        wnnn=w2/cs2

        do j = 1, ny
           do i = 1, nx
! non local Lennard-Jones (may put a density functional, not density itself)

              fnnx  = psi(i,j)*(psi(i+1,j)-psi(i-1,j))
              fnny  = psi(i,j)*(psi(i,j+1)-psi(i,j-1))
              fnnnx = psi(i,j)*((psi(i+1,j-1) + psi(i+1,j+1)    &
                               - psi(i-1,j-1) - psi(i-1,j+1)))
              fnnny = psi(i,j)*((psi(i+1,j+1) + psi(i-1,j+1)    &
                               - psi(i-1,j-1) - psi(i+1,j-1)))

! interactions
              f2x =-(gnnn * fnnx * w1 + gnnn * fnnnx * w2)- &  ! prima era solo -(gnn * fnnx * w1 + gnn * fnnnx * w2)
               (gnn * fnnx*1./9.+ gnn * fnnnx *1./36.)  ! aggiunte per completare la simmetria del reticolo (Chib)
              f2y =-(gnnn * fnny * w1 + gnnn * fnnny * w2)-  & ! prima era solo -(gnn * fnnx * w1 + gnn * fnnnx * w2)
               (gnn * fnny*1./9.+ gnn * fnnny *1./36.)  ! aggiunte per completare la simmetria del reticolo (Chib) 

! SHIFT Equilibrium
              u1(i,j)=u1(i,j)+f2x/(omega*rhod1(i,j))
              v1(i,j)=v1(i,j)+f2y/(omega*rhod1(i,j))
           end do
        end do

        return
        end subroutine force
