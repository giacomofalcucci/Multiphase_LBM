!
        module storage
!
        implicit double precision(a-h,o-z)
!
!       parameter (nx =  128, ny =  128, npop = 9)
!       parameter (nx =  192, ny =  192, npop = 9)
        parameter (nx =  256, ny =  256, npop = 9)
!       parameter (nx =  512, ny =  512, npop = 9)
!       parameter (nx = 1024, ny = 1024, npop = 9) 
!
        character*5 fileout

! /precision/
        integer, parameter :: sp = kind(1.0)
        integer,parameter:: dp = selected_real_kind(2*precision(1.0_sp))
!
! /constants/
        real(dp):: cs2,cs22,cssq,rhoin,omega,fpois,den,visc
        real(dp)::  w0,w1,w2,w4,w5,w8,gnn,gnnn,rhoaver,dinvrho
        real(dp):: rhopsi,dt,dx,dump,c1_2,c2_2,c4_2,c5_2,c8_2
!                
! /phys/   u0,uf,fom
        real(dp):: u0,uf,fom
!                
        common /arrays/ u1(0:nx+1,0:ny+1), &
                        v1(0:nx+1,0:ny+1), &
                        u2(0:nx+1,0:ny+1), & 
                        v2(0:nx+1,0:ny+1), &
                        psi(-1:nx+2,-1:ny+2), &
                        rhod1(-1:nx+2,-1:ny+2), &
                        rhod2(-1:nx+2,-1:ny+2), &
                        p(0:nx+1,0:ny+1), &
                        iflag(0:nx+1,0:ny+1),  &
                        param(1:nx,1:ny), &
                        feq(0:npop-1,0:nx+1,0:ny+1), &
                        f(0:npop-1,0:nx+1,0:ny+1), &
                        fp(0:npop-1,0:nx+1,0:ny+1), &
                        f_guo, &
                        w(0:npop-1), & 
                        u_ci(0:npop-1)
!                
! /count/ 
        integer:: istep,nout,ndiag,nsteps,nobst
!        
        common /ile/ fileout
!        
! /noutput/i
        integer:: nrhout
!        
! /logic/ 
        logical:: iforce,iobst
!        
        end module  storage

