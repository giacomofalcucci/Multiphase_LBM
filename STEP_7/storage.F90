!
        module storage
!
        implicit double precision(a-h,o-z)
!
        integer, parameter:: sp=kind(1.0)
        integer, parameter:: dp=selected_real_kind(2*precision(1.0_sp))
!
        integer :: nx
        integer :: ny
        integer, parameter :: npop = 9
!
        character*5 fileout
        integer::  nrhout
        logical iforce,iobst
!/phys/ 
        real(dp) u0,v0,uf,fom  
        real(dp) f_guo
!/constants/
        real(dp) cs2,cs22,cssq,rhoin,omega,fpois,den,visc
        real(dp) w0,w1,w2,w4,w5,w8,gnn,gnnn,rhoaver,dinvrho
        real(dp) rhopsi,dt,dx,dump,c1_2,c2_2,c4_2,c5_2,c8_2 
        real(dp), parameter :: cte04 = (4.d0/ 9.d0)
        real(dp), parameter :: cte09 = (1.d0/ 9.d0)
        real(dp), parameter :: cte36 = (1.d0/36.d0)
!/count/ 
        integer istep,nout,ndiag,nsteps,nobst,icond
!/arrays to allocate/
        real(dp), dimension (:,:), allocatable :: f0
        real(dp), dimension (:,:), allocatable :: f1
        real(dp), dimension (:,:), allocatable :: f2
        real(dp), dimension (:,:), allocatable :: f3
        real(dp), dimension (:,:), allocatable :: f4
        real(dp), dimension (:,:), allocatable :: f5
        real(dp), dimension (:,:), allocatable :: f6
        real(dp), dimension (:,:), allocatable :: f7
        real(dp), dimension (:,:), allocatable :: f8
!
        real(dp), dimension (:,:), allocatable :: fp0
        real(dp), dimension (:,:), allocatable :: fp1
        real(dp), dimension (:,:), allocatable :: fp2
        real(dp), dimension (:,:), allocatable :: fp3
        real(dp), dimension (:,:), allocatable :: fp4
        real(dp), dimension (:,:), allocatable :: fp5
        real(dp), dimension (:,:), allocatable :: fp6
        real(dp), dimension (:,:), allocatable :: fp7
        real(dp), dimension (:,:), allocatable :: fp8
!
        real(dp), dimension (:,:), allocatable :: u1
        real(dp), dimension (:,:), allocatable :: v1
        real(dp), dimension (:,:), allocatable :: u2
        real(dp), dimension (:,:), allocatable :: v2
        real(dp), dimension (:,:), allocatable :: psi
        real(dp), dimension (:,:), allocatable :: rhod1
        real(dp), dimension (:,:), allocatable :: rhod2
        real(dp), dimension (:,:), allocatable :: p
        real(dp), dimension (:,:), allocatable :: param
!       
        integer, dimension (:,:), allocatable ::  iflag

!/arrays/
        real(dp), dimension (0:npop-1) ::     w
        real(dp), dimension (0:npop-1) ::     u_ci
!                
!        
        end module  storage

