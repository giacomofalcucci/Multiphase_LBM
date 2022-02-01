!--------------------------------------------------
        subroutine alloca

        use storage
        implicit none
!
! -------vector allocation
!
        allocate( fp(0:npop-1,0:nx+1,0:ny+1))
        allocate(  f(0:npop-1,0:nx+1,0:ny+1))

        allocate(u1(0:nx+2,0:ny+2))
        allocate(v1(0:nx+2,0:ny+2))
        allocate(u2(0:nx+2,0:ny+2))
        allocate(v2(0:nx+2,0:ny+2))
        allocate(psi(-1:nx+2,-1:ny+2))
        allocate(rhod1(-1:nx+2,-1:ny+2))
        allocate(rhod2(-1:nx+2,-1:ny+2))
        allocate(p(0:nx+1,0:ny+1))
        allocate(iflag(0:nx+1,0:ny+1))
        allocate(param(1:nx,1:ny))

#ifdef DEBUG
        write(6,*) "Completed subroutine alloca"
#endif

        return
        end subroutine alloca

