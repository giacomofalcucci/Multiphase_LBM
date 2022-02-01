! --------------------------------------------------
        subroutine move
!---------------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)

!$acc kernels
!$acc loop independent collapse(2)
        do j = 1,ny
           do i = 1, nx
!              fp0(i,j) = f0(i,j)
              fp1(i,j) = f1(i-1,j)
              fp2(i,j) = f2(i,j-1)
              fp3(i,j) = f3(i+1,j)
              fp4(i,j) = f4(i,j+1)
              fp5(i,j) = f5(i-1,j-1)
              fp6(i,j) = f6(i+1,j-1)
              fp7(i,j) = f7(i+1,j+1)
              fp8(i,j) = f8(i-1,j+1)
           enddo
        enddo
!$acc end kernels

! border fix

        j = 0
!$acc kernels
!$acc loop independent 
        do i = 0, nx+1 
!           fp0(i,j)=f0(i,j)
           fp1(i,j)=f1(i,j)
           fp2(i,j)=f2(i,j)
           fp3(i,j)=f3(i,j)
           fp4(i,j)=f4(i,j)
           fp5(i,j)=f5(i,j)
           fp6(i,j)=f6(i,j)
           fp7(i,j)=f7(i,j)
           fp8(i,j)=f8(i,j)
        enddo
!$acc end kernels
!
        j = ny+1
!$acc kernels
!$acc loop independent 
        do i = 0, nx+1 
!           fp0(i,j)=f0(i,j)
           fp1(i,j)=f1(i,j)
           fp2(i,j)=f2(i,j)
           fp3(i,j)=f3(i,j)
           fp4(i,j)=f4(i,j)
           fp5(i,j)=f5(i,j)
           fp6(i,j)=f6(i,j)
           fp7(i,j)=f7(i,j)
           fp8(i,j)=f8(i,j)
        enddo
!$acc end kernels
!
        i = 0
!$acc kernels
!$acc loop independent 
        do j = 0, ny+1
!           fp0(i,j)=f0(i,j)
           fp1(i,j)=f1(i,j)
           fp2(i,j)=f2(i,j)
           fp3(i,j)=f3(i,j)
           fp4(i,j)=f4(i,j)
           fp5(i,j)=f5(i,j)
           fp6(i,j)=f6(i,j)
           fp7(i,j)=f7(i,j)
           fp8(i,j)=f8(i,j)
        enddo
!$acc end kernels
!
        i = nx+1
!$acc kernels
!$acc loop independent 
        do j = 0, ny+1
!           fp0(i,j)=f0(i,j)
           fp1(i,j)=f1(i,j)
           fp2(i,j)=f2(i,j)
           fp3(i,j)=f3(i,j)
           fp4(i,j)=f4(i,j)
           fp5(i,j)=f5(i,j)
           fp6(i,j)=f6(i,j)
           fp7(i,j)=f7(i,j)
           fp8(i,j)=f8(i,j)
        enddo
!$acc end kernels
!
#ifdef DEBUG
        write(6,*) "Completed subroutine move"
#endif
        end subroutine move
