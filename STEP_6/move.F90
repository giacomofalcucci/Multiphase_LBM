! --------------------------------------------------
        subroutine move
!---------------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)

!$OMP PARALLEL DEFAULT(NONE)                            &
!$OMP PRIVATE(i,j)                                      &
!$OMP SHARED(nx,ny)                                     &
!$OMP SHARED(f1,f2,f3,f4,f5,f6,f7,f8)                   &
!$OMP SHARED(fp1,fp2,fp3,fp4,fp5,fp6,fp7,fp8)           
!$OMP DO
        do j = 1,ny
           do i = 1, nx
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
!$OMP END PARALLEL
        

! border fix

        j = 0
        do i = 0, nx+1 
           fp1(i,j)=f1(i,j)
           fp2(i,j)=f2(i,j)
           fp3(i,j)=f3(i,j)
           fp4(i,j)=f4(i,j)
           fp5(i,j)=f5(i,j)
           fp6(i,j)=f6(i,j)
           fp7(i,j)=f7(i,j)
           fp8(i,j)=f8(i,j)
        enddo
!
        j = ny+1
        do i = 0, nx+1 
           fp1(i,j)=f1(i,j)
           fp2(i,j)=f2(i,j)
           fp3(i,j)=f3(i,j)
           fp4(i,j)=f4(i,j)
           fp5(i,j)=f5(i,j)
           fp6(i,j)=f6(i,j)
           fp7(i,j)=f7(i,j)
           fp8(i,j)=f8(i,j)
        enddo
!
        i = 0
        do j = 0, ny+1
           fp1(i,j)=f1(i,j)
           fp2(i,j)=f2(i,j)
           fp3(i,j)=f3(i,j)
           fp4(i,j)=f4(i,j)
           fp5(i,j)=f5(i,j)
           fp6(i,j)=f6(i,j)
           fp7(i,j)=f7(i,j)
           fp8(i,j)=f8(i,j)
        enddo
!
        i = nx+1
        do j = 0, ny+1
           fp1(i,j)=f1(i,j)
           fp2(i,j)=f2(i,j)
           fp3(i,j)=f3(i,j)
           fp4(i,j)=f4(i,j)
           fp5(i,j)=f5(i,j)
           fp6(i,j)=f6(i,j)
           fp7(i,j)=f7(i,j)
           fp8(i,j)=f8(i,j)
        enddo
!
#ifdef DEBUG
        write(6,*) "Completed subroutine move"
#endif
        end subroutine move
