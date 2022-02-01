c----------------------------------------------------------
        subroutine savepop

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

        write(6,*)'.....done!!! :^)'
   
        return
        end subroutine savepop
       
