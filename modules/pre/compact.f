!=======================================================================
!       Package of Subroutines for Compact Finite Difference  !
!=======================================================================
!**********************************************************************
! pred1m60:
!
! generic coefficients alpha,beta,delta y
!      a0,a1,a2 y a3 for the plan of 
!      compact finite difference with uniform mesh.
!            first derivative
!
! delta*(du(j+2)+du(j-2))beta*(du(j+2)+du(j-2))+
!      alpha*(du(j+1)+du(j-1))+du(j) =
!    = a1*(u(j+1)-u(j-1)) + a2(u(j+2)-u(j-2)) + a3(u(j+2)-u(j-3))
!
!    on the edges of the mesh:
!    alpha*du(2) + du(1) =
!     = a1*u(1) + a2*u(2) + a3*u(3) + a4*u(4) + a5*u(5)
!
!    alpha*(du(3)+du(1)) + du(2) =
!     = a1*(u(3)-u(1))
!
!*********************************************************************
      subroutine pred1m6o(ny,den,ben,aln,one,alp,bep,dep,a3n,
     &                    a2n,a1n,a0,a1p,a2p,a3p)
      implicit none
      integer ny
      real(8) alpha,beta,delta,aaa,bbb,ccc
      real(8) den(ny),ben(ny),aln(ny),one(ny),alp(ny),bep(ny),dep(ny)
      real(8) a3n(ny),a2n(ny),a1n(ny), a0(ny),a1p(ny),a2p(ny),a3p(ny)

      real(8) h,pi
      integer j

      pi = 4d0*atan(1d0)

      do j=1,ny
         one(j)=1d0
      enddo

      h=2d0/dble(ny-1)


!-----------------------------------------------------------------------
!     center of the control
!     impose the coefficients from external calculations
! ----------------------------------------------------------------------

      aaa   =  0.56926818603810d0                   !=    21d0/32d0  !
      bbb   =  0.31997939759919d0                   !=   231d0/1d3   !
      ccc   =  0.03751544855390d0                   !=    49d0/4d3   !
      alpha = -3d0/4d0  +(33d0/20d0  )*ccc          !=    9d0/16d0   !
     .        +(239d0/180d0)*aaa + (82d0/45d0 )*bbb !                !
      beta  =  3d0/10d0 +(24d0/25d0  )*ccc          !=    9d0/1d2    !
     .        -(88d0/225d0 )*aaa + (34d0/225d0)*bbb !                !
      delta = -1d0/2d1  +(39d0/100d0 )*ccc          !=     1d0/4d2   !
     .      +(19d0/300d0 )*aaa + (2d0/75d0  )*bbb   !                !

      do j=4,ny-3
         den(j)=delta
         ben(j)=beta
         aln(j)=alpha
         alp(j)=alpha
         bep(j)=beta
         dep(j)=delta
         a3n(j)=-ccc/h
         a2n(j)=-bbb/h
         a1n(j)=-aaa/h
         a0(j)=0d0
         a1p(j)= aaa/h
         a2p(j)= bbb/h
         a3p(j)= ccc/h
      enddo

!-----------------------------------------------------------------------
!
!       
!       using the information available
!       impose equations of different order
!
!-----------------------------------------------------------------------

      j=1
!aaf 2 points stencil used |ox-----...|
      den(j) =  0d0  !0d0          !0d0         !
      ben(j) =  0d0!0d0          !0d0         !
      aln(j) =  0d0!0d0          !0d0         !
      alp(j) =  1d0!9d0          !6d0         !
      bep(j) =  0d0!alp(j)       !3d0         !
      dep(j) =  0d0!1.d0         !0d0         !

      a3n(j) =  0d0 !0d0          !0d0         !
      a2n(j) =  0d0!0d0          !0d0         !
      a1n(j) =  0d0!0d0          !0d0         !
      a0 (j) =  -2d0/h !-11d0/3d0/h  !-10d0/3d0/h !
      a1p(j) =  2d0/h !-9d0/h       !-3d0/h      !
      a2p(j) =  0d0  !-a1p(j)      !6d0/h       !
      a3p(j) =  0d0  !-a0(j)       !1d0/3d0/h   !

      j=2
!aaf 3 points center stencil |xox-----...|
      den(j) =  0d0  !0d0
      ben(j) =  0d0  !0d0
      aln(j) =  1d0/4d0  !1d0/16d0
      alp(j) =  1d0/4d0  !9d0/4d0
      bep(j) =  0d0  !1d0
      dep(j) =  0d0  !aln(j)

      a3n(j) =  0d0  !0d0
      a2n(j) =  0d0  !0d0
      a1n(j) =  -3d0/4d0/h  !-25d0/96d0/h
      a0 (j) =  0d0  !-5d0/3d0/h
      a1p(j) =  3d0/4d0/h  !0d0
      a2p(j) =  0d0   !-a0(j)
      a3p(j) =  0d0   !-a1n(j)

!aaf center sentcil with 5 points |xxoxx----...|
      j=3
      den(j) =  0d0 !0d0
      ben(j) =  1d0/36d0 !1d-2          ! 
      aln(j) =  4d0/9d0 !1d0/4d0       ! 
      alp(j) =  4d0/9d0 !1d0           ! 
      bep(j) =  1d0/36d0 !aln(j)
      dep(j) =  0d0 !ben(j)

      a3n(j) =  0d0 !0d0
      a2n(j) =  -25d0/216d0/h !-137d0/3d3/h  !
      a1n(j) =  -20d0/27d0/h !-13d0/24d0/h  !
      a0 (j) =  0d0 !-2d0/3d0/h    !
      a1p(j) =  20d0/27d0/h !-a0(j)
      a2p(j) =  25d0/216d0/h !-a1n(j)
      a3p(j) =  0d0 !-a2n(j)

! -----------------  the other end --------------

      j=ny
      den(j) =   dep(1)
      ben(j) =   bep(1)
      aln(j) =   alp(1)
      alp(j) =   aln(1)
      bep(j) =   ben(1)
      dep(j) =   den(1)
      a3n(j) =   -a3p(1)
      a2n(j) =   -a2p(1)
      a1n(j) =   -a1p(1)
      a0 (j) =   -a0(1)
      a1p(j) =   -a1n(1)
      a2p(j) =   -a2n(1)
      a3p(j) =   -a3n(1)

      j=ny-1
      den(j) =   dep(2)
      ben(j) =   bep(2)
      aln(j) =   alp(2)
      alp(j) =   aln(2)
      bep(j) =   ben(2)
      dep(j) =   den(2)
      a3n(j) =   -a3p(2)
      a2n(j) =   -a2p(2)
      a1n(j) =   -a1p(2)
      a0 (j) =   -a0(2)
      a1p(j) =   -a1n(2)
      a2p(j) =   -a2n(2)
      a3p(j) =   -a3n(2)
  
      j=ny-2
      den(j) =    dep(3)
      ben(j) =    bep(3)
      aln(j) =    alp(3)
      alp(j) =    aln(3)
      bep(j) =    ben(3)
      dep(j) =    den(3)
      a3n(j) =    -a3p(3)
      a2n(j) =    -a2p(3)
      a1n(j) =    -a1p(3)
      a0 (j) =    -a0(3)
      a1p(j) =    -a1n(3)
      a2p(j) =    -a2n(3)
      a3p(j) =    -a3n(3)

      return
      end subroutine pred1m6o

!***********************************************************************
! prederiv2:
!
! generate the coefficients alpha,beta,a0,a1 y a2 for the compact
! finite difference plan. second derivative
!
! beta*(d2u(j+2)+d2u(j-2))+alpha*(d2u(j+1)+d2u(j-1))+d2u(j) =
!  = a0*u(j)+a1*(u(j+1)+u(j-1)) + a2(u(j+2)+u(j-2))
!
!
!
!--------------------------------------------------------------------
      subroutine prederiv2(ny,ben,aln,aaa,alp,bep,a2n,a1n,a0,a1p,a2p,x)

      implicit none

      integer ny,n
      real(8) x(ny)
      real(8) ben(ny),aln(ny),aaa(ny),alp(ny),bep(ny)
      real(8) a2n(ny),a1n(ny),a0(ny),a1p(ny),a2p(ny)
      real(8) c(9,9),b(9)
   
      integer j

      do j=1,ny
      aaa(j)=1d0
      enddo

!		Ahora calculo los coeficientes del centro del dominio
  
      do j=3,ny-2
      c(1,1)=0d0
      c(1,2)=0d0
      c(1,3)=0d0
      c(1,4)=0d0
      c(1,5)=1d0
      c(1,6)=1d0
      c(1,7)=1d0
      c(1,8)=1d0
      c(1,9)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=0d0
      c(2,4)=0d0
      c(2,5)=(x(j-2)-x(j))
      c(2,6)=(x(j-1)-x(j))
      c(2,7)=0d0
      c(2,8)=(x(j+1)-x(j))
      c(2,9)=(x(j+2)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=2d0
      c(3,3)=2d0
      c(3,4)=2d0
      c(3,5)=-(x(j-2)-x(j))**2
      c(3,6)=-(x(j-1)-x(j))**2
      c(3,7)=0d0
      c(3,8)=-(x(j+1)-x(j))**2
      c(3,9)=-(x(j+2)-x(j))**2
      b(3)=-2d0

      do n=3,8
      c(n+1,1)=n*(n-1)*(x(j-2)-x(j))**(n-2)
      c(n+1,2)=n*(n-1)*(x(j-1)-x(j))**(n-2)
      c(n+1,3)=n*(n-1)*(x(j+1)-x(j))**(n-2)
      c(n+1,4)=n*(n-1)*(x(j+2)-x(j))**(n-2)
      c(n+1,5)=-(x(j-2)-x(j))**n
      c(n+1,6)=-(x(j-1)-x(j))**n
      c(n+1,7)=0d0
      c(n+1,8)=-(x(j+1)-x(j))**n
      c(n+1,9)=-(x(j+2)-x(j))**n
      b(n+1)=0d0
      enddo
      
!		call LSLRG(9,c,9,b,1,sol)
      call gaussj(c,9,9,b,1,1)
      

      ben(j)=b(1)
      aln(j)=b(2)
      alp(j)=b(3)
      bep(j)=b(4)
      a2n(j)=b(5)
      a1n(j)=b(6)
      a0(j)=b(7)
      a1p(j)=b(8)
      a2p(j)=b(9)


      enddo


!			fronteras del dominio  

!*************************************************************************
!*************************************************************************
      j=1
      c(1,1)=0d0
      c(1,2)=0d0
      c(1,3)=1d0
      c(1,4)=1d0
      c(1,5)=1d0
      b(1)=0d0
      
      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=0d0
      c(2,4)=(x(j+1)-x(j))
      c(2,5)=(x(j+2)-x(j))
!      c(2,5)=(x(j+3)-x(j))
!      c(2,6)=(x(j+4)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=2d0
      c(3,3)=0d0
      c(3,4)=-(x(j+1)-x(j))**2
      c(3,5)=-(x(j+2)-x(j))**2
!      c(3,5)=-(x(j+3)-x(j))**2
!      c(3,6)=-(x(j+4)-x(j))**2
      b(3)=-2d0

!      do n=3,5
      do n=3,4
      c(n+1,1)=n*(n-1)*(x(j+1)-x(j))**(n-2)
      c(n+1,2)=n*(n-1)*(x(j+2)-x(j))**(n-2)
      c(n+1,3)=0d0
      c(n+1,4)=-(x(j+1)-x(j))**n
      c(n+1,5)=-(x(j+2)-x(j))**n
!      c(n+1,5)=-(x(j+3)-x(j))**n
!      c(n+1,6)=-(x(j+4)-x(j))**n
      b(n+1)=0d0      
      enddo

!		call LSLRG(5,c,9,b,1,sol)
!      call gaussj(c,6,9,b,1,1)
      call gaussj(c,5,9,b,1,1)

      ben(j)=0d0
      aln(j)=0d0
      alp(j)=b(1)
      bep(j)=b(2)
!      a2n(j)=b(6)
!      a1n(j)=b(5)
      a2n(j)=0d0
      a1n(j)=0d0
      a0(j)=b(3)
      a1p(j)=b(4)
      a2p(j)=b(5)
!************************************************************************
!************************************************************************
      j=2
      c(1,1)=0d0
      c(1,2)=0d0
      c(1,3)=0d0
      c(1,4)=1d0
      c(1,5)=1d0
      c(1,6)=1d0
      c(1,7)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=0d0
      c(2,4)=(x(j-1)-x(j))
      c(2,5)=0d0
      c(2,6)=(x(j+1)-x(j))
      c(2,7)=(x(j+2)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=2d0
      c(3,3)=2d0
      c(3,4)=-(x(j-1)-x(j))**2
      c(3,5)=0d0
      c(3,6)=-(x(j+1)-x(j))**2
      c(3,7)=-(x(j+2)-x(j))**2
      b(3)=-2d0

      do n=3,6
      c(n+1,1)=n*(n-1)*(x(j-1)-x(j))**(n-2)
      c(n+1,2)=n*(n-1)*(x(j+1)-x(j))**(n-2)
      c(n+1,3)=n*(n-1)*(x(j+2)-x(j))**(n-2)
      c(n+1,4)=-(x(j-1)-x(j))**n
      c(n+1,5)=0d0
      c(n+1,6)=-(x(j+1)-x(j))**n
      c(n+1,7)=-(x(j+2)-x(j))**n
      b(n+1)=0d0
      enddo

      call gaussj(c,7,9,b,1,1)
      
      ben(j)=0d0
      aln(j)=b(1)
      alp(j)=b(2)
      bep(j)=b(3)
      a2n(j)=0d0
      a1n(j)=b(4)
      a0(j)=b(5)
      a1p(j)=b(6)
      a2p(j)=b(7)

!*************************************************************************
!*************************************************************************
      j=ny-1
      c(1,1)=0d0
      c(1,2)=0d0
      c(1,3)=0d0
      c(1,4)=1d0
      c(1,5)=1d0
      c(1,6)=1d0
      c(1,7)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=0d0
      c(2,4)=(x(j-2)-x(j))
      c(2,5)=(x(j-1)-x(j))
      c(2,6)=0d0
      c(2,7)=(x(j+1)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=2d0
      c(3,3)=2d0
      c(3,4)=-(x(j-2)-x(j))**2
      c(3,5)=-(x(j-1)-x(j))**2
      c(3,6)=0d0
      c(3,7)=-(x(j+1)-x(j))**2
      b(3)=-2d0

      do n=3,6
      c(n+1,1)=n*(n-1)*(x(j-2)-x(j))**(n-2)
      c(n+1,2)=n*(n-1)*(x(j-1)-x(j))**(n-2)
      c(n+1,3)=n*(n-1)*(x(j+1)-x(j))**(n-2)
      c(n+1,4)=-(x(j-2)-x(j))**n
      c(n+1,5)=-(x(j-1)-x(j))**n
      c(n+1,6)=0d0
      c(n+1,7)=-(x(j+1)-x(j))**n
      b(n+1)=0d0
      enddo
      
      
      call gaussj(c,7,9,b,1,1)
      
      ben(j)=b(1)
      aln(j)=b(2)
      alp(j)=b(3)
      bep(j)=0d0
      a2n(j)=b(4)
      a1n(j)=b(5)
      a0(j)=b(6)
      a1p(j)=b(7)
      a2p(j)=0d0
      
!*************************************************************************
!*************************************************************************	
      j=ny
      c(1,1)=0d0
      c(1,2)=0d0
      c(1,3)=1d0
      c(1,4)=1d0
      c(1,5)=1d0
!      c(1,5)=1d0
!      c(1,6)=1d0
      b(1)=0d0
      
      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=0d0
      c(2,4)=(x(j-1)-x(j))
      c(2,5)=(x(j-2)-x(j))
!      c(2,5)=(x(j-3)-x(j))
!      c(2,6)=(x(j-4)-x(j))
      b(2)=0d0
      
      c(3,1)=2d0
      c(3,2)=2d0
      c(3,3)=0d0
      c(3,4)=-(x(j-1)-x(j))**2
      c(3,5)=-(x(j-2)-x(j))**2
!      c(3,5)=-(x(j-3)-x(j))**2
!      c(3,6)=-(x(j-4)-x(j))**2
      b(3)=-2d0

!      do n=3,5
      do n=3,4
      c(n+1,1)=n*(n-1)*(x(j-1)-x(j))**(n-2)
      c(n+1,2)=n*(n-1)*(x(j-2)-x(j))**(n-2)
      c(n+1,3)=0d0
      c(n+1,4)=-(x(j-1)-x(j))**n
      c(n+1,5)=-(x(j-2)-x(j))**n
!      c(n+1,5)=-(x(j-3)-x(j))**n
!      c(n+1,6)=-(x(j-4)-x(j))**n
      b(n+1)=0d0
      enddo
      
!      call gaussj(c,6,9,b,1,1)
      call gaussj(c,5,9,b,1,1)
      
      ben(j)=b(2)
      aln(j)=b(1)
      alp(j)=0d0
      bep(j)=0d0
      a2n(j)=b(5)
      a1n(j)=b(4)
      a0(j)=b(3)
!      a1p(j)=b(5)
!      a2p(j)=b(6)
      a1p(j)=0d0
      a2p(j)=0d0
      return

      end subroutine prederiv2

 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c prepares cfdiff pieces for using cfdiff subroutines and laps later
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derivadas(y,my,prem1,prem3,dt11,dt12,dt21,dt22,fmap)
      implicit none

      integer,intent(in) :: my
      double precision, intent(out),dimension(7,my):: prem1,
     .                  dt11,dt12
      double precision, intent(out),dimension(5,my):: prem3,
     .                  dt21,dt22
      double precision, intent(in), dimension(my) :: y
      double precision, intent(out), dimension(my) :: fmap

      real(8) :: d11(my,7),d12(my,7),d21(my,5),d22(my,5)
      integer i,j
!   ---  start doing something -------------

      call pred1m6o(my,d11(1,1),d11(1,2),d11(1,3),d11(1,4),d11(1,5), 
     $      d11(1,6),d11(1,7),d12(1,1),d12(1,2),d12(1,3), 
     $      d12(1,4),d12(1,5),d12(1,6),d12(1,7))

      call prederiv2(my,d21(1,1),d21(1,2),d21(1,3),d21(1,4),d21(1,5),
     $      d22(1,1),d22(1,2),d22(1,3),d22(1,4),d22(1,5),y)

!       Copiamos, transponemos, y calculamos LU

      do j=1,my
         do i=1,7
            prem1(i,j) = d11(j,i)
            dt11(i,j) = d11(j,i) !aaf
            dt12(i,j)  = d12(j,i)
         enddo
      enddo

      do j=1,my
         do i=1,5
            prem3(i,j) = d21(j,i)
            dt21(i,j)  = d21(j,i)
            dt22(i,j)  = d22(j,i)
         enddo
      enddo

!     Calculo de traspuestas

      call bandec7(prem1,my)
      call bandec(prem3,my) 

!     Calculo de fmap
       fmap(1) =dt12(4,1)*y(1)+dt12(5,1)*y(2) + 
     $          dt12(6,1)*y(3)+dt12(7,1)*y(4)
       fmap(2) =dt12(3,2)*y(1)+dt12(4,2)*y(2) + 
     $          dt12(5,2)*y(3)+dt12(6,2)*y(4) + 
     $          dt12(7,2)*y(5)
       fmap(3) =dt12(2,3)*y(1)+dt12(3,3)*y(2) + 
     $          dt12(4,3)*y(3)+dt12(5,3)*y(4) + 
     $          dt12(6,3)*y(5)+dt12(7,3)*y(6)
      do j=4,my-3
      fmap(j)=dt12(1,j)*y(j-3)
         do i=2,7
               fmap(j) = fmap(j) + dt12(i,j)*y(i+j-4)
         enddo
      enddo
      fmap(my-2)=dt12(1,my-2)*y(my-5)+dt12(2,my-2)*y(my-4)+ 
     $           dt12(3,my-2)*y(my-3)+dt12(4,my-2)*y(my-2)+  
     $           dt12(5,my-2)*y(my-1)+dt12(6,my-2)*y(my)  
      fmap(my-1)=dt12(1,my-1)*y(my-4)+dt12(2,my-1)*y(my-3)+ 
     $           dt12(3,my-1)*y(my-2)+dt12(4,my-1)*y(my-1)+ 
     $           dt12(5,my-1)*y(my) 
      fmap(my)=  dt12(1,my  )*y(my-3)+dt12(2,my  )*y(my-2)+ 
     $           dt12(3,my  )*y(my-1)+dt12(4,my  )*y(my)  
      call banbks7(prem1,my,fmap)
      do j=1,my
         fmap(j) = 1d0/fmap(j)
      enddo

      !!!call endderiv()

      end subroutine derivadas
 

!------------------------ Linear Algebra Routines ----------------------

!***********************************************************************
!    Subroutines that work on N-diagonal 
!    (Numerical Recipes in FORTRAN)
!      bandec, banbks
!***********************************************************************

!     Adapted to liso

      subroutine banbks(a,n,b)
      INTEGER, intent(in) :: n
      double precision, intent(in) :: a(5,n)
      double precision, intent(inout) :: b(n)
      INTEGER i,k
 
      do k=1,n-2
         b(k+1) = b(k+1)-a(4,k)*b(k)
         b(k+2) = b(k+2)-a(5,k)*b(k)
      enddo
      b(n) = b(n)- a(4,n-1)*b(n-1)

!     back substitution

      b(n) = b(n)*a(1,n)
      b(n-1) = (b(n-1)-a(2,n-1)*b(n))*a(1,n-1)
      do i=n-2,1,-1
         b(i) = (b(i)-a(2,i)*b(1+i)-a(3,i)*b(2+i))*a(1,i)
      enddo

      return
      end subroutine banbks


      subroutine bandec(a,n)
      INTEGER,intent(in) :: n
      double precision,intent(inout):: a(5,n)
      INTEGER j,k

      do j=1,3
         a(j,1)=a(j+2,1)
      enddo
      do j=1,4
         a(j,2)=a(j+1,2)
      enddo


      do k=1,n-2
         a(1,k)   = 1d0/a(1,k)

         a(4,k)   = a(1,k+1)*a(1,k)

         a(1,k+1) = a(2,k+1)-a(4,k)*a(2,k)
         a(2,k+1) = a(3,k+1)-a(4,k)*a(3,k)
         a(3,k+1) = a(4,k+1)
       !
         a(5,k)   = a(1,k+2)*a(1,k)

         a(1,k+2) = a(2,k+2)-a(5,k)*a(2,k)
         a(2,k+2) = a(3,k+2)-a(5,k)*a(3,k)
         a(3,k+2) = a(4,k+2)
         a(4,k+2) = a(5,k+2)
      enddo

!     k=n-1

      a(1,n-1) = 1d0/a(1,n-1)

      a(4,n-1)=a(1,n)*a(1,n-1)

      a(1,n) = a(2,n)-a(4,n-1)*a(2,n-1)
      a(2,n) = a(3,n)-a(4,n-1)*a(3,n-1)
      a(3,n) = a(4,n)

!     the next loop will be used in banbk

      a(1,n)=1d0/a(1,n)


      end subroutine bandec


      subroutine banbks7(a,n,b)
      integer,intent(in) :: n
      double precision,intent(in) :: a(7,n)
      double precision,intent(inout):: b(n)
      integer i,k
  
      do k=1,n-3
         b(k+1) = b(k+1)-a(5,k)*b(k)
         b(k+2) = b(k+2)-a(6,k)*b(k)
         b(k+3) = b(k+3)-a(7,k)*b(k)
      enddo

!     n-2

      b(n-1) = b(n-1)-a(5,n-2)*b(n-2)
      b(n)   = b(n)  -a(6,n-2)*b(n-2)

!     n-1

      b(n) = b(n)    -a(5,n-1)*b(n-1)

!     back substitution

      b(n) = b(n)*a(1,n)
      b(n-1) = (b(n-1)-a(2,n-1)*b(n))*a(1,n-1)
      b(n-2) = (b(n-2)-a(2,n-2)*b(n-1)-a(3,n-2)*b(n))*a(1,n-2)

      do i=n-3,1,-1
         b(i) = (b(i)-a(2,i)*b(1+i)-a(3,i)*b(2+i)-a(4,i)*b(3+i))*a(1,i)
      enddo

      return
      end subroutine banbks7

      subroutine bandec7(a,n)
      integer,intent(in) :: n
      double precision,intent(inout) :: a(7,n)
      integer j,k
!
      do j=1,4
         a(j,1)=a(j+3,1)
      enddo

      do j=1,5
         a(j,2)=a(j+2,2)
      enddo

      do j=1,6
         a(j,3)=a(j+1,3)
      enddo


!     LU

      do k=1,n-3
         a(1,k)   = 1d0/a(1,k)

         a(5,k)   = a(1,k+1)*a(1,k)
         a(1,k+1) = a(2,k+1)-a(5,k)*a(2,k)
         a(2,k+1) = a(3,k+1)-a(5,k)*a(3,k)
         a(3,k+1) = a(4,k+1)-a(5,k)*a(4,k)
         a(4,k+1) = a(5,k+1)

         a(6,k)   = a(1,k+2)*a(1,k)
         a(1,k+2) = a(2,k+2)-a(6,k)*a(2,k)
         a(2,k+2) = a(3,k+2)-a(6,k)*a(3,k)
         a(3,k+2) = a(4,k+2)-a(6,k)*a(4,k)
         a(4,k+2) = a(5,k+2)
         a(5,k+2) = a(6,k+2)

         a(7,k)   = a(1,k+3)*a(1,k)
         a(1,k+3) = a(2,k+3)-a(7,k)*a(2,k)
         a(2,k+3) = a(3,k+3)-a(7,k)*a(3,k)
         a(3,k+3) = a(4,k+3)-a(7,k)*a(4,k)
         a(4,k+3) = a(5,k+3)
         a(5,k+3) = a(6,k+3)
         a(6,k+3) = a(7,k+3)

      enddo

!     k=n-2
      a(1,n-2) = 1d0/a(1,n-2)

      a(5,n-2) = a(1,n-1)*a(1,n-2)

      a(1,n-1) = a(2,n-1)-a(5,n-2)*a(2,n-2)
      a(2,n-1) = a(3,n-1)-a(5,n-2)*a(3,n-2)
      a(3,n-1) = a(4,n-1)-a(5,n-2)*a(4,n-2)
      a(4,n-1) = a(5,n-1)

      a(6,n-2) = a(1,n)*a(1,n-2)

      a(1,n)   = a(2,n)-a(6,n-2)*a(2,n-2)
      a(2,n)   = a(3,n)-a(6,n-2)*a(3,n-2)
      a(3,n)   = a(4,n)-a(6,n-2)*a(4,n-2)
      a(4,n)   = a(5,n)
      a(5,n)   = a(6,n)
  
!     k=n-1

      a(1,n-1) = 1d0/a(1,n-1)
      a(5,n-1) = a(1,n)*a(1,n-1)

      a(1,n)   = a(2,n)-a(5,n-1)*a(2,n-1)
      a(2,n)   = a(3,n)-a(5,n-1)*a(3,n-1)
      a(3,n)   = a(4,n)-a(5,n-1)*a(4,n-1)
      a(4,n)   = a(5,n)

!     the next loop will be used in banbk

      a(1,n)=1d0/a(1,n)

      return
      end subroutine bandec7


!**************************************************
!*    Subroutine to resolve linear system
!*    (Numerical Recipes in FORTRAN)
!**************************************************
      subroutine gaussj(a,n,np,b,m,mp)
      integer m,mp,n,np,NMAX
      real(8) :: a(np,np),b(np,mp)
      parameter (NMAX=50)
      integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX),ii
      real(8) :: big,dum,pivinv

      do j=1,n
         ipiv(j)=0
      enddo

      do i=1,n
         big=0d0
         do j=1,n
            if(ipiv(j).ne.1)then
               do k=1,n
                  if (ipiv(k).eq.0) then
                     if (abs(a(j,k)).ge.big)then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  else if (ipiv(k).gt.1) then
                     pause 'singular matrix in gaussj'
                  endif
               enddo
            endif
         enddo
  
         ipiv(icol)=ipiv(icol)+1
         if (irow.ne.icol) then
            do l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
            enddo
            do l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
            enddo
         endif
  
         indxr(i)=irow
         indxc(i)=icol
         if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
         pivinv=1d0/a(icol,icol)
         a(icol,icol)=1d0
         do l=1,n
            a(icol,l)=a(icol,l)*pivinv
         enddo
         do l=1,m
            b(icol,l)=b(icol,l)*pivinv
         enddo
         do ll=1,n
            if(ll.ne.icol)then
               dum=a(ll,icol)
               a(ll,icol)=0d0
               do l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
               enddo
               do l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
               enddo
            endif
         enddo
      enddo !done with big loop at top
  
      do l=n,1,-1
         if(indxr(l).ne.indxc(l))then
            do k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
            enddo
         endif
      enddo
  
      return
      end subroutine gaussj

!************************************************************************
!                                                                       !
!    Subroutine deryr                                                   !
!                                                                       !
!    computes the first derivative in y of u, with size(u)=(ny)         !
!                                                                       !
!    Input:                                                             !
!         u: Field to derivate. Assumed to be real of size ny           !
!    Output:                                                            !
!        du: First derivative of u                                      !
!                                                                       !
!************************************************************************
      subroutine deryr(u,du,my,dt12,prem1,fmap)
      implicit none

      integer, intent(in) :: my
      double precision, intent(in),dimension(my) :: fmap,u
      double precision, intent(in),dimension(7,my):: prem1,
     .                  dt12

      double precision, intent(out),dimension(my) :: du
      
      integer i,j
      double precision, dimension(my) :: fwk

      fwk(1)=   dt12(4,1)*u(1)+dt12(5,1)*u(2) + 
     $    dt12(6,1)*u(3)+dt12(7,1)*u(4)
    
      fwk(2)=   dt12(3,2)*u(1)+dt12(4,2)*u(2) + 
     $    dt12(5,2)*u(3)+dt12(6,2)*u(4) + 
     $    dt12(7,2)*u(5)
    
      fwk(3)=   dt12(2,3)*u(1)+dt12(3,3)*u(2) + 
     $    dt12(4,3)*u(3)+dt12(5,3)*u(4) + 
     $    dt12(6,3)*u(5)+dt12(7,3)*u(6)

      do j=4,my-3
       fwk(j)=dt12(1,j)*u(j-3)
       do i=2,7
          fwk(j)=fwk(j) + dt12(i,j)*u(i+j-4)
       enddo
      enddo

      fwk(my-2)=dt12(1,my-2)*u(my-5)+dt12(2,my-2)*u(my-4)+ 
     $     dt12(3,my-2)*u(my-3)+dt12(4,my-2)*u(my-2)+ 
     $     dt12(5,my-2)*u(my-1)+dt12(6,my-2)*u(my)
    
      fwk(my-1)=dt12(1,my-1)*u(my-4)+dt12(2,my-1)*u(my-3)+ 
     $     dt12(3,my-1)*u(my-2)+dt12(4,my-1)*u(my-1)+ 
     $     dt12(5,my-1)*u(my)

      fwk(my)=  dt12(1,my  )*u(my-3)+dt12(2,my  )*u(my-2)+ 
     $     dt12(3,my  )*u(my-1)+dt12(4,my  )*u(my)

      call banbks7(prem1,my,fwk)

      do j=1,my
         du(j)=fwk(j)*fmap(j)
      enddo

      end subroutine deryr


!************************************************************************
!                                                                       !
!    Subroutine deryyr                                                  !
!                                                                       !
!    Computes the second derivative in y of u, with size(u)=(ny)        !
!                                                                       !
!    Input:                                                             !
!         u: Field to derivate. Assumed to be real of size ny           !
!                                                                       !
!    Output:                                                            !
!        fwk: nonlinear term of Modes00 equation                        !
!                                                                       !
!************************************************************************

      subroutine deryyr(u,du,my,dt22,prem3)
      implicit none

      integer, intent(in) :: my
      double precision, intent(in),dimension(my) :: u
      double precision, intent(in),dimension(5,my):: prem3,
     .                  dt22

      double precision, intent(out),dimension(my) :: du
 
      integer i,j
      real(8),dimension(my) :: fwk

      fwk(1) =  dt22(3,1)*u(1) + dt22(4,1)*u(2) + 
     $          dt22(5,1)*u(3)

      fwk(2) =  dt22(2,2)*u(1) + dt22(3,2)*u(2) + 
     $          dt22(4,2)*u(3) + dt22(5,2)*u(4)

      do j=3,my-2
         fwk(j)=dt22(1,j)*u(j-2)
         do i=2,5
            fwk(j) = fwk(j) + dt22(i,j)*u(i+j-3)
         enddo
      enddo

      fwk(my-1)=dt22(1,my-1)*u(my-3) + dt22(2,my-1)*u(my-2)+ 
     $          dt22(3,my-1)*u(my-1) + dt22(4,my-1)*u(my  )
    
      fwk(my)=  dt22(1,my  )*u(my-2) + dt22(2,my  )*u(my-1)+ 
     $          dt22(3,my  )*u(my)

      call banbks(prem3,my,fwk)
! ---- first derivative ----
      do j=1,my
          du(j) = fwk(j)
      enddo

      end subroutine deryyr


!***********************************************************************!
! aaf 2014-05
! Integrates (once) du in order to give u 
! all in real 8
!***********************************************************************!
      subroutine inty8(f,du,my,fmap,dt11,dt12,prem1)
      implicit none

      integer, intent(in) :: my
      double precision, intent(in),dimension(my) :: fmap,du
      double precision, intent(in),dimension(7,my):: dt11,dt12,prem1

      double precision, intent(out),dimension(my) :: f
 


! ----------------------- workspace ---------------------------
      integer i,j
      real(8),dimension(my) :: vpr
      real(8)               :: zero, suma
      real(8),dimension(7,my)   :: wk2

!-----------------Program------------------------------------!

      zero = 0d0
! calculate v :
      !do j=1,my
      do j=1,my
       do i=1,7
          wk2(i,j)=dt12(i,j)
       enddo
      enddo
!      !BC @ j=1 --> f=0.0
!      wk2(1,1) = zero
!      wk2(2,1) = zero
!      wk2(3,1) = zero
!      wk2(4,1) = 1d0
!      wk2(5,1) = zero
!      wk2(6,1) = zero
!      wk2(7,1) = zero
!     
      !create bandwitch matrix 
      call bandec7(wk2,my)

      !save the rhs into f
      do j=1,my
       !f(j)=du(j)
       f(j)=du(j)/fmap(j)
      enddo

      !prepare rhs multiplying by dt11
      vpr(1   ) =   dt11(4,1   )*f(1   ) + dt11(5,1   )*f(2) 
     $            + dt11(6,1   )*f(3   ) + dt11(7,1   )*f(4)
      !setting value of the integral
!      !setting value of the integral
!      vpr(1   ) =   0d0 
!                     
      vpr(2   ) =   dt11(3,2   )*f(1   ) + dt11(4,2   )*f(2) 
     $            + dt11(5,2   )*f(3   ) + dt11(6,2   )*f(4) 
     $            + dt11(7,2   )*f(5   )

      vpr(3   ) =   dt11(2,3   )*f(1   ) + dt11(3,3   )*f(2)  
     $            + dt11(4,3   )*f(3   ) + dt11(5,3   )*f(4) 
     $            + dt11(6,3   )*f(5   ) + dt11(7,3   )*f(6)

      do j=4,my-3
         vpr(j)=dt11(1,j)*f(j-3)
         do i=2,7
            vpr(j) = vpr(j) + dt11(i,j)*f(i+j-4)
         enddo
      enddo
    
      vpr(my-2) =    dt11(1,my-2)*f(my-5) + dt11(2,my-2)*f(my-4)  
     $             + dt11(3,my-2)*f(my-3) + dt11(4,my-2)*f(my-2)
     $             + dt11(5,my-2)*f(my-1) + dt11(6,my-2)*f(my)    
  
      vpr(my-1) =    dt11(1,my-1)*f(my-4) + dt11(2,my-1)*f(my-3) 
     $             + dt11(3,my-1)*f(my-2) + dt11(4,my-1)*f(my-1) 
     $             + dt11(5,my-1)*f(my)
          
      vpr(my)   =    dt11(1,my  )*f(my-3) + dt11(2,my  )*f(my-2) 
     $             + dt11(3,my  )*f(my-1) + dt11(4,my  )*f(my)

!Inverting dt12 and premultiplying by dt11*f gives integration
      !call banbks7(prem1,my,vpr)     
      call banbks7(wk2,my,vpr)     
!Now we put zero on y=-Ly and make rest half of rhov00 @ Ly
       do j=1,my
          f(j)  = vpr(j)
         !u(j)  = vpr(j) - (vpr(1)+vpr(my))/2.0d0
         !u(j)  = vpr(j) - vpr(1) - (vpr(my)-vpr(1))/2.0d0
      enddo 
!This way we make -DM/2 @ y=-Ly, and DM/2 @ y = Ly
      

      end subroutine inty8
!---------------------------------------------------------------! 






!***********************************************************************!
! aaf 2014-05
! Integrates (once) du in order to give u 
! all in real 8
!***********************************************************************!
      subroutine inty80(f,du,my,fmap,dt11,dt12,prem1)
      implicit none

      integer, intent(in) :: my
      double precision, intent(in),dimension(my) :: fmap,du
      double precision, intent(in),dimension(7,my):: dt11,dt12,prem1

      double precision, intent(out),dimension(my) :: f
 


! ----------------------- workspace ---------------------------
      integer i,j
      real(8),dimension(my) :: vpr
      real(8)               :: zero, suma
      real(8),dimension(7,my)   :: wk2

!-----------------Program------------------------------------!

      zero = 0d0
! calculate v :
      !do j=1,my
      do j=1,my
       do i=1,7
          wk2(i,j)=dt12(i,j)
       enddo
      enddo

!      !BC @ j=1 --> f=0.0
!      wk2(1,1) = zero
!      wk2(2,1) = zero
!      wk2(3,1) = zero
!      wk2(4,1) = 1d0
!      wk2(5,1) = zero
!      wk2(6,1) = zero
!      wk2(7,1) = zero

      call bandec7(wk2,my)
     
      !create bandwitch matrix 

      !save the rhs into f
      do j=1,my
       !f(j)=du(j)
       f(j)=du(j)/fmap(j)
      enddo

      !prepare rhs multiplying by dt11
      vpr(1   ) =   dt11(4,1   )*f(1   ) + dt11(5,1   )*f(2) 
     $            + dt11(6,1   )*f(3   ) + dt11(7,1   )*f(4)
!                     
      vpr(2   ) =   dt11(3,2   )*f(1   ) + dt11(4,2   )*f(2) 
     $            + dt11(5,2   )*f(3   ) + dt11(6,2   )*f(4) 
     $            + dt11(7,2   )*f(5   )

      vpr(3   ) =   dt11(2,3   )*f(1   ) + dt11(3,3   )*f(2)  
     $            + dt11(4,3   )*f(3   ) + dt11(5,3   )*f(4) 
     $            + dt11(6,3   )*f(5   ) + dt11(7,3   )*f(6)

      do j=4,my-3
         vpr(j)=dt11(1,j)*f(j-3)
         do i=2,7
            vpr(j) = vpr(j) + dt11(i,j)*f(i+j-4)
         enddo
      enddo
    
      vpr(my-2) =    dt11(1,my-2)*f(my-5) + dt11(2,my-2)*f(my-4)  
     $             + dt11(3,my-2)*f(my-3) + dt11(4,my-2)*f(my-2)
     $             + dt11(5,my-2)*f(my-1) + dt11(6,my-2)*f(my)    
  
      vpr(my-1) =    dt11(1,my-1)*f(my-4) + dt11(2,my-1)*f(my-3) 
     $             + dt11(3,my-1)*f(my-2) + dt11(4,my-1)*f(my-1) 
     $             + dt11(5,my-1)*f(my)
          
      vpr(my)   =    dt11(1,my  )*f(my-3) + dt11(2,my  )*f(my-2) 
     $             + dt11(3,my  )*f(my-1) + dt11(4,my  )*f(my)

!Inverting dt12 and premultiplying by dt11*f gives integration
      !call banbks7(prem1,my,vpr)     
      call banbks7(wk2,my,vpr)     
!Now we put zero on y=-Ly and make rest half of rhov00 @ Ly
       do j=1,my
          !f(j)  = vpr(j)
          f(j)  = vpr(j) - vpr(1)
         !u(j)  = vpr(j) - (vpr(1)+vpr(my))/2.0d0
         !u(j)  = vpr(j) - vpr(1) - (vpr(my)-vpr(1))/2.0d0
      enddo 
!This way we make -DM/2 @ y=-Ly, and DM/2 @ y = Ly
      

      end subroutine inty80
!---------------------------------------------------------------! 





