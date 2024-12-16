!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2022 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!
      implicit none
!
      character*80 cmname
!
      integer ndi,nshr,ntens,nstatv,nprops,noel,npt,layer,kspt,
     &  kstep,kinc
!
      real*8 stress(ntens),statev(nstatv),
     &  ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     &  stran(ntens),dstran(ntens),time(2),celent,
     &  props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),
     &  sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,predef,dpred,
     &  pnewdt
!
!     START EXAMPLE LINEAR ELASTIC MATERIAL
!
      integer i,j
      real*8 E,nu,lamda,mu,am1,am2
!
      write(*,*) 'noel,npt ',noel,npt
      write(*,*) 'stress ',(stress(i),i=1,6)
      write(*,*) 'stran ',(stran(i),i=1,6)
      write(*,*) 'dstran ',(dstran(i),i=1,6)
      write(*,*) 'drot ',((drot(i,j),i=1,3),j=1,3)
      E=props(1)
      nu=props(2)
      lamda=nu*E/(1.d0+nu)/(1.d0-2.d0*nu)
      mu=E/2.d0/(1.d0+nu)
      am1=lamda+2.d0*mu
      am2=mu
!
!     stress
!      
      stress(1)=stress(1)+am1*dstran(1)+lamda*(dstran(2)+dstran(3))
      stress(2)=stress(2)+am1*dstran(2)+lamda*(dstran(1)+dstran(3))
      stress(3)=stress(3)+am1*dstran(3)+lamda*(dstran(1)+dstran(2))
      stress(4)=stress(4)+am2*dstran(4)
      stress(5)=stress(5)+am2*dstran(5)
      stress(6)=stress(6)+am2*dstran(6)
!
!     stiffness
!
      print *,"-------Here in umat-------"	
      do i=1,6
         do j=1,6
            ddsdde(i,j)=0.d0
         enddo
      enddo
      ddsdde(1,1)=lamda+2.d0*mu
      ddsdde(1,2)=lamda
      ddsdde(2,1)=lamda
      ddsdde(2,2)=lamda+2.d0*mu
      ddsdde(1,3)=lamda
      ddsdde(3,1)=lamda
      ddsdde(2,3)=lamda
      ddsdde(3,2)=lamda
      ddsdde(3,3)=lamda+2.d0*mu
      ddsdde(4,4)=mu
      ddsdde(5,5)=mu
      ddsdde(6,6)=mu
!
!     END EXAMPLE LINEAR ELASTIC MATERIAL
!
      return
      end
