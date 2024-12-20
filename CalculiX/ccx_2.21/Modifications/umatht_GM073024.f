C>>  CalculiX UMATHT Subroutine
C>>  -Guru Madireddy 07/17/2024
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
!
!
!     heat transfer material subroutine
!
!     INPUT:
!
!     statev(nstatv)      internal state variables at the start
!                         of the increment
!     temp                temperature at the start of the increment
!     dtemp               increment of temperature
!     dtemdx(ntgrd)       current values of the spatial gradients of the 
!                         temperature
!     time(1)             step time at the beginning of the increment
!     time(2)             total time at the beginning of the increment
!     dtime               time increment
!     predef              not used
!     dpred               not used
!     cmname              material name
!     ntgrd               number of spatial gradients of temperature
!     nstatv              number of internal state variables as defined
!                         on the *DEPVAR card
!     props(nprops)       user defined constants defined by the keyword
!                         card *USER MATERIAL,TYPE=THERMAL
!     nprops              number of user defined constants, as specified
!                         on the *USER MATERIAL,TYPE=THERMAL card
!     coords              global coordinates of the integration point
!     pnewd               not used
!     noel                element number
!     npt                 integration point number
!     layer               not used
!     kspt                not used
!     kstep               not used
!     kinc                not used
!     vold(0..4,1..nk)    solution field in all nodes
!                         0: temperature
!                         1: displacement in global x-direction
!                         2: displacement in global y-direction
!                         3: displacement in global z-direction
!                         4: static pressure
!     co(3,1..nk)         coordinates of all nodes
!                         1: coordinate in global x-direction
!                         2: coordinate in global y-direction
!                         3: coordinate in global z-direction
!     lakonl              element label
!     konl(1..20)         nodes belonging to the element
!     ipompc(1..nmpc))   ipompc(i) points to the first term of
!                        MPC i in field nodempc
!     nodempc(1,*)       node number of a MPC term
!     nodempc(2,*)       coordinate direction of a MPC term
!     nodempc(3,*)       if not 0: points towards the next term
!                                  of the MPC in field nodempc
!                        if 0: MPC definition is finished
!     coefmpc(*)         coefficient of a MPC term
!     nmpc               number of MPC's
!     ikmpc(1..nmpc)     ordered global degrees of freedom of the MPC's
!                        the global degree of freedom is
!                        8*(node-1)+direction of the dependent term of
!                        the MPC (direction = 0: temperature;
!                        1-3: displacements; 4: static pressure;
!                        5-7: rotations)
!     ilmpc(1..nmpc)     ilmpc(i) is the MPC number corresponding
!                        to the reference number in ikmpc(i)   
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedomm per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!
!     OUTPUT:
!     
!     u                   not used
!     dudt                not used
!     dudg(ntgrd)         not used
!     flux(ntgrd)         heat flux at the end of the increment
!     dfdt(ntgrd)         not used
!     dfdg(ntgrd,ntgrd)   variation of the heat flux with respect to the 
!                         spatial temperature gradient
!     statev(nstatv)      internal state variables at the end of the
!                         increment
!
      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,
     &  statev,temp,dtemp,dtemdx,time,dtime,predef,dpred,
     &  cmname,ntgrd,nstatv,props,nprops,coords,pnewdt,
     &  noel,npt,layer,kspt,kstep,kinc,vold,co,lakonl,konl,
     &  ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,mi)
      use, intrinsic :: iso_c_binding
	  implicit real*8(a-h,o-z), integer*4(i-n)
!
      character*8 lakonl
      character*80 cmname
!
      integer ntgrd,nstatv,nprops,noel,npt,layer,kspt,kstep,kinc,
     &  konl(20),ipompc(*),nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),mi(*)
!
      real*8 u,dudt,dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     &  statev(nstatv),pnewdt,temp,dtemp,dtemdx(ntgrd),time(2),dtime,
     &  predef,dpred,props(nprops),coords(3),dfdg(ntgrd,ntgrd),
     &  vold(0:mi(2),*),co(3,*),coefmpc(*), i, cond, j
!
!     insert here your code
!     
C>> Material's thermal properties common block
      COMMON /THM_DATA/ CPA(2),AKA(2),specht
C
!      CHARACTER*256 FNAME
      
!	  FNAME = 'C:\\Users\\gmadireddy\\CalculiX\\Charlie\\A.THM'
! 	  A.THM is used as short version of STEEL_THERM_DATA.THM
!	  because the longer name is not compiling due to line charecters limit
!     need to figure out how to split the longer file locations - 07/17/24


 	  call abaqus_sub(u,dudt,dudg,flux,dfdt,dfdg,statev,temp,
     $     dtemp,dtemdx,time,dtime,predef,dpred,cmname,ntgrd,nstatv,
     $     props,nprops,coords,pnewdt,noel,npt,layer,kspt,kstep,kinc)
	 
!	  OPEN(UNIT=16, FILE=FNAME, STATUS='OLD')
C	  WRITE(*,*) 'File opened correctly'
C	  
!	  DO I=1,3
!		READ(16,*)
!	  ENDDO
!	  READ(16,*)
!	  READ(16,*) (CPA(J), J=1,2)
c
!	  DO I=1, 3
!         READ(16,*)
!	  END DO
!	  READ(16,*)
!	  READ(16,*) (AKA(J), J=1,2)
!	  CLOSE(16)
C
C	  WRITE(*,*) 'Material Data is Read ikn Successfully!'

C>>
!	  cond = AKA(1)+AKA(2)*temp    !thermal conductivity
!      specht = CPA(1)+CPA(2)*temp  !specific heat
!	  print *, specht, 'umatht.f'
c
!      dudt = specht
!      du = dudt*dtemp
!      u = u+du
!	  print *, 'UMATHT', specht
	  
c
!      DO i=1, ntgrd
!         flux(i) = -cond*dtemdx(i)
!      END DO
	  
!	  DO i=1, ntgrd
!         DO j=1, ntgrd
!            dfdg(i,j)=0.d0   !The non-diagonal elements need to specified as 0 explicitly
!         END DO
!         dfdg(i,i)=cond      !Here is the conductivity, in abaqus it is negative conductvity
!      END DO
C
      RETURN
      END
