!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine film(h,sink,temp,kstep,kinc,time,noel,npt,
     &  coords,jltyp,field,nfield,loadtype,node,area,vold,mi,
     &  ipkon,kon,lakon,iponoel,inoel,ielprop,prop,ielmat,
     &  shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon,
     &  ipobody,xbody,ibody,heatnod,heatfac)
!
!     user subroutine film
!
!
!     INPUT:
!
!     sink               most recent sink temperature
!     temp               current temperature value
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     noel               element number
!     npt                integration point number
!     coords(1..3)       global coordinates of the integration point
!     jltyp              loading face kode:
!                        11 = face 1 
!                        12 = face 2 
!                        13 = face 3 
!                        14 = face 4 
!                        15 = face 5 
!                        16 = face 6
!     field              currently not used
!     nfield             currently not used (value = 1)
!     loadtype           load type label
!     node               network node (only for forced convection)
!     area               area covered by the integration point
!     vold(0..4,1..nk)   actual solution field in all nodes; 
!                        for structural nodes:
!                        0: temperature
!                        1: displacement in global x-direction
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!                        4: static pressure
!                        for network nodes
!                        0: total temperature (at end nodes)
!                           = static temperature for liquids
!                        1: mass flow (at middle nodes)
!                        2: total pressure (at end nodes)
!                           = static pressure for liquids
!                        3: static temperature (at end nodes; only for gas)
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedom per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!     mi(3)              max # of layers in any element
!     ipkon(i)           points to the location in field kon preceding
!                        the topology of element i
!     kon(*)             contains the topology of all elements. The
!                        topology of element i starts at kon(ipkon(i)+1)
!                        and continues until all nodes are covered. The
!                        number of nodes depends on the element label
!     lakon(i)           contains the label of element i
!     iponoel(i)         the network elements to which node i belongs
!                        are stored in inoel(1,iponoel(i)),
!                        inoel(1,inoel(2,iponoel(i)))...... until
!                        inoel(2,inoel(2,inoel(2......)=0
!     inoel(1..2,*)      field containing the network elements
!     ielprop(i)         points to the location in field prop preceding
!                        the properties of element i
!     prop(*)            contains the properties of all network elements. The
!                        properties of element i start at prop(ielprop(i)+1)
!                        and continues until all properties are covered. The
!                        appropriate amount of properties depends on the
!                        element label. The kind of properties, their
!                        number and their order corresponds
!                        to the description in the user's manual,
!                        cf. the sections "Fluid Section Types"
!     ielmat(j,i)        contains the material number for element i
!                        and layer j
!     shcon(0,j,i)       temperature at temperature point j of material i
!     shcon(1,j,i)       specific heat at constant pressure at the
!                        temperature point j of material i
!     shcon(2,j,i)       dynamic viscosity at the temperature point j of
!                        material i
!     shcon(3,1,i)       specific gas constant of material i
!     nshcon(i)          number of temperature data points for the specific
!                        heat of material i
!     rhcon(0,j,i)       temperature at density temperature point j of 
!                        material i
!     rhcon(1,j,i)       density at the density temperature point j of
!                        material i
!     nrhcon(i)          number of temperature data points for the density
!                        of material i
!     ntmat_             maximum number of temperature data points for 
!                        any material property for any material
!     ncocon(1,i)        number of conductivity constants for material i
!     ncocon(2,i)        number of temperature data points for the 
!                        conductivity coefficients of material i
!     cocon(0,j,i)       temperature at conductivity temperature point
!                        j of material i
!     cocon(k,j,i)       conductivity coefficient k at conductivity
!                        temperature point j of material i
!     ipobody(1,i)       points to an entry in fields ibody and xbody 
!                        containing the body load applied to element i, 
!                        if any, else 0
!     ipobody(2,i)       index referring to the line in field ipobody
!                        containing a pointer to the next body load
!                        applied to element i, else 0
!     ibody(1,i)         code identifying the kind of body load i:
!                        -1,1=centrifugal, 2=gravity, 3=generalized gravity
!     ibody(2,i)         amplitude number for load i
!     ibody(3,i)         load case number for load i
!     xbody(1,i)         size of body load i
!     xbody(2..4,i)      for centrifugal loading: point on the axis,
!                        for gravity loading with known gravity vector:
!                          normalized gravity vector
!     xbody(5..7,i)      for centrifugal loading: normalized vector on the
!                          rotation axis
!
!     OUTPUT:
!
!     h(1)               magnitude of the film coefficient
!     h(2)               not used; please do NOT assign any value
!     sink               (updated) sink temperature (need not be
!                        defined for forced convection)
!     heatnod            extra heat flow going to the network node
!                        (zero if not specified)
!     heatfac            extra heat flow going to the structural face
!                        (zero if not specified)
!           
      implicit none
!
      character*8 lakon(*)
      character*20 loadtype
!
      integer kstep,kinc,noel,npt,jltyp,nfield,node,mi(*),ipkon(*),
     &  kon(*),iponoel(*),inoel(2,*),ielprop(*),ielmat(mi(3),*),ntmat_,
     &  nshcon(*),nrhcon(*),ncocon(2,*),nodem,indexprop,indexe,
     &  iel1,iel2,ielup,iit,imat,icase,itherm,ipobody(2,*),
     &  ibody(3,*)
!
      real*8 h(2),sink,time(2),coords(3),temp,field(nfield),area,
     &  vold(0:mi(2),*),prop(*),shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),
     &  cocon(0:6,ntmat_,*),rho,r,pt,Pr,xl,Tt,Ts,xflow,Tsold,Re,um,
     &  xks,xkappa,xlambda,f,cp,A,D,form_fact,xbody(7,*),heatnod,
     &  heatfac
!
!	  
	  
      if (time(2) .LE. 3600) then
		sink = 875
		h(1) = (7.142857142857143e-8) * temp + 0.00005
	  else if (time(2) .GT. 3610 .and. time(2) .LE. 4210) then
		sink = 65
		if (temp .LE. 20) then
			h(1) = 0.0005
		else if (temp .GT. 20 .and. temp .LE. 150) then
			h(1) = 0.0017
		else if (temp .GT. 150 .and. temp .LE. 300) then
			h(1) = 0.0036
		else if (temp .GT. 300 .and. temp .LE. 350) then
			h(1) = 0.0048
		else if (temp .GT. 350 .and. temp .LE. 400) then
			h(1) = 0.006
		else if (temp .GT. 400 .and. temp .LE. 450) then
			h(1) = 0.00675
		else if (temp .GT. 450 .and. temp .LE. 500) then
			h(1) = 0.008
		else if (temp .GT. 500 .and. temp .LE. 550) then
			h(1) = 0.008
		else if (temp .GT. 550 .and. temp .LE. 600) then
			h(1) = 0.0076
		else if (temp .GT. 600 .and. temp .LE. 650) then
			h(1) = 0.0051
		else if (temp .GT. 650 .and. temp .LE. 700) then
			h(1) = 0.0044
		else if (temp .GT. 700 .and. temp .LE. 750) then
			h(1) = 0.0033
		else if (temp .GT. 750 .and. temp .LE. 800) then
			h(1) = 0.0026
		else if (temp .GT. 800 .and. temp .LE. 850) then
			h(1) = 0.00145
		else if (temp .GT. 850 .and. temp .LE. 1000) then
			h(1) = 0.00176
		end if
	  end if

!
      return
      end

