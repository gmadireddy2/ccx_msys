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
      subroutine specificheats(inpc,textpart,shcon,nshcon,
     &  nmat,ntmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,ier)
!
!     reading the input deck: *SPECIFIC HEAT
!	  ntmat --> number of the temperature dependent specific props
!	  nmat  --> number of different materials
!     shcon(0,j,i)       temperature at temperature point j of material i
!     shcon(1,j,i)       specific heat at constant pressure at the
!                        temperature point j of material i
!     shcon(2,j,i)       dynamic viscosity at the temperature point j of
!                        material i
!     shcon(3,1,i)       specific gas constant of material i
!     nshcon(i)          number of temperature data points for the specific
!                        heat of material i

      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nshcon(*),nmat,ntmat,ntmat_,istep,istat,n,ipoinpc(0:*),
     &  key,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),i,ier
!
      real*8 shcon(0:3,ntmat_,*)
!
      ntmat=0
!
!      if((istep.gt.0).and.(irstrt(1).ge.0)) then
!         write(*,*) 
!     &    '*ERROR reading *SPECIFIC HEAT: *SPECIFIC HEAT should be'
!         write(*,*) '  placed before all step definitions'
!         ier=1
!         return
!      endif
!
!      if(nmat.eq.0) then
!         write(*,*) 
!     &    '*ERROR reading *SPECIFIC HEAT: *SPECIFIC HEAT should be'
!         write(*,*) '  preceded by a *MATERIAL card'
!         ier=1
!         return
!      endif
!
      do i=2,n
         write(*,*) 
     &      '*WARNING reading *SPECIFIC HEAT: parameter not recognized:'
		 write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*SPECIFIC HEAT%")
      enddo
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
		 if((istat.lt.0).or.(key.eq.1)) return
         ntmat=ntmat+1
         nshcon(nmat)=ntmat
         if(ntmat.gt.ntmat_) then
			write(*,*) '*ERROR reading *SPECIFIC HEAT: increase ntmat_'
            ier=1
            return
         endif
         read(textpart(1)(1:20),'(f20.0)',iostat=istat) 
     &        shcon(1,ntmat,nmat)
		 if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*SPECIFIC HEAT%",ier)
            return
         endif
         read(textpart(2)(1:20),'(f20.0)',iostat=istat) 
     &        shcon(0,ntmat,nmat)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*SPECIFIC HEAT%",ier)
            return
         endif
      enddo
      return
      end

