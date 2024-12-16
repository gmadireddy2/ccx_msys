!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine beamgeneralsections(inpc,textpart,set,istartset,
     &  iendset,ialset,nset,ielmat,matname,nmat,ielorien,orname,norien,
     &  thicke,ipkon,iponor,xnor,ixfree,offset,lakon,irstrt,istep,
     &  istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,mi,ielprop,nprop,
     &  nprop_,prop,nelcon)
!
!     reading the input deck: *BEAM GENERAL SECTION
!
      implicit none
!
      character*1 inpc(*)
      character*4 section
      character*8 lakon(*)
      character*80 matname(*),orname(*),material,orientation
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),mi(*),ielmat(mi(3),*),
     &  ipoinpc(0:*),ielorien(mi(3),*),ipkon(*),iline,ipol,inl,lprop,
     &  ipoinp(2,*),inp(3,*),nset,nmat,norien,istep,istat,n,key,i,j,k,l,
     &  imaterial,iorientation,ipos,m,iponor(2,*),ixfree,indexx,indexe,
     &  irstrt,ielprop(*),nprop_,nprop,npropstart,ndprop,ndpropread,
     &  nelcon(2,*)
!
      real*8 thicke(mi(3),*),thickness1,thickness2,p(3),xnor(*),
     &  offset(2,*),offset1,offset2,dd,prop(*)
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR reading *BEAM GENERAL SECTION:'
         write(*,*) '       *BEAM GENERAL SECTION should'
         write(*,*) '       be placed before all step definitions'
         call exit(201)
      endif
!
      offset1=0.d0
      offset2=0.d0
      orientation='                    
     &                           '
      section='    '
      ipos=1
      npropstart=nprop
!
      do i=2,n
         if(textpart(i)(1:9).eq.'MATERIAL=') then
            material=textpart(i)(10:89)
         elseif(textpart(i)(1:12).eq.'ORIENTATION=') then
            orientation=textpart(i)(13:92)
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(21:21)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         elseif(textpart(i)(1:8).eq.'SECTION=') then
            if(textpart(i)(9:12).eq.'PIPE') then
               section='PIPE'
               ndprop=2
            elseif(textpart(i)(9:11).eq.'BOX') then
               section='BOX'
               ndprop=6
            else
               write(*,*) 
     &           '*ERROR reading *BEAM GENERAL SECTION: unknown section'
               call exit(201)
            endif
         elseif(textpart(i)(1:8).eq.'OFFSET1=') then
            read(textpart(i)(9:28),'(f20.0)',iostat=istat) offset1
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*BEAM GENERAL SECTION%")
         elseif(textpart(i)(1:8).eq.'OFFSET2=') then
            read(textpart(i)(9:28),'(f20.0)',iostat=istat) offset2
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*BEAM GENERAL SECTION%")
         else
            write(*,*) '*WARNING reading *BEAM GENERAL SECTION:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*BEAM GENERAL SECTION%")
         endif
      enddo
!
!     check whether a sections was defined
!
      if(section.eq.'    ') then
         write(*,*) '*ERROR reading *BEAM GENERAL SECTION:'
         write(*,*) '       no section defined'
         call exit(201)
      endif
!
!     check for the existence of the set,the material and orientation
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         write(*,*) '*ERROR reading *BEAM GENERAL SECTION:'
         write(*,*) '       nonexistent material'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &"*BEAM GENERAL SECTION%")
         call exit(201)
      endif
      imaterial=i
!
      if(orientation.eq.'                    
     &                                 ') then
         iorientation=0
      elseif(nelcon(1,i).eq.2) then
         write(*,*) '*INFO reading *SOLID SECTION: an orientation'
         write(*,*) '      is for isotropic materials irrelevant'
         call inputinfo(inpc,ipoinpc,iline,
     &"*SOLID SECTION%")
         iorientation=0
      else
         do i=1,norien
            if(orname(i).eq.orientation) exit
         enddo
         if(i.gt.norien) then
            write(*,*)
     &   '*ERROR reading *BEAM GENERAL SECTION: nonexistent orientation'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &"*BEAM GENERAL SECTION%")
            call exit(201)
         endif
         iorientation=i
      endif
!
      do i=1,nset
         if(set(i).eq.elset) exit
      enddo
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*)'*ERROR reading *BEAM GENERAL SECTION: element set ',
     &      elset(1:ipos)
         write(*,*) '  has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &"*BEAM GENERAL SECTION%")
         call exit(201)
      endif
!
!     assigning the elements of the set the appropriate material,
!     orientation number, section and offset(s)
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            if(lakon(ialset(j))(1:4).ne.'B32R') then
               write(*,*) '*ERROR reading *BEAM GENERAL SECTION:'
               write(*,*) '       *BEAM GENERAL SECTION can'
               write(*,*) '       only be used for B32R elements.'
               write(*,*) '       Element ',ialset(j),' is not a B32R el
     &ement.'
               call exit(201)
            endif
            ielmat(1,ialset(j))=imaterial
            ielorien(1,ialset(j))=iorientation
            offset(1,ialset(j))=offset1
            offset(2,ialset(j))=offset2
            if(section.eq.'PIPE') then
               lakon(ialset(j))(8:8)='P'
            elseif(section.eq.'BOX') then
               lakon(ialset(j))(8:8)='B'
            endif
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if(lakon(k)(1:1).ne.'B') then
                  write(*,*) '*ERROR reading *BEAM GENERAL SECTION:'
                  write(*,*) '       *BEAM GENERAL SECTION can'
                  write(*,*) '       only be used for beam elements.'
                  write(*,*) '       Element ',k,' is not a beam element
     &.'
                  call exit(201)
               endif
               ielmat(1,k)=imaterial
               ielorien(1,k)=iorientation
               offset(1,k)=offset1
               offset(2,k)=offset2
               if(section.eq.'PIPE') then
                  lakon(k)(8:8)='P'
               elseif(section.eq.'BOX') then
                  lakon(k)(8:8)='B'
               endif
            enddo
         endif
      enddo
!
!     reading the properties
!
      lprop=0
      ndpropread=ndprop
      do j=1,(ndpropread-1)/8+1
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         do k=1,8
            lprop=lprop+1
            if(lprop.gt.ndpropread) exit
            read(textpart(k),'(f40.0)',iostat=istat)
     &           prop(nprop+lprop)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*BEAM GENERAL SECTION%")
         enddo
      enddo
      nprop=nprop+ndprop
!
      if(nprop.gt.nprop_) then
         write(*,*) 
     &       '*ERROR reading *BEAM GENERAL SECTION: increase nprop_'
         call exit(201)
      endif
!
!     calculating the dimensions of the rectangular parent beam
!
      if(section.eq.'PIPE') then
         thickness1=2.d0*prop(npropstart+1)
         thickness2=thickness1
      elseif(section.eq.'BOX') then
         thickness1=prop(npropstart+1)
         thickness2=prop(npropstart+2)
      endif
!
!     assigning the thickness and the properties to the elements
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            indexe=ipkon(ialset(j))
            do l=1,8
               thicke(1,indexe+l)=thickness1
               thicke(2,indexe+l)=thickness2
            enddo
            ielprop(ialset(j))=npropstart
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               indexe=ipkon(k)
               do l=1,8
                  thicke(1,indexe+l)=thickness1
                  thicke(2,indexe+l)=thickness2
               enddo
               ielprop(k)=npropstart
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) return
!
!     assigning normal direction 1 for the beam
!
      indexx=-1
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) p(1)
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*BEAM GENERAL SECTION%")
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) p(2)
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*BEAM GENERAL SECTION%")
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) p(3)
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*BEAM GENERAL SECTION%")
      dd=dsqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
      if(dd.lt.1.d-10) then
         write(*,*) 
     &    '*ERROR reading *BEAM GENERAL SECTION: normal in direction 1'
         write(*,*) '       has zero size'
         call exit(201)
      endif
      do j=1,3
         p(j)=p(j)/dd
      enddo
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            indexe=ipkon(ialset(j))
            do l=1,8
               if(indexx.eq.-1) then
                  indexx=ixfree
                  do m=1,3
                     xnor(indexx+m)=p(m)
                  enddo
                  ixfree=ixfree+6
               endif
               iponor(1,indexe+l)=indexx
            enddo
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               indexe=ipkon(k)
               do l=1,8
               if(indexx.eq.-1) then
                  indexx=ixfree
                  do m=1,3
                     xnor(indexx+m)=p(m)
                  enddo
                  ixfree=ixfree+6
               endif
               iponor(1,indexe+l)=indexx
               enddo
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

