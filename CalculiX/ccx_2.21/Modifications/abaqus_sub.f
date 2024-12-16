! This is a sample subroutine for creating the framework to run
! Dante's subroutine as dll.
! This subroutine works as an example of how the Dante's subroutine need
! to be written for allowing it to be integrated with CalculiX.
! This subroutine need to be compiled with gfortran using cygwin. 
! Untested if the dll works if compiled by any other software
! -- Guru Madireddy 08/07/2024
c>> Abaqus UMATHT subroutine 

      subroutine abaqus_sub(u,dudt,temp,dtemp,ntgrd,dfdg,flux,dtemdx,
     &  statev,dudg,time,dtime,predef,dpred,cmname,nstatv,
     &  props,nprops,coords,pnewdt,noel,npt,layer,kspt,kstep,
     &  kinc,vold,co,lakonl,konl,ipompc,nodempc,coefmpc,nmpc,
     &  ikmpc,ilmpc,mi,specht,cond)
  
      implicit real*8(a-h,o-z), integer*4(i-n)
	  character*8 lakonl
	  character*80 cmname
	  integer ntgrd,nstatv,nprops,noel,npt,layer,kspt,kstep,kinc,
     &  konl(20),ipompc(*),nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),mi(*)

      real*8 u,dudt,dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     &  statev(nstatv),pnewdt,temp,dtemp,dtemdx(ntgrd),time(2),dtime,
     &  predef,dpred,props(nprops),coords(3),dfdg(ntgrd,ntgrd),
     &  vold(0:mi(2),*),co(3,*),coefmpc(*),cond,CPA(2),AKA(2)
	 
      character*256 FNAME
	  FNAME = 'C:\\Users\\gmadireddy\\CalculiX\\Charlie\\A.THM'
	  
	  OPEN(UNIT=16, FILE=FNAME, STATUS='OLD')
	  
	  DO I=1, 3
         READ(16,*)
	  END DO
	  READ(16,*)
      READ(16,*) (CPA(J), J=1,2)

	  DO I=1, 3
         READ(16,*)
	  END DO
	  READ(16,*)
	  READ(16,*) (AKA(J), J=1,2)
	  CLOSE(16)
	  
      cond = AKA(1)+AKA(2)*temp
	  specht = CPA(1)+CPA(2)*temp

      dudt = specht
      du = dudt*dtemp
      u = u+du
	  
	  
      do i=1, ntgrd
         flux(i) = -cond*dtemdx(i)
      end do
	  
	  DO i=1, ntgrd
         DO j=1, ntgrd
            dfdg(i,j)=0.d0   !The non-diagonal elements need to specified as 0 explicitly
         END DO
         dfdg(i,i)=cond      !Here is the conductivity, in abaqus it is negative conductvity
      END DO
	  return
      end
!
