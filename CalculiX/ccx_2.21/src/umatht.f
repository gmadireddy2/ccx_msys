      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,statev,temp,
     &  dtemp,dtemdx,time,dtime,predef,dpred,cmname,ntgrd,
     &  nstatv,props,nprops,coords,pnewdt,noel,npt,layer,
     &  kspt,kstep,kinc,vold,co,lakonl,konl,ipompc,nodempc,
     &  coefmpc,nmpc,ikmpc,ilmpc,mi)

	  COMMON /THM_DATA/ specht,cond
!	  Below calls abqsub.dll where the specht and other properties
!	  are calculated from Dante's material library. The dll returns
!	  all the variables calculated from abqsub subroutine.

	  integer ntgrd,nstatv,nprops
	  real*8 props(nprops),temp(2),statev(nstatv)
	  
	  call abqsub(u,dudt,temp,dtemp,ntgrd,dfdg,flux,dtemdx,
     &  statev,dudg,time,dtime,predef,dpred,cmname,nstatv,
     &  props,nprops,coords,pnewdt,noel,npt,layer,kspt,kstep,
     &  kinc,vold,co,lakonl,konl,ipompc,nodempc,coefmpc,nmpc,
     &  ikmpc,ilmpc,mi,specht,cond,dfdt)
	  

	  return
      end
