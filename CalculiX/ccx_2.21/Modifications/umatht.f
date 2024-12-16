      subroutine umatht(u, dudt, dudg, flux, dfdt, dfdg, statev, temp, 
     &  dtemp, dtemdx, time, dtime, predef, dpred, cmname, ntgrd, 
     &  nstatv, props, nprops, coords, pnewdt, noel, npt, layer, 
     &  kspt, kstep, kinc, vold, co, lakonl, konl, ipompc, nodempc, 
     &  coefmpc, nmpc, ikmpc, ilmpc, mi)



	 
      COMMON /THM_DATA/ specht,cond
!	  external abaqus_sub
	  WRITE(*,*) 'Calling abaqus_sub from umatht'
!	  call abaqus_sub(u,dudt,temp,dtemp,ntgrd,dfdg,flux,dtemdx)
	  WRITE(*,*) 'Returned from abaqus_sub in umatht'
	  	  

	  return
      end
