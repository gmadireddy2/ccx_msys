/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2021 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static ITG *nface1,*ielmatf1,*ntmat1_,*mi1,*ielfa1,*ipnei1,*nef1,*num_cpus1,
    *nface1,*neij1;

static double *vfa1,*shcon1,*vel1,*flux1,*xxi1,*xle1,*gradpel1,*gradtel1;

void hrr_udmain(ITG *nface,double *vfa,double *shcon,ITG *ielmatf,ITG *ntmat_,
			  ITG *mi,ITG *ielfa,ITG *ipnei,double *vel,ITG *nef,
			  double *flux,ITG *num_cpus,double *xxi,double *xle,
                          double *gradpel,double *gradtel,ITG *neij){

    ITG i;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;
    
    pthread_t tid[*num_cpus];

    /* calculation of the density at the cell centers */
    
    vfa1=vfa;shcon1=shcon;ielmatf1=ielmatf;ntmat1_=ntmat_;mi1=mi;ielfa1=ielfa;
    ipnei1=ipnei;vel1=vel;nef1=nef;flux1=flux;num_cpus1=num_cpus;nface1=nface;
    xxi1=xxi;xle1=xle;gradpel1=gradpel;gradtel1=gradtel;neij1=neij;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)hrr_ud1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of calcgammav1 */

void *hrr_ud1mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(hrr_ud,(vfa1,shcon1,ielmatf1,ntmat1_,mi1,ielfa1,
			      ipnei1,vel1,nef1,flux1,&nfacea,&nfaceb,
                              xxi1,xle1,gradpel1,gradtel1,neij1));

    return NULL;
}