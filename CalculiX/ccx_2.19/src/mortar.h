/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2021 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; either version 2 of    */
/*     the License,or (at your option) any later version.               */

/*     This program is distributed in the hope that it will be useful,  */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not,write to the Free Software       */
/*     Foundation,Inc.,675 Mass Ave,Cambridge,MA 02139,USA.         */

void add_rect(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1,ITG m_1,
	      double *au_2,ITG * irow_2,ITG * jq_2,ITG n_2,ITG m_2,
	      double **au_rp,ITG **irow_rp,ITG * jq_r,ITG *nzs);
       
void bdfill(ITG **irowbdp,ITG *jqbd,double **aubdp,ITG *nzsbd,
	    ITG **irowbdtilp,ITG *jqbdtil,double **aubdtilp,ITG *nzsbdtil,
	    ITG **irowbdtil2p,ITG *jqbdtil2,double **aubdtil2p,ITG *nzsbdtil2,
	    ITG **irowddp,ITG *jqdd,double **auddp,
	    ITG **irowddtilp,ITG *jqddtil,double **auddtilp,
	    ITG **irowddtil2p,ITG *jqddtil2,double **auddtil2p,
	    ITG **irowddinvp,ITG *jqddinv,double **auddinvp,
	    ITG *irowtloc,ITG *jqtloc,double *autloc,
	    ITG *irowtlocinv,ITG *jqtlocinv,double *autlocinv,
	    ITG *ntie,ITG *ipkon,ITG *kon,
	    char *lakon,ITG *nslavnode,ITG *nmastnode,ITG *imastnode,
	    ITG *islavnode,ITG *islavsurf,ITG *imastsurf,double *pmastsurf,
	    ITG *itiefac,char *tieset,ITG *neq,ITG *nactdof,double *co,double *vold,
	    ITG *iponoels,ITG *inoels,ITG *mi,double *gapmints,double *gap,
	    double* pslavsurf,double* pslavdual,double* pslavdualpg,ITG *nintpoint,double *slavnor,ITG *nk,
	    ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
	    ITG *ikmpc,ITG *ilmpc,
	    ITG *nmpc2,ITG *ipompc2,ITG *nodempc2,double *coefmpc2,
	    ITG *ikmpc2,ITG *ilmpc2,
	    ITG *nslavspc,ITG *islavspc,ITG *nsspc,ITG *nslavmpc,ITG *islavmpc,ITG *nsmpc,
	    ITG *nslavmpc2,ITG *islavmpc2,ITG *nsmpc2,
	    ITG *nmastspc,ITG *imastspc,ITG *nmspc,ITG *nmastmpc,ITG *imastmpc,ITG *nmmpc,
	    ITG *nmastmpc2,ITG *imastmpc2,ITG *nmmpc2,
	    ITG *iit,ITG *iinc,ITG *islavactdof,ITG *islavact,ITG *islavnodeinv,
	    double **Bdp,ITG **irowbp,ITG *jqb,
	    double **Bdhelpp,ITG **irowbhelpp,ITG *jqbhelp,
	    double **Ddp,ITG **irowdp,ITG *jqd,
	    double **Ddtilp,ITG **irowdtilp,ITG *jqdtil,
	    double **Bdtilp,ITG **irowbtilp,ITG *jqbtil,
	    double **Bpgdp,ITG **irowbpgp,ITG *jqbpg,
	    double **Dpgdp,ITG **irowdpgp,ITG *jqdpg,
	    double **Dpgdtilp,ITG **irowdpgtilp,ITG *jqdpgtil,
	    double **Bpgdtilp,ITG **irowbpgtilp,ITG *jqbpgtil,
	    ITG *iflagdualquad,ITG *ithermal);  
    
void buildtquad(ITG *ntie,ITG *ipkon,ITG *kon,ITG *nk,
		char *lakon,ITG *nslavnode,
		ITG *itiefac,char *tieset,
		ITG *islavnode,ITG *islavsurf,
		ITG **irowtlocp,ITG *jqtloc,double **autlocp,
		ITG **irowtlocinvp,ITG *jqtlocinv,double **autlocinvp,
		ITG *iflagdualquad);
    
void FORTRAN(checkspcmpc,(ITG *ntie,char *tieset,ITG *islavnode,ITG *imastnode,
			  ITG *nslavnode,ITG *nmastnode,double *slavnor,
			  ITG *islavact,ITG *nboun,ITG *ndirboun,double *xboun,
			  ITG *nodempc,double *coefmpc,ITG *ikboun,ITG *ilboun,
			  ITG *nmpc2,ITG *ipompc2,ITG *nodempc2,ITG *nslavspc,
			  ITG *islavspc,ITG *nsspc,ITG *nslavmpc,ITG *islavmpc,
			  ITG *nsmpc,ITG *nmspc,ITG *nmastmpc,ITG *imastmpc,
			  ITG *nmmpc));
    
void contactmortar(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,char *set,
		   ITG *istartset,ITG *iendset,ITG *ialset,ITG *itietri,
		   char *lakon,ITG *ipkon,ITG *kon,ITG *koncont,ITG *ne,
		   double *cg,double *straight,double *co,
		   double *vold,ITG *ielmat, double *elcon,
		   ITG *istep,ITG *iinc,ITG *iit,ITG *ncmat_,ITG *ntmat_,
		   ITG *ne0,double *vini,
		   ITG *nmethod,ITG *neq,ITG *nzs,ITG *nactdof,ITG *itiefac,
		   ITG *islavsurf,ITG *islavnode,ITG *imastnode,
		   ITG *nslavnode,ITG *nmastnode,double *ad,
		   double **aup,double *b,ITG **irowp,ITG *icol,ITG *jq,
		   ITG *imastop,
		   ITG *iponoels,ITG *inoels,ITG *nzsc,double **aucp,
		   double *adc,ITG **irowcp,ITG *jqc,ITG *islavact,
		   double *gap,
		   double *slavnor,double *slavtan,
		   double *bhat,
		   ITG **irowbdp,ITG *jqbd,double **aubdp,
		   ITG **irowbdtilp,ITG *jqbdtil ,double **aubdtilp,
		   ITG **irowbdtil2p,ITG *jqbdtil2,double **aubdtil2p,
		   ITG **irowddp,ITG *jqdd,double **auddp,
		   ITG **irowddtilp,ITG *jqddtil,double **auddtilp,
		   ITG **irowddtil2p,ITG *jqddtil2,double **auddtil2p,
		   ITG **irowddinvp,ITG *jqddinv,double **auddinvp,
		   ITG *irowtloc,ITG *jqtloc,double *autloc,  
		   ITG *irowtlocinv,ITG *jqtlocinv,double *autlocinv,   
		   ITG *mi,ITG *ipe,ITG *ime,double *tietol,ITG *iflagact,
		   double *cstress,
		   double *cstressini,double *bp_old,ITG *iflag_fric,ITG *nk,
		   ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		   ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		   ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		   ITG *nboun2,ITG *ndirboun2,ITG *nodeboun2,double *xboun2,
		   ITG *nmpc2,ITG *ipompc2,ITG *nodempc2,double *coefmpc2,
		   ITG *ikboun2,ITG *ilboun2,ITG *ikmpc2,ITG *ilmpc2,
		   ITG *nslavspc,ITG *islavspc,ITG *nsspc,ITG *nslavmpc,
		   ITG *islavmpc,ITG *nsmpc,
		   ITG *nslavspc2,ITG *islavspc2,ITG *nsspc2,ITG *nslavmpc2,
		   ITG *islavmpc2,ITG *nsmpc2,
		   ITG *nmastspc,ITG *imastspc,ITG *nmspc,ITG *nmastmpc,
		   ITG *imastmpc,ITG *nmmpc,
		   ITG *nmastmpc2,ITG *imastmpc2,ITG *nmmpc2,
		   double *pslavdual,double* pslavdualpg,
		   ITG *islavactdof,
		   ITG *islavactdoftie,
		   double *plicon,ITG *nplicon,ITG *npmat_,ITG *nelcon,
		   double *dtime,
		   ITG *islavnodeinv,
		   double **Bdp,ITG **irowbp,ITG *jqb,
		   double **Bdhelpp,ITG **irowbhelpp,ITG *jqbhelp,
		   double **Ddp,ITG **irowdp,ITG *jqd,
		   double **Ddtilp,ITG **irowdtilp,ITG *jqdtil,
		   double **Bdtilp,ITG **irowbtilp,ITG *jqbtil,
		   double **Bpgdp,ITG **irowbpgp,ITG *jqbpg,
		   double **Dpgdp,ITG **irowdpgp,ITG *jqdpg,
		   double **Dpgdtilp,ITG **irowdpgtilp,ITG *jqdpgtil,
		   double **Bpgdtilp,ITG **irowbpgtilp,ITG *jqbpgtil,
		   double *lambdaiwan,double *lambdaiwanini,double *bet,
		   ITG *iflagdualquad,
		   char *labmpc2,double *cfsinitil,double *reltime,
		   ITG *ithermal,double *plkcon,ITG *nplkcon);

void stressmortar(double *bhat,double *adc,double *auc,ITG *jqc,
		  ITG *irowc,ITG *neq,double *gap,double *b,ITG *islavact,
		  ITG *irowddinv,ITG *jqddinv,double *auddinv,
		  ITG *irowtloc,ITG *jqtloc,double *autloc, 
		  ITG *irowtlocinv,ITG *jqtlocinv,double *autlocinv,
		  ITG *ntie,ITG *nslavnode,
		  ITG *islavnode,ITG *nmastnode,ITG *imastnode,double *slavnor,double *slavtan,
		  ITG *nactdof,ITG *iflagact,double *cstress,double *cstressini,ITG *mi,
		  double *cdisp,double *f_cs,double *f_cm,ITG *iit,ITG *iinc,
		  double *vold,double *vini,double* bp,ITG *nk,
		  ITG *nboun2,ITG *ndirboun2,ITG *nodeboun2,double *xboun2,
		  ITG *nmpc2,ITG *ipompc2,ITG *nodempc2,double *coefmpc2,
		  ITG *ikboun2,ITG *ilboun2,ITG *ikmpc2,ITG *ilmpc2,
		  ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		  ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		  ITG *nslavspc2,ITG *islavspc2,ITG *nsspc2,ITG *nslavmpc2,ITG *islavmpc2,ITG *nsmpc2,
		  ITG *nmastspc,ITG *imastspc,ITG *nmspc,ITG *nmastmpc,ITG *imastmpc,ITG *nmmpc,
		  char *tieset,
		  double  *elcon,double *tietol,ITG *ncmat_,ITG *ntmat_,
		  double *plicon,ITG *nplicon,ITG *npmat_,ITG *nelcon,
		  double *dtime,double *cfs,double *cfm,ITG *islavnodeinv,
		  double *Bd,ITG *irowb,ITG *jqb,
		  double *Dd,ITG *irowd,ITG *jqd,
		  double *Ddtil,ITG *irowdtil,ITG *jqdtil,
		  double *Bdtil,ITG *irowbtil,ITG *jqbtil,
		  double *Bpgd,ITG *irowbpg,ITG *jqbpg,
		  double *Dpgd,ITG *irowdpg,ITG *jqdpg,
		  double *lambdaiwan,double *lambdaiwanini,ITG *nmethod,double *bet,
		  ITG *iflagdualquad,ITG *ithermal,ITG *iperturb,
		  char *labmpc,char *labmpc2,double *cam,double *veold,double *accold,double *gam,ITG *nk2,
		  double *cfsini,double *cfstil,double *plkcon,ITG *nplkcon,char *filab,double *f,double *fn,
		  double *qa,ITG *nprint,char *prlab,double *xforc,
		  ITG *nforc);
    
void FORTRAN(catsmpcslavno,(ITG *ntie,ITG *islavnode,ITG *imastnode,
			    ITG *nslavnode,ITG *nmastnode,ITG *nboun,
			    ITG *ndirboun,ITG *nodeboun,ITG *nmpc,ITG *ipompc,
			    ITG *nodempc,ITG *ikbou,ITG *ilbou,ITG *ikmpc,
			    ITG *ilmpc,ITG *nboun2,ITG *nmpc2,ITG *ipompc2,
			    ITG *nodempc2,ITG *ikbou2,ITG *ilbou2,ITG *ikmpc2,
			    ITG *ilmpc2,ITG *nslavspc,ITG *islavspc,ITG *nsspc,
			    ITG *nslavmpc,ITG *islavmpc,ITG *nsmpc,
			    ITG *nslavspc2,ITG *islavspc2,ITG *nsspc2,
			    ITG *nslavmpc2,ITG *islavmpc2,ITG *nsmpc2,
			    ITG *nmastspc,ITG *imastspc,ITG *nmspc,
			    ITG *nmastmpc,ITG *imastmpc,ITG *nmmpc,
			    ITG *nmastmpc2,ITG *imastmpc2,ITG *nmmpc2));      
       
void FORTRAN(createbd,(ITG *ict,ITG *l,ITG *ipkon,ITG *kon,char *lakon,
		       double *co,double *vold,double* gapmints,ITG *islavsurf,
		       ITG *imastsurf,double *pmastsurf,double *contr,
		       ITG *isconstr,ITG *imcontr,double *dcontr,ITG *idcontr1,
		       ITG *idcontr2,double *gcontr,ITG *igcontr,ITG *mi,
		       double* pslavsurf,double* pslavdual,ITG *nslavnode,
		       ITG *islavnode,ITG *nmastnode,ITG *imastnode,
		       ITG *icounter,ITG *icounter2,ITG *islavact,
		       ITG *iflagdualquad));
     
void FORTRAN(createbd_pg,(ITG *ict,ITG *l,ITG *ipkon,ITG *kon,char *lakon,double *co,double *vold,
			  double* gapmints,ITG *islavsurf,ITG *imastsurf,
			  double *pmastsurf,double *contr,ITG *isconstr,ITG *imcontr,
			  double *dcontr,ITG *idcontr1,ITG *idcontr2,double *gcontr,ITG *igcontr,
			  ITG *iponoels,ITG *inoels,ITG *mi,double* pslavsurf,
			  double* pslavdual,ITG *nslavnode,ITG *islavnode,ITG *nmastnode,ITG *imastnode,
			  ITG *icounter,ITG *icounter2,ITG *islavact,ITG *iflagdualquad));
     
void FORTRAN(createtele,(ITG *ipkon,ITG *kon,char *lakon,ITG *islavsurf,
			 double *dcontr,ITG *idcontr1,ITG *idcontr2,
			 ITG *icounter,ITG *l));

void FORTRAN(createteleinv,(ITG *ipkon,ITG *kon,char *lakon,ITG *islavsurf,
			    double *dcontr,ITG *idcontr1,ITG *idcontr2,
			    ITG *icounter,ITG *l));
      
void FORTRAN(createtele_lin,(ITG *ipkon,ITG *kon,char *lakon,ITG *islavsurf,
			     double *dcontr,ITG *idcontr1,ITG *idcontr2,
			     ITG *icounter,ITG *l));

void FORTRAN(createteleinv_lin,(ITG *ipkon,ITG *kon,char *lakon,ITG *islavsurf,
				double *dcontr,ITG *idcontr1,ITG *idcontr2,
				ITG *icounter,ITG *l));

void decascade_mortar(ITG *nmpc,ITG *ipompc,ITG **nodempcp,double **coefmpcp,
		      ITG *ikmpc,ITG *ilmpc,ITG *memmpc_,ITG *mpcfree);

void FORTRAN(nortanslav,(char *tieset,ITG *ntie,ITG *ipkon,ITG *kon,
			 char *lakon,char *set,double *co,double *vold,
			 ITG *nset,ITG *islavsurf,ITG *itiefac,ITG *islavnode,
			 ITG *nslavnode,double *slavnor,double *slavtan,
			 ITG *mi));

void FORTRAN(gendualcoeffs,(char *tieset,ITG *ntie,ITG *ipkon,ITG *kon,
			    char *lakon,double *co,double *vold,ITG *islavact,
			    ITG *islavsurf,ITG *itiefac,ITG *islavnode,
			    ITG *nslavnode,ITG *mi,double *pslavsurf,
			    double* pslavdual,double* pslavdualpg,
			    ITG *iflagdualquad));

void FORTRAN(genfirstactif,(char *tieset,ITG *ntie,ITG *itietri,ITG *ipkon,
			    ITG *kon,char *lakon,double *cg,double *straight,
			    double *co,double *vold,double *xo,double *yo,
			    double *zo,double *x,double *y,double *z,ITG *nx,
			    ITG *ny,ITG *nz,ITG *istep,ITG *iinc,ITG *iit,
			    ITG *mi,ITG *imastop,
			    ITG *nslavnode,ITG *islavnode,ITG *islavsurf,
			    ITG *itiefac,double *areaslav,char *set,ITG *nset,
			    ITG *istartset,ITG *iendset,ITG *ialset,
			    ITG *islavact,ITG *ifree,double *tietol));

void FORTRAN(genislavactdof,(ITG *ntie,char *tieset,ITG *nactdof,
			     ITG *nslavnode,ITG *nmastnode,ITG *imastnode,
			     ITG *islavactdof,ITG *islavnode,ITG *mi,
			     ITG *ithermal));
     
void FORTRAN(genislavelinv,(ITG *islavelinv,ITG *jqtloc,char *lakon,
			    ITG *ipkon,ITG *kon,ITG *ne,ITG *nasym));

void FORTRAN(getcontactparams,(double *mu,ITG *regmode,ITG *regmodet,double *fkninv,double *fktauinv,
			       double *p0,double *beta,double *tietol,double *elcon,ITG *itie,
			       ITG *ncmat_,ITG *ntmat_,ITG *niwan));
     
ITG FORTRAN(getlocno,(ITG *m,ITG *jfaces,ITG *nope));

void FORTRAN(getnumberofnodes,(ITG *nelems,ITG *jfaces,char *lakon,ITG *nope,
			       ITG *nopes,ITG *idummy)); 

void inimortar(double **enerp,ITG *mi,ITG *ne ,ITG *nslavs,ITG *nk,ITG *nener,
	       ITG **ipkonp,char **lakonp,ITG **konp,ITG *nkon,
	       ITG *maxprevcontel,double **xstatep,ITG *nstate_,
	       ITG **islavactdoftiep,double **bpp,ITG **islavactp,
	       double **gapp,double **slavnorp,double **slavtanp,double **cdispp,
	       double **cstressp,double **cfsp,double **cfmp,double **cfsinip,double **cfsinitilp,double **cfstilp,
	       double **bpinip,ITG **islavactinip,double **cstressinip,
	       ITG *niwan,ITG *ntie,char *tieset,
	       ITG *nslavnode,ITG *islavnode,
	       ITG **islavnodeinvp,ITG **islavelinvp,double **pslavdualp,double **pslavdualpgp,
	       double **autlocp,ITG **irowtlocp,ITG **jqtlocp,
	       double **autlocinvp,ITG **irowtlocinvp,ITG **jqtlocinvp,
	       double **Bdp,ITG **irowbp,ITG **jqbp,
	       double **Bdhelpp,ITG **irowbhelpp,ITG **jqbhelpp,
	       double **Ddp,ITG **irowdp,ITG **jqdp,
	       double **Ddtilp,ITG **irowdtilp,ITG **jqdtilp,
	       double **Bdtilp,ITG **irowbtilp,ITG **jqbtilp,
	       double **Bpgdp,ITG **irowbpgp,ITG **jqbpgp,
	       double **Dpgdp,ITG **irowdpgp,ITG **jqdpgp,
	       double **Dpgdtilp,ITG **irowdpgtilp,ITG **jqdpgtilp,
	       double **Bpgdtilp,ITG **irowbpgtilp,ITG **jqbpgtilp,
	       ITG *iflagdualquad,ITG *itiefac,ITG *islavsurf,
	       ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
	       ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
	       ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
	       ITG *nboun2,ITG **ndirboun2p,ITG **nodeboun2p,double **xboun2p,
	       ITG *nmpc2,ITG **ipompc2p,ITG **nodempc2p,double **coefmpc2p,char **labmpc2p,
	       ITG **ikboun2p,ITG **ilboun2p,ITG **ikmpc2p,ITG **ilmpc2p,
	       ITG **nslavspcp,ITG **islavspcp,ITG **nslavmpcp,ITG **islavmpcp,
	       ITG **nslavspc2p,ITG **islavspc2p,ITG **nslavmpc2p,ITG **islavmpc2p,
	       ITG **nmastspcp,ITG **imastspcp,ITG **nmastmpcp,ITG **imastmpcp,
	       ITG **nmastmpc2p,ITG **imastmpc2p,ITG *nmmpc2,
	       ITG *nsspc,ITG *nsspc2,ITG *nsmpc,ITG *nsmpc2,
	       ITG *imastnode,ITG *nmastnode,ITG *nmspc,ITG *nmmpc,
	       ITG *iponoels,ITG *inoels,
	       double *tietol,double *elcon,ITG *ncmat_,ITG *ntmat_,ITG *nasym,
	       ITG *iflag_fric,double **lambdaiwanp,double **lambdaiwaninip,
	       ITG *nk2,double *vold,ITG *nset,char *set,ITG *mortar,ITG *memmpc_,
	       ITG **ielmatp,ITG **ielorienp,ITG *norien,ITG *nmethod,
	       ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	       ITG **nodeforc2p,ITG **ndirforc2p,double **xforc2p,ITG *nforc2);

void insertas(ITG **irowp,ITG **mast1p,ITG *i1,
	      ITG *i2,ITG *ifree,ITG *nzs_,double *contribution,double **bdp);

void insertas_ws(ITG **irowp, ITG *i1,ITG *i2, ITG *ifree, ITG *nzs_,
		 double *contribution, double **bdp);

void matrixsort(double *au,ITG *mast1,ITG *irow,ITG *jq,
		ITG *nzs,ITG *dim);
      
void mortar_prefrd(ITG *ne,ITG *nslavs,ITG *mi,ITG *nk,ITG *nkon,
		   double **stxp,double *cdisp,
		   double *fn,double *cfs,double *cfm);
   
void mortar_postfrd(ITG *ne,ITG *nslavs,ITG *mi,ITG *nk,ITG *nkon,
		    double *fn,double *cfs,double *cfm);
       
void multimortar(double **aup,double *ad,ITG **irowp,ITG *jq,ITG *nzs,
		 double **aucp,double *adc,ITG **irowcp,ITG *jqc,ITG *nzsc,
		 double *aubd,ITG *irowbd,ITG *jqbd,
		 double *aubdtil,ITG *irowbdtil,ITG *jqbdtil,
		 double *aubdtil2,ITG *irowbdtil2,ITG *jqbdtil2,
		 ITG *irowdd,ITG *jqdd,double *audd,
		 ITG *irowddtil2,ITG *jqddtil2,double *auddtil2,
		 ITG *irowddinv,ITG *jqddinv,double *auddinv,
		 double *Bd,ITG *irowb,ITG *jqb,
		 double *Dd,ITG *irowd,ITG *jqd,
		 double *Ddtil,ITG *irowdtil,ITG *jqdtil,
		 ITG *neq,double *b,double *bhat,ITG *islavnode,ITG *imastnode,
		 ITG *nslavnode,ITG *nmastnode,
		 ITG *islavact,ITG *islavactdof,
		 double *gap,
		 double *slavnor,double *slavtan,
		 double *vold,double *vini,double *cstress,double *cstressini,
		 double *bp_old,ITG *nactdof,ITG *ntie,ITG *mi,ITG *nk,
		 ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		 ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		 ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		 ITG *nslavspc,ITG *islavspc,ITG *nsspc,ITG *nslavmpc,ITG *islavmpc,ITG *nsmpc,
		 ITG *nmastspc,ITG *imastspc,ITG *nmspc,ITG *nmastmpc,ITG *imastmpc,ITG *nmmpc,
		 char *tieset,
		 ITG *islavactdoftie,ITG *nelcon,double  *elcon,double *tietol,ITG *ncmat_,ITG *ntmat_,
		 double *plicon,ITG *nplicon,ITG *npmat_,double *dtime,
		 ITG *irowtloc,ITG *jqtloc,double *autloc,
		 ITG *irowtlocinv,ITG *jqtlocinv,double *autlocinv,
		 ITG *islavnodeinv,double *lambdaiwan,double *lambdaiwanini,ITG *iit,ITG *nmethod,double *bet,ITG *ithermal,
		 double *plkcon,ITG *nplkcon
		 );

void multi_rect(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1,ITG m_1,
		double *au_2,ITG * irow_2,ITG * jq_2,ITG n_2,ITG m_2,
		double **au_rp,ITG **irow_rp,ITG * jq_r,ITG *nzs);

void multi_rectv(double *au_1,ITG * irow_1,ITG * jq_1,ITG n_1,ITG m_1,
		 double * b,double ** v_rp);

void multi_scal(double *au_1,ITG * irow_1,ITG * jq_1,
		double *au_2,ITG * irow_2,ITG * jq_2,
		ITG m,ITG n,double*value,ITG *flag);
       
void FORTRAN(opnonsym,(ITG *neq,double *aux,double *b,double *bhat,
		       double *bdd,double*bdu,ITG *jqbd,ITG *irowbd));

void FORTRAN(opnonsymt,(ITG *neq,double *aux,double *b,double *bhat,
			double *bdd,double*bdu,ITG *jqbd,ITG *irowbd));

void premortar(ITG *iflagact,ITG *ismallsliding,ITG *nzs,ITG *nzsc2,
	       double **auc2p,double **adc2p,ITG **irowc2p,ITG **icolc2p,
	       ITG **jqc2p,
	       double **aubdp,ITG **irowbdp,ITG **jqbdp,
	       double **aubdtilp,ITG **irowbdtilp,ITG **jqbdtilp,
	       double **aubdtil2p,ITG **irowbdtil2p,ITG **jqbdtil2p,
	       double **auddp,ITG **irowddp,ITG **jqddp,
	       double **auddtilp,ITG **irowddtilp,ITG **jqddtilp,
	       double **auddtil2p,ITG **irowddtil2p,ITG **jqddtil2p,
	       double **auddinvp,ITG **irowddinvp,ITG **jqddinvp,
	       ITG **jqtempp,ITG **irowtempp,ITG **icoltempp,ITG *nzstemp,
	       ITG *iit,double *slavnor,double *slavtan,
	       ITG *icol,ITG *irow,ITG *jq,
	       ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
	       ITG *nboun2,ITG **ndirboun2p,ITG **nodeboun2p,
	       double **xboun2p,
	       ITG *nmpc2,ITG **ipompc2p,ITG **nodempc2p,double **coefmpc2p,
	       char **labmpc2p,
	       ITG **ikboun2p,ITG **ilboun2p,ITG **ikmpc2p,ITG **ilmpc2p,
	       ITG **nslavspcp,ITG **islavspcp,ITG **nslavmpcp,
	       ITG **islavmpcp,
	       ITG **nslavspc2p,ITG **islavspc2p,ITG **nslavmpc2p,
	       ITG **islavmpc2p,
	       ITG **nmastspcp,ITG **imastspcp,ITG **nmastmpcp,
	       ITG **imastmpcp,
	       ITG **nmastmpc2p,ITG **imastmpc2p,ITG *nmmpc2,
	       ITG *nsspc,ITG *nsspc2,ITG *nsmpc,ITG *nsmpc2,
	       ITG *imastnode,ITG *nmastnode,ITG *nmspc,ITG *nmmpc,
	       double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,double *stn,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,double *prestr,
	       ITG *iprestr,char *filab,double *eme,double *emn,
	       double *een,ITG *iperturb,double *f,ITG *nactdof,
	       ITG *iout,double *qa,
	       double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
	       double *xbounact,double *xboun,ITG *nboun,ITG *ipompc,
	       ITG *nodempc,
	       double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
	       ITG *neq,double *veold,double *accold,
	       double *dtime,double *time,
	       double *ttime,double *plicon,
	       ITG *nplicon,double *plkcon,ITG *nplkcon,
	       double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
	       char *matname,ITG *mi,ITG *ielas,
	       ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
	       double *vini,double *ener,
	       double *enern,double *emeini,double *xstaten,double *eei,
	       double *enerini,double *cocon,ITG *ncocon,char *set,
	       ITG *nset,ITG *istartset,
	       ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
	       char *prset,double *qfx,double *qfn,double *trab,
	       ITG *inotr,ITG *ntrans,ITG *nelemload,
	       ITG *nload,ITG *istep,ITG *iinc,
	       double *springarea,double *reltime,ITG *ne0,double *xforc,
	       ITG *nforc,double *thicke,
	       double *shcon,ITG *nshcon,char *sideload,double *xload,
	       double *xloadold,ITG *icfd,ITG *inomat,
	       ITG *islavelinv,ITG *islavsurf,
	       ITG *iponoels,ITG *inoels,
	       ITG *mortar,ITG *nslavnode,ITG *islavnode,ITG *nslavs,
	       ITG *ntie,
	       double *autloc,ITG *irowtloc,ITG *jqtloc,
	       double *autlocinv,ITG *irowtlocinv,ITG *jqtlocinv,
	       ITG *nk2,ITG *iflagdualquad,
	       char *tieset,ITG *itiefac  ,ITG *rhsi,
	       double *au,double *ad,double **f_cmp,double **f_csp,
	       double *t1act,double *cam,double *bet,double *gam,
	       double *epn,
	       double *xloadact,ITG *nodeforc,ITG *ndirforc,double *xforcact,
	       double *xbodyact,ITG *ipobody,ITG *nbody,double *cgr,
	       ITG *nzl,double *sti,ITG *iexpl,ITG *mass,ITG *buckling,
	       ITG *stiffness,
	       ITG *intscheme,double *physcon,ITG *coriolis,ITG *ibody,
	       ITG *integerglob,double *doubleglob,ITG *nasym,
	       double *alpham,double *betam,double *auxtil2,
	       double *pslavsurf,double *pmastsurf,
	       double *clearini,ITG *ielprop,double *prop,
	       ITG *islavact,double *cdn,ITG *memmpc_,
	       double *cvinitil,double *cvtil,ITG *idamping,
	       ITG *ilin,ITG *iperturb_sav,double *adb,double *aub,
	       ITG **nodeforc2p,ITG **ndirforc2p,double **xforc2p,
	       ITG *nforc2,
	       ITG *itietri,double *cg,double *straight,ITG *koncont,
	       double *energyini,
	       double *energy,ITG *kscale,ITG *iponoel,ITG *inoel,ITG *nener,
	       char *orname,ITG *network,
	       char *typeboun,ITG *num_cpus,double *t0g,double *t1g,
	       double *smscale,ITG *mscalmethod);
       
void FORTRAN(regularization_gn_c,(double *lambdap,ITG *divmode,ITG *regmode,
				  double *gnc,double *aninvloc,double *p0,
				  double *beta,double *elcon,ITG *nelcon,ITG *itie,ITG *ntmat_,
				  double *plicon,ITG *nplicon,ITG *npmat_,ITG *ncmat_,
				  double *tietol,double *scal));
     
void FORTRAN(regularization_gt_c,(double *lambdatt,ITG *divmode,ITG *regmode,
				  double *gtc,double *atauinvloc));
     
void FORTRAN(regularization_slip_iwan,(double *lambdan,double *ut,
				       double *bp,double *atau2,double *resreg,
				       ITG *divmode,ITG *regmode,double *lambdaiwan,double *lambdaiwanini,
				       ITG *inode,double *n,double *t,double *mu,double *rslip,double *ltslip,
				       double *ltu,ITG *yielded,ITG *iit,ITG *debug,ITG *niwan,double *dut));
 
void FORTRAN(regularization_slip_lin,(double *utilt,double *bp,double *atauinv,
				      double *resreg,ITG *divmode,
				      ITG *islavact,double *lambdat,
				      double *lambdatilt,double *constantt,
				      ITG *debug,ITG *inode,double *n2,
				      double *t,double *that,double *mu,
				      double *rslip,double *ltslip,
				      double *ltu));
     
void FORTRAN(resultsini_mortar,(int *nk,double *v,int *ithermal,
				int *iperturb,int *nactdof,int *iout,
				double *vold,double *b,int *nodeboun,int *ndirboun,
				double *xboun,int *nboun,int *ipompc,int *nodempc,double *coefmpc,
				char *labmpc,int *nmpc,int *nmethod,double *cam,
				double *bet,double *gam,double *dtime,
				int *mi));     

void FORTRAN(slavintmortar,(ITG *ntie,ITG *itietri,ITG *ipkon,ITG *kon,
			    char *lakon,double *straight,ITG *nintpoint,
			    ITG *koncont,double *co,double *vold,double *xo,
			    double *yo,double *zo,double *x,double *y,
			    double *z,ITG *nx,ITG *ny,ITG *nz,ITG *iinc,
			    ITG *islavsurf,ITG *imastsurf,double *pmastsurf,
			    ITG *islavnode,ITG *nslavnode,ITG *imastop,
			    double *gap,ITG *islavact,ITG *mi,ITG *ncont,
			    ITG *ipe,ITG *ime,double *pslavsurf,ITG *i,ITG *l,
			    ITG *ntri,double *tietol,double *reltime,
			    ITG *nmethod));

void transformspcsmpcs_quad(ITG *nboun,ITG *ndirboun,ITG *nodeboun,
			    double *xboun,ITG *nmpc,ITG *ipompc,ITG *nodempc,
			    double *coefmpc,char *labmpc,ITG *ikboun,
			    ITG *ilboun,ITG *ikmpc,ITG *ilmpc,ITG *nboun2,
			    ITG **ndirboun2p,ITG **nodeboun2p,double **xboun2p,
			    ITG *nmpc2,ITG **ipompc2p,ITG **nodempc2p,
			    double **coefmpc2p,char **labmpc2p,ITG **ikboun2p,
			    ITG **ilboun2p,ITG **ikmpc2p,ITG **ilmpc2p,
			    ITG *irowtlocinv,ITG *jqtlocinv,double *autlocinv,
			    ITG *nk,ITG *nk2,ITG *iflagdualquad,ITG *ntie,
			    char *tieset,ITG *itiefac,ITG *islavsurf,
			    char *lakon,ITG *ipkon,ITG *kon,ITG *mi,
			    ITG *memmpc_,ITG *nodeforc,ITG *ndirforc,
			    double *xforc,ITG *nforc,ITG **nodeforc2p,
			    ITG **ndirforc2p,double **xforc2p,ITG *nforc2);
    
void transpose(double *au,ITG *jq,ITG *irow,ITG *dim,
	       double *au_t,ITG *jq_t,ITG *irow_t);

void trafontmortar2(ITG *neq,ITG *nzs,ITG *islavactdof,ITG *islavact,ITG *nslavnode,ITG *nmastnode,
		    double *f_da,double *f_atil,
		    double *au_dan,ITG *irow_dan,ITG *jq_dan,
		    double *au_dam,ITG *irow_dam,ITG *jq_dam,
		    double *au_dai,ITG *irow_dai,ITG *jq_dai,
		    double *au_daa,ITG *irow_daa,ITG *jq_daa,
		    double **au_antilp,ITG **irow_antilp,ITG *jq_antil,
		    double **au_amtilp,ITG **irow_amtilp,ITG *jq_amtil,
		    double **au_aitilp,ITG **irow_aitilp,ITG *jq_aitil,
		    double **au_aatilp,ITG **irow_aatilp,ITG *jq_aatil,
		    double *gap,
		    double *Bd,ITG *irowb,ITG *jqb,
		    double *Dd,ITG *irowd,ITG *jqd,
		    double *Ddtil,ITG *irowdtil,ITG *jqdtil,
		    double *au_bdtil2,ITG *irow_bdtil2,ITG *jq_bdtil2,
		    double *au_ddtil2i,ITG *irow_ddtil2i,ITG *jq_ddtil2i,
		    double *au_ddtil2a,ITG *irow_ddtil2a,ITG *jq_ddtil2a,
		    ITG *m_flagr,ITG *i_flagr,ITG *a_flagr,ITG *a_flag,ITG *i_flag,ITG *m_flag,
		    ITG *row_ln,ITG *row_lm,ITG *row_li,ITG *row_la,
		    double *slavnor,double *slavtan,
		    double *vold,double *vini,double *cstress,double *cstressini,
		    double *bp_old,ITG *nactdof,ITG *islavnode,ITG *imastnode,ITG *ntie,ITG *mi,ITG *nk,
		    ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		    ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		    ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		    ITG *nslavspc,ITG *islavspc,ITG *nsspc,ITG *nslavmpc,ITG *islavmpc,ITG *nsmpc,
		    ITG *nmastspc,ITG *imastspc,ITG *nmspc,ITG *nmastmpc,ITG *imastmpc,ITG *nmmpc,
		    char *tieset,
		    ITG *islavactdoftie,ITG *nelcon,double  *elcon,double *tietol,ITG *ncmat_,ITG *ntmat_,
		    double *plicon,ITG *nplicon,ITG *npmat_,double *dtime,
		    ITG *irowtloc,ITG *jqtloc,double *autloc,
		    ITG *irowtlocinv,ITG *jqtlocinv,double *autlocinv,
		    ITG *islavnodeinv,double *lambdaiwan,double *lambdaiwanini,ITG *iit,ITG *nmethod,double *beta,ITG *ithermal,
		    double *plkcon,ITG *nplkcon);
    
void trafontspcmpc( double *n,double *t,double *n2,double *that,ITG *islavnodeentry,
		    ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		    ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		    ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		    ITG *nslavspc,ITG *islavspc,ITG *nsspc,ITG *nslavmpc,ITG *islavmpc,ITG *nsmpc,
		    ITG *nmastspc,ITG *imastspc,ITG *nmspc,ITG *nmastmpc,ITG *imastmpc,ITG *nmmpc,
		    ITG *debug,ITG *node);

void FORTRAN(updatecont,(ITG *koncont,ITG *ncont,double *co,double *vold,
			 double *cg,double *straight,ITG *mi));
 
void FORTRAN(writematrix,(double *au,double *ad,ITG *irow,ITG *jq,ITG *neq,
			  ITG *number));
 
void FORTRAN(writematrix2,(double *au,double *ad,ITG *irow,ITG *jq,ITG *neq,
			   ITG *number));

void FORTRAN(writevector,(double *ad,ITG *neq,ITG *number));