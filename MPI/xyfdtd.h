#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <malloc.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>
#define Min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define Max(X, Y) (((X) > (Y)) ? (X) : (Y))

int nproc,ierr;

/*--------------------------------------------------------------------*/
int i,j,n;

const double c=2.99795e8,pi=3.14159265,eps0=8.854e-12,xmu0=1.55063706e-6,qe=1.602176487e-19,cmasse=9.10938215e-31,akb=1.3806503e-23;

double xd1,xd2,yd1,yd2,xxpr,timp,xxpr1;
double t,dt,ds,dte,dtm,dteds,dtmds,dtsi,DTAC;
double c1,c2,c3,c4;
double E0,FREQ,OMEG;
double om2snu,QSM,crec,dnma,xfront;
double FNUM,RECOMB,EMOB,EDIF,ETEM,nu2omeg2;
double PARC,ACCEL,CORNUI,ndifmax;
double TIMD,PRESSURE,TEMP0;
double DENG0,fnucor,amu,nu_omeg,e2_epsm;
double powt,fnuimax,cene;
double difa,dife;
int IORDER,IABSOR,naccel;
int nmaxwell,nx,ny,nperdt,nx1,ny1,nec;
int isig,nlamb;
int ani;
int nstep,nini,inul;
double inv_nperdt,inv_c,inv_cmasse, inv_ds2;
double AA1,BB1,AA2,BB2,AA3,BB3,EE3;
double **ERMS;
double **erms2;
double **exrms;
double **exrms2;
double **eyrms;
double **eyrms2;
double **erms;
double **den;
double **dent;
double **exi;
double **eyi;
double **exi1;
double **eyi1;
double **exs;
double **eys;
double **hzi;
double **hzs;
double **vx;
double **vy;
double **ext;
double **eyt;
double **exs1;
double **exs2;
double **eys1;
double **eys2;
double *xmid;
double *ymid;
double *sgdx0;
double *sgdy0;
double *DINI;
double **DIFFUSION;
double **SIO;
double **SIOT;
double **pow1;
double **frqio;
double **puis;
double **denp;
double **frq;
int *done;
int *imid;
int *jmid;
