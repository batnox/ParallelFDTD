#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <malloc.h>
#include <sys/time.h>
#include <omp.h>
#define Min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define Max(X, Y) (((X) > (Y)) ? (X) : (Y))

int nproc,ierr;
int tag=1,llim[101],ulim[101],rank=0;

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
float **ERMS;
float **erms2;
float **exrms;
float **exrms2;
float **eyrms;
float **eyrms2;
float **erms;
float **den;
float **dent;
float **exi;
float **eyi;
float **exi1;
float **eyi1;
float **exs;
float **eys;
float **hzi;
float **hzs;
float **vx;
float **vy;
float **ext;
float **eyt;
float **exs1;
float **exs2;
float **eys1;
float **eys2;
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
