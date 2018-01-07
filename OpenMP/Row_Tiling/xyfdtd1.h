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
__declspec(align(64)) float ERMS[2500][2500];
__declspec(align(64)) float erms2[2500][2500];
__declspec(align(64)) float exrms[2500][2500];
__declspec(align(64)) float exrms2[2500][2500];
__declspec(align(64)) float eyrms[2500][2500];
__declspec(align(64)) float eyrms2[2500][2500];
__declspec(align(64)) float erms[2500][2500];
__declspec(align(64)) float den[2500][2500];
__declspec(align(64)) float dent[2500][2500];
__declspec(align(64)) float exi[2500][2500];
__declspec(align(64)) float eyi[2500][2500];
__declspec(align(64)) float exi1[2500][2500];
__declspec(align(64)) float eyi1[2500][2500];
__declspec(align(64)) float exs[2500][2500];
__declspec(align(64)) float eys[2500][2500];
__declspec(align(64)) float hzi[2500][2500];
__declspec(align(64)) float hzs[2500][2500];
__declspec(align(64)) float vx[2500][2500];
__declspec(align(64)) float vy[2500][2500];
__declspec(align(64)) float ext[2500][2500];
__declspec(align(64)) float eyt[2500][2500];
__declspec(align(64)) float exs1[2500][2500];
__declspec(align(64)) float exs2[2500][2500];
__declspec(align(64)) float eys1[2500][2500];
__declspec(align(64)) float eys2[2500][2500];
double t,dt,ds,dte,dtm,dteds,dtmds,dtsi,DTAC;
double c1,c2,c3,c4;
__declspec(align(64)) double xmid[2500];
__declspec(align(64)) double ymid[2500];
__declspec(align(64)) double sgdx0[2500];
__declspec(align(64)) double sgdy0[2500];
double E0,FREQ,OMEG;
double om2snu,QSM,crec,dnma,xfront;
double FNUM,RECOMB,EMOB,EDIF,ETEM,nu2omeg2;
double PARC,ACCEL,CORNUI,ndifmax;
__declspec(align(64)) double DINI[2500];
double TIMD,PRESSURE,TEMP0;
double DENG0,fnucor,amu,nu_omeg,e2_epsm;
__declspec(align(64)) double DIFFUSION[2500][2500];
__declspec(align(64)) double SIO[2500][2500];
__declspec(align(64)) double SIOT[2500][2500];
__declspec(align(64)) double pow1[2500][2500];
double powt;
__declspec(align(64)) double frqio[2500][2500];
__declspec(align(64)) double puis[2500][2500];
double fnuimax,cene;
__declspec(align(64)) double denp[2500][2500];
double difa,dife;
__declspec(align(64)) double frq[2500][2500];
__declspec(align(64)) int done[2500];
int IORDER,IABSOR,naccel;
int nmaxwell,nx,ny,nperdt,nx1,ny1,nec;
int isig,nlamb;
int ani;
int nstep,nini;
__declspec(align(64)) double imid[2500];
double inul;
__declspec(align(64)) double jmid[2500];
double inv_nperdt,inv_c,inv_cmasse, inv_ds2;
double AA1,BB1,AA2,BB2,AA3,BB3,EE3;
