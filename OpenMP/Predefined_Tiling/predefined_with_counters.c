#include "xyfdtd.h"

/*
***********************************************************
boht sare changes wala
***********************************************************
*/

int KELEC,KRMS,nmod,NECRIR,icpling,K;
    double output[421][2012];   //ulimp -> llimp:ulimp ??
    double slambx,slamby,tstop, radius = 0;
    int length,istart,iend,icc;
    time_t itime;   //long
    int start, end;
    float sine, sine1, x, x0;
    float sineb,sine1b,x0b;
    int i0, section_size, total_threads;
    double aa,betax,betay,alpha,gamma1,extk,eytk,const7,const8;
    double omp2x,omp2y,qmdt,const3,const4,const5x,const6x,const5y,const6y;
    long long add=0,mul=0,div_flo=0,sqrt_flo=0,pow_flo=0,trig_flo=0,exp_flo=0,assn=0,comp=0,reciprocal_flo=0,max_min_flo=0;

    FILE *fptr,*fptr1,*fptr2,*fptr3, *file_a, *file_xelec, *file_yelec, *file_eyi;

    void SETUP();
    void INC_EFIELD();
    void ADV_EFIELD_vel();
    void ADV_HFIELD();
    void RMS(int k);
    void MR_MUR(int row);
    // void MR_MUR();
    double FIONIZ(double EE,double PR);
    void ELEC_DENS();
    void ZERO();
    void anim();
    void free_all();
    void EFIELD_HFIELD();
    int main(){
    ERMS = malloc(sizeof(float *) * 6200);
if (ERMS){
	for (i = 0; i < 6200; i++){
		ERMS[i] = malloc(sizeof(float) * 6200);
	}
}
erms2 = malloc(sizeof(float *) * 6200);
if (erms2){
	for (i = 0; i < 6200; i++){
		erms2[i] = malloc(sizeof(float) * 6200);
	}
}
exrms = malloc(sizeof(float *) * 6200);
if (exrms){
	for (i = 0; i < 6200; i++){
		exrms[i] = malloc(sizeof(float) * 6200);
	}
}
exrms2 = malloc(sizeof(float *) * 6200);
if (exrms2){
	for (i = 0; i < 6200; i++){
		exrms2[i] = malloc(sizeof(float) * 6200);
	}
}
eyrms = malloc(sizeof(float *) * 6200);
if (eyrms){
	for (i = 0; i < 6200; i++){
		eyrms[i] = malloc(sizeof(float) * 6200);
	}
}
eyrms2 = malloc(sizeof(float *) * 6200);
if (eyrms2){
	for (i = 0; i < 6200; i++){
		eyrms2[i] = malloc(sizeof(float) * 6200);
	}
}
erms = malloc(sizeof(float *) * 6200);
if (erms){
	for (i = 0; i < 6200; i++){
		erms[i] = malloc(sizeof(float) * 6200);
	}
}
den = malloc(sizeof(float *) * 6200);
if (den){
	for (i = 0; i < 6200; i++){
		den[i] = malloc(sizeof(float) * 6200);
	}
}
dent = malloc(sizeof(float *) * 6200);
if (dent){
	for (i = 0; i < 6200; i++){
		dent[i] = malloc(sizeof(float) * 6200);
	}
}
exi = malloc(sizeof(float *) * 6200);
if (exi){
	for (i = 0; i < 6200; i++){
		exi[i] = malloc(sizeof(float) * 6200);
	}
}
eyi = malloc(sizeof(float *) * 6200);
if (eyi){
	for (i = 0; i < 6200; i++){
		eyi[i] = malloc(sizeof(float) * 6200);
	}
}
exi1 = malloc(sizeof(float *) * 6200);
if (exi1){
	for (i = 0; i < 6200; i++){
		exi1[i] = malloc(sizeof(float) * 6200);
	}
}
eyi1 = malloc(sizeof(float *) * 6200);
if (eyi1){
	for (i = 0; i < 6200; i++){
		eyi1[i] = malloc(sizeof(float) * 6200);
	}
}
exs = malloc(sizeof(float *) * 6200);
if (exs){
	for (i = 0; i < 6200; i++){
		exs[i] = malloc(sizeof(float) * 6200);
	}
}
eys = malloc(sizeof(float *) * 6200);
if (eys){
	for (i = 0; i < 6200; i++){
		eys[i] = malloc(sizeof(float) * 6200);
	}
}
hzi = malloc(sizeof(float *) * 6200);
if (hzi){
	for (i = 0; i < 6200; i++){
		hzi[i] = malloc(sizeof(float) * 6200);
	}
}
hzs = malloc(sizeof(float *) * 6200);
if (hzs){
	for (i = 0; i < 6200; i++){
		hzs[i] = malloc(sizeof(float) * 6200);
	}
}
vx = malloc(sizeof(float *) * 6200);
if (vx){
	for (i = 0; i < 6200; i++){
		vx[i] = malloc(sizeof(float) * 6200);
	}
}
vy = malloc(sizeof(float *) * 6200);
if (vy){
	for (i = 0; i < 6200; i++){
		vy[i] = malloc(sizeof(float) * 6200);
	}
}
ext = malloc(sizeof(float *) * 6200);
if (ext){
	for (i = 0; i < 6200; i++){
		ext[i] = malloc(sizeof(float) * 6200);
	}
}
eyt = malloc(sizeof(float *) * 6200);
if (eyt){
	for (i = 0; i < 6200; i++){
		eyt[i] = malloc(sizeof(float) * 6200);
	}
}
exs1 = malloc(sizeof(float *) * 6200);
if (exs1){
	for (i = 0; i < 6200; i++){
		exs1[i] = malloc(sizeof(float) * 6200);
	}
}
exs2 = malloc(sizeof(float *) * 6200);
if (exs2){
	for (i = 0; i < 6200; i++){
		exs2[i] = malloc(sizeof(float) * 6200);
	}
}
eys1 = malloc(sizeof(float *) * 6200);
if (eys1){
	for (i = 0; i < 6200; i++){
		eys1[i] = malloc(sizeof(float) * 6200);
	}
}
eys2 = malloc(sizeof(float *) * 6200);
if (eys2){
	for (i = 0; i < 6200; i++){
		eys2[i] = malloc(sizeof(float) * 6200);
	}
}
xmid = malloc(sizeof(double) * 6200);
ymid = malloc(sizeof(double) * 6200);
sgdx0 = malloc(sizeof(double) * 6200);
sgdy0 = malloc(sizeof(double) * 6200);
DINI = malloc(sizeof(double) * 6200);
DIFFUSION = malloc(sizeof(double *) * 6200);
if (DIFFUSION){
	for (i = 0; i < 6200; i++){
		DIFFUSION[i] = malloc(sizeof(double) * 6200);
	}
}
SIO = malloc(sizeof(double *) * 6200);
if (SIO){
	for (i = 0; i < 6200; i++){
		SIO[i] = malloc(sizeof(double) * 6200);
	}
}
SIOT = malloc(sizeof(double *) * 6200);
if (SIOT){
	for (i = 0; i < 6200; i++){
		SIOT[i] = malloc(sizeof(double) * 6200);
	}
}
pow1 = malloc(sizeof(double *) * 6200);
if (pow1){
	for (i = 0; i < 6200; i++){
		pow1[i] = malloc(sizeof(double) * 6200);
	}
}
frqio = malloc(sizeof(double *) * 6200);
if (frqio){
	for (i = 0; i < 6200; i++){
		frqio[i] = malloc(sizeof(double) * 6200);
	}
}
puis = malloc(sizeof(double *) * 6200);
if (puis){
	for (i = 0; i < 6200; i++){
		puis[i] = malloc(sizeof(double) * 6200);
	}
}
denp = malloc(sizeof(double *) * 6200);
if (denp){
	for (i = 0; i < 6200; i++){
		denp[i] = malloc(sizeof(double) * 6200);
	}
}
frq = malloc(sizeof(double *) * 6200);
if (frq){
	for (i = 0; i < 6200; i++){
		frq[i] = malloc(sizeof(double) * 6200);
	}
}
done = malloc(sizeof(int) * 6200);
imid = malloc(sizeof(int) * 6200);
jmid = malloc(sizeof(int) * 6200);
        struct timeval begin, end, total_start, total_end;
        double t1, t2, t_inc_efield, t_adv_efield_vel, t_mr_mur, t_hfield, t_efield_hfield,t_elec_dens,t_anim;
        gettimeofday(&begin, NULL);
        t1 = begin.tv_usec;
        gettimeofday(&total_start, NULL);
        nini=1;
        fptr=fopen("xstart_800.dat","r");

        char *line = NULL;
        size_t len = 0;
        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        tstop=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        nlamb=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        slambx=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        slamby=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        E0=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        PRESSURE=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        FREQ=atof(line);
        len=0;


        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        IABSOR=atof(line);
        len=0;

        if(nini<=0)
        nini=1;  
        
        K=1;
        for(K;K<=nini;K++){
        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        xmid[K]=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        ymid[K]=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        sgdx0[K]=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        sgdy0[K]=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        DINI[K]=atof(line);
        len=0; 
        }
          
        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        icpling=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        nmaxwell=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        ndifmax=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        naccel=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        NECRIR=atof(line);
        len=0; 
          
        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        crec=atof(line);
        len=0;

        getline(&line, &len, fptr);
        getline(&line, &len, fptr);
        len=10;
        getline(&line, &len, fptr);
        cene=atof(line);
        len=0; 

        fclose(fptr);

        gettimeofday(&end, NULL);
        printf("Reading Input time: %f s\n", ((end.tv_sec - begin.tv_sec) + 
        ((end.tv_usec - begin.tv_usec)/1000000.0)));


        //======================Caculating Constants=====================
           gettimeofday(&begin, NULL);
           t1 = begin.tv_usec;

           nx=nlamb*slambx;
           ny=nlamb*slamby;
           nperdt=2.0*nlamb;
           nx1=nx-1;
           ny1=ny-1;
           t=0;
           nec=0;
           ani=0;
           nstep=0;
        //!=============================================================

        //!============================================================
        ZERO();
        SETUP();

        total_threads=atoi(getenv("OMP_NUM_THREADS"));//omp_get_num_threads();
        printf("%s threads\n",getenv("OMP_NUM_THREADS"));
        // exit(0);
        done[total_threads]=1;
        printf("total threads %d\n", total_threads);
        section_size=nx/total_threads;

        printf("totallll and section %d %d\n",total_threads,section_size);

        gettimeofday(&end, NULL);

        printf("Initialization time: %f s\n", ((end.tv_sec - begin.tv_sec) + 
          ((end.tv_usec - begin.tv_usec)/1000000.0)));
        ELEC_DENS();

        //*----
        printf("computing up to time  %f \n",tstop);

        //!=============================================================
        fptr=fopen("output.txt","w");      

        //C--------------- Initial parameters ouput-----
        fprintf(fptr,"nx=%d   ny=%d\n",nx,ny);
        fprintf(fptr,"DS=%f   DT=%f\n",ds,dt);
        fprintf(fptr,"Freq=%f   Omega=%f\n",FREQ,OMEG);
        fprintf(fptr,"Time period=%f\n",1.0/FREQ);
        fprintf(fptr,"Lambda=%f\n",3.0e8/FREQ);
           
        //!put 0    
        fprintf(fptr,"Collision Freq=%f\n",FNUM);
        fprintf(fptr,"Recombination Coef=%f\n",RECOMB);
        fprintf(fptr,"Mobility=%f\n",EMOB);
        fprintf(fptr,"Electron Temp= (eV)%f\n",ETEM);
        fprintf(fptr,"DIFFUSION Coef=%f\n",EDIF);
        fprintf(fptr,"Initial Gas/Neutral density=%f\n",DENG0);
        fprintf(fptr,"Cutoff-density=%f\n",(eps0*cmasse/pow(qe,2))*(pow(OMEG,2)+pow(FNUM,2)));// check formula
        fprintf(fptr, "Tstop %.10lf\n",tstop); 
        fclose(fptr);
        
        //!=============================================================
        KELEC=0;
        // n=0;
        n=0;
        int f=0;

        fptr1=fopen("eyi1.dat","w");
        fptr2=fopen("eyt.dat","w");
        fptr3=fopen("hzs.dat","w");
        file_eyi=fopen("file_eyi.dat","w");
        //!=============================================================

       
        // INC_EFIELD();
        // ADV_EFIELD_vel();
        // MR_MUR();

       

        do{
        n=n+1;
        // printf("Inside while\n");
        if(n%2000==0)
            printf("%d \n", n);

        nstep += 1;
        //if(n>=500000)
           //break;
       
        if(TIMD>tstop){
            printf("%e %d %d\n", dt,nx,ny); 
            printf("OUTPUT 1\n");
            gettimeofday(&total_end, NULL);
            //total_end = end.tv_usec;
            printf("EFIELD_HFIELD time: %f s\n", t_efield_hfield);
            printf("ELEC_DENS time: %f s\n", t_elec_dens);
            printf("ANIM time: %f s\n",t_anim);
            printf("Total time: %f s\n", ((total_end.tv_sec - total_start.tv_sec) + 
                    ((total_end.tv_usec - total_start.tv_usec)/1000000.0)));
            printf("Flops : %lld\n",add+mul+assn+2*div_flo+reciprocal_flo+2*pow_flo+3*trig_flo+sqrt_flo+3*exp_flo+comp+3*max_min_flo); 
            free_all();
            exit(0);
        }
        // INC_EFIELD(); 
        gettimeofday(&begin, NULL);
        EFIELD_HFIELD();
        gettimeofday(&end, NULL);
        t_efield_hfield += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
        // t = t + dt/2.0;
        // printf("Efield hfield done\n");
        t = t + dt;
        add++;
        assn++;
        KRMS=1;
        //!=============================================================
        if(fmod(n,nperdt)==0) {
            KRMS=2;
        }
        RMS(KRMS);
        //!=============================================================
        //printf("%d\n", icpling);
        comp++;
        if(fmod(n,nperdt*nmaxwell)==0){
            if(icpling!=0){
               gettimeofday(&begin, NULL);
                ELEC_DENS();
                gettimeofday(&end, NULL);
      t_elec_dens += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));

            }
            else{
                DTAC=1.0/FREQ;
                TIMD=(float)(n)*inv_nperdt*nmaxwell*DTAC;
                reciprocal_flo++;
                mul+=3;
                assn+=2;
            }
             //if(nx/2>=llim[rank] && nx/2<=ulim[rank])
                 //front();// ! check
         
            //!=============================================================XXX
            KELEC=KELEC+1;
           
            if(fmod(KELEC,NECRIR)==0){
                nec=nec+1;
                //ECRIR();
                // printf("huahihfldkfj l %d",n);
                ani++;
                gettimeofday(&begin, NULL);
                //anim();              
                gettimeofday(&end, NULL);
                t_anim += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
            }

             //!=============================================================XXXXX
           
           
               //printf("%d %f %f %f %f",n,TIMD,DTAC,ACCEL,dnma);

            if(KELEC==1)
                fptr=fopen("RES.out","w");
            if(KELEC>1){ 
                fptr=fopen("RES.out","a");
                if(fptr!=NULL)
                    fprintf(fptr,"%d %f %f %f %f",n,TIMD,DTAC,ACCEL,dnma);
                else
                    goto l100;
                fclose(fptr);
            }
        }

        l100:
        //!=============================================================XXXXX
        nmod=2*nperdt;
        
    }while(1);


     fclose(fptr1);
    fclose(fptr2);
    fclose(fptr3);
    fclose(file_eyi);
    //fptr=fopen("OUTPUT2.txt","w");
    //time(&itime);
    //date=ctime(&itime);
    //fprintf(fptr,"%s %c","Program was ends at",date);
    //fclose(fptr);

    //exit(0);
    printf("OUTPUT 2\n");
    gettimeofday(&total_end, NULL);
    //total_end = end.tv_usec;
    printf("INC_EFIELD time: %f s\n", t_inc_efield);
    printf("ADV_EFIELD_vel time: %f s\n", t_adv_efield_vel);
    printf("MR_MUR time: %f s\n", t_mr_mur);
    printf("ADV_HFIELD time: %f s\n", t_hfield);
    printf("Total time: %f s\n", ((total_end.tv_sec - total_start.tv_sec) + 
        ((total_end.tv_usec - total_start.tv_usec)/1000000.0)));
        free_all();      
    return 0;
}

void EFIELD_HFIELD()
{
    qmdt=qe/cmasse*dt;  
    aa=FNUM*dt/2.0; 
    alpha=(1.0-aa)/(1.0+aa);
    gamma1=1+aa;
    const3=.50*qe*qe/eps0/cmasse;
    const4=dt*dt/4.0/gamma1;
    const7=.25*dte*(1.0+alpha);
    const8=1.0/2.0/gamma1;
    i0=2;
    x0=(i0-1)*ds;
    x0b=(nx-1)*ds;
    assn+=11;
    add+=4;
    mul+=7;
    div_flo+=9;
     for(i=0;i<total_threads;i++){
        done[i]=0;
    }
    #pragma omp parallel private(start,end,i,j,x,sine,sine1,omp2x, betax, const5x, const6x, extk,eytk,omp2y, betay, const5y, const6y)
    {
        int thread_id=omp_get_thread_num();
       
        start=section_size*thread_id;
        if(thread_id==total_threads-1)
        {
            end=nx;
        }
        else
        {
            end=start+section_size;
        }
      
        //---------------------INC EFIELD-------------------------------
        // INC_EFIELD();
        int tempstart=start;
        int tempend=end;
        if(thread_id==0)
        {
            tempstart+=1;
        }
        if(thread_id==total_threads-1)
        {
            tempend+=1;
        }
        
         for(i=tempstart;i<tempend;i++){
            x=(i-1)*ds;
           
            sine=0.0;
            sine1=0.0;
            mul+=2;
            assn+=3;
            comp+=1;
            add++;
            // sineb=0.0;
            // sine1b=0.0;
            if(x<=(x0+c*t))
            {     
                sine = sin(OMEG*(t-(x-x0)/c));
                add+=2;
                mul+=2;
                trig_flo++;
                assn++;
            }
            add+=2;
            mul++;
            comp++;
            if(x<=(x0+c*(t+dt))) 
            {   
                sine1 = sin(OMEG*(t+dt-(x-x0)/c));
                add+=3;
                mul+=2;
                trig_flo++;
                assn++;
            }
            
           #pragma simd
            for(j=1;j<=ny;j++) {
               
                eyi[i][j] =   E0*(sine);
                eyi1[i][j] =  E0*(sine1);
                mul+=2;
                assn+=2;
                // if(n<10)
                // {
                //     // printf("badum %d %0.10lf %lf %lf %lf %lf %0.10lf %0.10lf %0.10lf \n",n,t,eyi[i][j],eyi1[i][j],E0,sine,t,x,x0);
                // }
            }
        }
        //---------------------INC EFIELD-------------------------------

        //---------------------ADV_EFIELD-------------------------------
        // int i;
       
    for(i=start;i<end;i++){
    
        #pragma simd  
        for(j=1;j<=ny-1;j++){
            omp2x=(den[i][j]+den[i+1][j])*const3;
            betax=omp2x*const4;
            const5x=1.0/(1.0+betax);
            const6x=1.0-betax;
            add+=3;
            mul+=2;
            reciprocal_flo++;
            comp++;            
            if(j>=2){
                extk=ext[i][j];

                exs[i][j]=const5x*( exs[i][j]*(const6x)+qe*(den[i][j]+den[i+1][j])*vx[i][j]*const7
                                     -(exi[i][j]+exi1[i][j])*betax+(hzs[i][j]-hzs[i][j-1])*dteds); 

                ext[i][j]=exs[i][j]+exi1[i][j];
                vx[i][j]=vx[i][j]*alpha - qmdt*(ext[i][j]+extk)*const8;
                assn+=4;
                add+=9;
                mul+=9;
            }
            comp++;
            if(i>1){
                omp2y=(den[i][j]+den[i][j+1])*const3;
                betay=omp2y*const4;
                const5y=1.0/(1.0+betay);
                const6y=1.0-betay;

                eytk=eyt[i][j];

                eys[i][j]=const5y*(eys[i][j]*(const6y)+qe*(den[i][j]+den[i][j+1])*vy[i][j]*const7
                                 -(eyi[i][j]+eyi1[i][j])*betay-(hzs[i][j]-hzs[i-1][j])*dteds);

               
                eyt[i][j]=eys[i][j]+eyi1[i][j];
                vy[i][j]=vy[i][j]*alpha - qmdt*(eyt[i][j]+eytk)*const8;
                assn+=8;
                add+=12;
                mul+=11;
                reciprocal_flo++;
               
            }
        }
        MR_MUR(i);
    }
    done[thread_id]=1;

        //---------------------ADV_EFIELD-------------------------------
        //----------------------HFIELD------------------------------------------
       

        for(i=start;i<end-1;i++)
        {
       
            #pragma simd //private(i)
            for(j=1;j<ny;j++)
            {
                double hzsp=hzs[i][j];
                hzs[i][j]=hzs[i][j]+(-(eys[i+1][j]-eys[i][j])+(exs[i][j+1]-exs[i][j]))*dtmds;
                assn+=2;
                add+=5;
                mul++;            
            }
        }
       	#pragma omp barrier
        #pragma simd
        for(j=1;j<ny;j++)
        {
            double hzsp=hzs[i][j];
            hzs[i][j]=hzs[i][j]+(-(eys[i+1][j]-eys[i][j])+(exs[i][j+1]-exs[i][j]))*dtmds;
            assn+=2;
            add+=5;
            mul++;
        }

        //----------------------HFIELD------------------------------------------

        //----------------------inc EFIELD------------------------------------------
        // i0=2;
        // x0=(i0-1)*ds;
        // x0b=(nx-1)*ds;
        // // #pragma omp parallel for
        // //PRAGMA j - private
        // int tempstart=start;
        // if(tempstart==0)
        //     tempstart+=1;
        // for(i=tempstart;i<=end;i++){
        //     x=(i-1)*ds;
        //     sine=0.0;
        //     sine1=0.0;
        //     // sineb=0.0;
        //     // sine1b=0.0;
        //     if(x<=(x0+c*t))  
        //         sine = sin(OMEG*(t-(x-x0)/c));
        //     if(x<=(x0+c*(t+dt)))  
        //         sine1 = sin(OMEG*(t+dt-(x-x0)/c));
            
        //     //    if(x>=x0b-c*t)  sineb = sin(OMEG*(t+(x-x0b)/c));
        //     //    if(x>=(x0b-c*(t+dt)))  sine1b = sin(OMEG*(t+dt+(x-x0b)/c));
        //     #pragma simd
        //     for(j=1;j<=ny;j++) {
        //         //      eyi[i][j] =   E0*(sine+sineb);
        //         //eyi1[i][j] =  E0*(sine1+sine1b);
        //         //eyi1[i][j] =  E0*(sine1b);
        //         //file_xelec=fopen("file_xelec.dat","w");

        //         eyi[i][j] =   E0*(sine);
        //         eyi1[i][j] =  E0*(sine1);
        //     }
        // }

        // INC_EFIELD();
        //----------------------inc EFIELD------------------------------------------
        //----------------------adv EFIELD------------------------------------------
 //        // double aa,betax,betay,alpha,gamma1,extk,eytk,const7,const8;
 //        // double omp2x,omp2y,qmdt,const3,const4,const5x,const6x,const5y,const6y;


 //        qmdt=qe/cmasse*dt;  
 //        aa=FNUM*dt/2.0; 
 //        alpha=(1.0-aa)/(1.0+aa);
 //        gamma1=1+aa;
 //        const3=.50*qe*qe/eps0/cmasse;
 //        const4=dt*dt/4.0/gamma1;
 //        const7=.25*dte*(1.0+alpha);
 //        const8=1.0/2.0/gamma1;

 //        //printf("const5: %f \n", const5);

 //        //          const5=0.25d0*(1.0d0+alpha)*dte
 //        //          const6=qmdt*.5d0/gamma1
 //       int i;
 //    // #pragma omp parallel for private(i,omp2x, betax, const5x, const6x, extk,eytk,omp2y, betay, const5y, const6y)
 //     //for private(omp2, beta, const5, const6, extk,eytk,i) schedule(dynamic)
 //    for(i=start;i<end;i++){
 //        #pragma simd  //private(omp2, beta, const5, const6, extk,eytk)  //private(omp2, beta, const5, const6, extk,i)
 // //       #pragma omp parallel for private(omp2, beta, const5, const6, extk,i) schedule(dynamic)
 //        //PRAGMA
 //        for(j=1;j<=ny-1;j++){
 //            omp2x=(den[i][j]+den[i+1][j])*const3;
 //            betax=omp2x*const4;
 //            const5x=1.0/(1.0+betax);
 //            const6x=1.0-betax;
            
 //            if(j>=2){
 //                extk=ext[i][j];

 //                exs[i][j]=const5x*( exs[i][j]*(const6x)+qe*(den[i][j]+den[i+1][j])*vx[i][j]*const7
 //                                     -(exi[i][j]+exi1[i][j])*betax+(hzs[i][j]-hzs[i][j-1])*dteds); 

 //                ext[i][j]=exs[i][j]+exi1[i][j];
 //                vx[i][j]=vx[i][j]*alpha - qmdt*(ext[i][j]+extk)*const8;
 //            }
 //            if(i>=1){
 //                omp2y=(den[i][j]+den[i][j+1])*const3;
 //                betay=omp2y*const4;
 //                const5y=1.0/(1.0+betay);
 //                const6y=1.0-betay;

 //                eytk=eyt[i][j];

 //                eys[i][j]=const5y*(eys[i][j]*(const6y)+qe*(den[i][j]+den[i][j+1])*vy[i][j]*const7
 //                                 -(eyi[i][j]+eyi1[i][j])*betay-(hzs[i][j]-hzs[i-1][j])*dteds);

 //                 if(n==1)
 //                {
 //                    // printf("blalala %lf %lf %lf %lf %lf %lf %lf %lf %lf\n\n",eys[i][j],const6y,den[i][j],den[i][j+1],vy[i][j],eyi[i][j],eyi[i][j+1],hzs[i][j],hzs[i-1][j]);
 //                }
 //                eyt[i][j]=eys[i][j]+eyi1[i][j];
 //                vy[i][j]=vy[i][j]*alpha - qmdt*(eyt[i][j]+eytk)*const8;
               
 //            }
 //        }
 //    }
    
        //----------------------adv EFIELD------------------------------------------

        





    }
    // MR_MUR();
}
/*

void INC_EFIELD(){
        //printf("In INC_EFIELD\n");
        // float sine,sine1,x,x0;
        // float sineb,sine1b,x0b;
        // int i0;
        i0=2;
        x0=(i0-1)*ds;
        // x0b=(nx-1)*ds;
        // #pragma omp parallel for
        //PRAGMA j - private
        #pragma omp parallel for private(x,sine,sine1,j)
        for(i=1;i<=nx;i++){
            x=(i-1)*ds;
            sine=0.0;
            sine1=0.0;
            // sineb=0.0;
            // sine1b=0.0;
            if(x<=(x0+c*t))  
                sine = sin(OMEG*(t-(x-x0)/c));
            if(x<=(x0+c*(t+dt)))  
                sine1 = sin(OMEG*(t+dt-(x-x0)/c));
            
            //    if(x>=x0b-c*t)  sineb = sin(OMEG*(t+(x-x0b)/c));
            //    if(x>=(x0b-c*(t+dt)))  sine1b = sin(OMEG*(t+dt+(x-x0b)/c));
            #pragma simd
            for(j=1;j<=ny;j++) {
                //      eyi[i][j] =   E0*(sine+sineb);
                //eyi1[i][j] =  E0*(sine1+sine1b);
                //eyi1[i][j] =  E0*(sine1b);
                //file_xelec=fopen("file_xelec.dat","w");

                eyi[i][j] =   E0*(sine);
                eyi1[i][j] =  E0*(sine1);
            }
        }
    }
    */

    void INC_EFIELD(){
        //printf("In INC_EFIELD\n");
        float sine,sine1,x,x0;
        float sineb,sine1b,x0b;
        int i0;
        i0=2;
        x0=(i0-1)*ds;
        x0b=(nx-1)*ds;
        // #pragma omp parallel for
        //PRAGMA j - private
        for(i=1;i<=nx;i++){
            x=(i-1)*ds;
            // printf("x= %lf\n",x);
            sine=0.0;
            sine1=0.0;
            sineb=0.0;
            sine1b=0.0;
            if(x<=(x0+c*t))
            { 
                
                sine = sin(OMEG*(t-(x-x0)/c));
               
            }
            if(x<=(x0+c*(t+dt))) 
            {   
                sine1 = sin(OMEG*(t+dt-(x-x0)/c));

                
            }
            
            for(j=1;j<=ny;j++) {
              

                eyi[i][j] =   E0*(sine);
                eyi1[i][j] =  E0*(sine1);
                
            }
        }
    }



//     void ADV_EFIELD_vel(){
//         //printf("In ADV_EFIELD_vel\n");
//         double aa,betax,betay,alpha,gamma1,extk,eytk,const7,const8;
//         double omp2x,omp2y,qmdt,const3,const4,const5x,const6x,const5y,const6y;


//         qmdt=qe/cmasse*dt;  
//         aa=FNUM*dt/2.0; 
//         alpha=(1.0-aa)/(1.0+aa);
//         gamma1=1+aa;
//         const3=.50*qe*qe/eps0/cmasse;
//         const4=dt*dt/4.0/gamma1;
//         const7=.25*dte*(1.0+alpha);
//         const8=1.0/2.0/gamma1;

//         //printf("const5: %f \n", const5);

//         //          const5=0.25d0*(1.0d0+alpha)*dte
//         //          const6=qmdt*.5d0/gamma1
//        int i;
//     #pragma omp parallel for private(i,omp2x, betax, const5x, const6x, extk,eytk,omp2y, betay, const5y, const6y)
//      //for private(omp2, beta, const5, const6, extk,eytk,i) schedule(dynamic)
//     for(i=0;i<=nx-1;i++){
//         #pragma simd  //private(omp2, beta, const5, const6, extk,eytk)  //private(omp2, beta, const5, const6, extk,i)
//  //       #pragma omp parallel for private(omp2, beta, const5, const6, extk,i) schedule(dynamic)
//         //PRAGMA
//         for(j=1;j<=ny-1;j++){
//             omp2x=(den[i][j]+den[i+1][j])*const3;
//             betax=omp2x*const4;
//             const5x=1.0/(1.0+betax);
//             const6x=1.0-betax;
            
//             if(j>=2){
//                 extk=ext[i][j];

//                 exs[i][j]=const5x*( exs[i][j]*(const6x)+qe*(den[i][j]+den[i+1][j])*vx[i][j]*const7
//                                      -(exi[i][j]+exi1[i][j])*betax+(hzs[i][j]-hzs[i][j-1])*dteds); 

//                 ext[i][j]=exs[i][j]+exi1[i][j];
//                 vx[i][j]=vx[i][j]*alpha - qmdt*(ext[i][j]+extk)*const8;
//             }
//             if(i>=1){
//                 omp2y=(den[i][j]+den[i][j+1])*const3;
//                 betay=omp2y*const4;
//                 const5y=1.0/(1.0+betay);
//                 const6y=1.0-betay;

//                 eytk=eyt[i][j];

//                 eys[i][j]=const5y*(eys[i][j]*(const6y)+qe*(den[i][j]+den[i][j+1])*vy[i][j]*const7
//                                  -(eyi[i][j]+eyi1[i][j])*betay-(hzs[i][j]-hzs[i-1][j])*dteds);

//                 eyt[i][j]=eys[i][j]+eyi1[i][j];
//                 vy[i][j]=vy[i][j]*alpha - qmdt*(eyt[i][j]+eytk)*const8;
//             }
//         }
//     }

// }

//     void ADV_HFIELD(){
//         //printf("In ADV_HFIELD\n");
//         /*
//         Subroutine to step 
//         currently for vacuum
//         */

//         int i;
//         #pragma omp parallel for private(j)//for private(i) schedule(dynamic)
//         for(i=0;i<nx;i++)
//         {
//         //PRAGMA
//     //        #pragma omp parallel for private(i) schedule(dynamic)
//             #pragma simd //private(i)
//             for(j=1;j<ny;j++)
//                 hzs[i][j]=hzs[i][j]+(-(eys[i+1][j]-eys[i][j])+(exs[i][j+1]-exs[i][j]))*dtmds;
//         }
        
//     }

    void RMS(int k){
        //printf("In RMS\n");
        double z1,z2;

        int i;
        #pragma omp parallel for private(z1,z2,j)
        for(i=1;i<nx;i++)
        {        
        // #pragma omp parallel for private(z1,z2,i) schedule(dynamic)
        //PRAGMA
            #pragma simd
            for(j=1;j<ny;j++){
                z1=(ext[i][j]*ext[i][j]+ext[i-1][j]*ext[i-1][j])*.5f;
                z2=(eyt[i][j]*eyt[i][j]+eyt[i][j-1]*eyt[i][j-1])*.5f;
                erms2[i][j]=erms2[i][j]+(z1+z2)*inv_nperdt;
                assn+=3;
                add+=4;
                mul+=7;
                // exrms2[i][j]=exrms2[i][j]+(z1)/(float)(nperdt);
                // eyrms2[i][j]=eyrms2[i][j]+(z2)/(float)(nperdt);
                //erms[i][j] = (z1+z2)/(float)(nperdt);
            }
        }

        if(k==2){
        int j;
        #pragma omp parallel for private(j)
        for(i=0;i<nx;i++){
            // #pragma omp parallel for private(j) schedule(dynamic)
            //PRAGMA
                #pragma simd
                for(j=0;j<ny;j++){
                    // printf("%e 2\n", (float)erms2[i][j]);
                    // printf("%e 3\n", (float)exrms2[i][j]);
                    // printf("%e 4\n", (float)eyrms2[i][j]);
                    ERMS[i][j] = sqrt(erms2[i][j]);  
                    erms2[i][j]=0.0f;
                    sqrt_flo++;
                    assn++;
                    // exrms[i][j] = sqrt(exrms2[i][j]);
                    // exrms2[i][j]=0.0;
                    // eyrms[i][j] = sqrt(eyrms2[i][j]);  
                    // eyrms2[i][j]=0.0;
                } 
        }
        }
    }

    // void MR_MUR(){
    //     //printf("In MR_MUR\n");
    //     /*
    //        Subroutine to apply Mr. Murs absorbing boundary conditions
    //     */
    //     double csym;
    //     // double aaa[2012],bbb[2012],ccc[2012],ddd[2012];
    //     // double fff[2012],ggg[2012],hhh[2012],jjj[2012];

    //     csym=1.0;
    //     if(IABSOR==2)
    //     csym=0.0;

    //     int j,i;
    //     // #pragma omp parallel for private(j) schedule(dynamic)
    //     #pragma simd
    //     for(j=1;j<=ny;j++){
    //         eys[1][j]=eys1[2][j]+c1*(eys[2][j]-eys1[1][j]);          
    //         eys[nx][j]=eys1[nx-1][j]+c1*(eys[nx-1][j]-eys1[nx][j]);
    //     }

    //     // #pragma omp parallel for private(i) schedule(dynamic)
    //     //PRAGMA
    //     #pragma simd
    //     for(i=0;i<=nx;i++){
    //         exs[i][1]=exs1[i][2] +csym*c1*(exs[i][2]-exs1[i][1]);
    //         exs[i][ny]=exs1[i][ny-1] +csym*c1*(exs[i][ny-1]-exs1[i][ny]);
    //     }

    //     #pragma omp parallel for private(j)
    //     for(i=0;i<=nx;i++){
    //     // #pragma omp parallel for private(i) schedule(dynamic)
    //     //PRAGMA
    //         #pragma simd
    //         for(j=1;j<=ny;j++){
    //             exs1[i][j]=exs[i][j];
    //             eys1[i][j]=eys[i][j];
    //         }
    //     }
    // }

    void MR_MUR(int row)
    {
        double csym;
        csym=1.0;
        assn++;
        if(IABSOR==2){
            csym=0.0;
            assn++;
        }

        int column=0;
        if(row==2)
        {   
            for(column=1;column<=ny;column++){
                eys[1][column]=eys1[2][column]+c1*(eys[2][column]-eys1[1][column]);
                eys1[1][column]=eys[1][column];  
                eys1[2][column]=eys[2][column];
                assn+=3;
                add+=2;
                mul++;      

            // eys[nx][j]=eys1[nx-1][j]+c1*(eys[nx-1][j]-eys1[nx][j]);
            } 
        }
        if(row==nx-1)
        {
            for(column=1;column<=ny;column++){
                // eys[1][j]=eys1[2][j]+c1*(eys[2][j]-eys1[1][j]);          
            eys[nx][column]=eys1[nx-1][column]+c1*(eys[nx-1][column]-eys1[nx][column]);
            eys1[nx][column]=eys[nx][column];
            eys1[nx-1][column]=eys[nx-1][column];
            assn+=3;
            add+=2;
            mul++;
            }   
        }
        // else{
            exs[row][1]=exs1[row][2] +csym*c1*(exs[row][2]-exs1[row][1]);
            exs1[row][1]=exs[row][1];
            exs1[row][2]=exs[row][2];
            exs[row][ny]=exs1[row][ny-1] +csym*c1*(exs[row][ny-1]-exs1[row][ny]);
            exs1[row][ny]=exs[row][ny];
            exs1[row][ny-1]=exs[row][ny-1];
            assn+=6;
            add+=4;
            mul+=4;

        // }
    }

    double FIONIZ(double EE,double PR){
        //printf("In FIONIZ\n");
        double fioniz,ARG,VD;
        double AA1,BB1,AA2,BB2,AA3,BB3,EE3;
        //c-------------------------------
        AA1=15.0*1.0e2;
        BB1=365.0*1.0e2;
        AA2=8.8050*1.0e2;
        BB2=258.450*1.0e2;
        AA3=0.0050*1.0e2;
        BB3=200.0*1.0e2;
        EE3=3.125e3;
        assn+=7;
        mul+=6;
        //c-------------------------------
        
        amu=QSM/FNUM*PRESSURE/PR;
        VD=amu*EE*PR;
        assn+=2;
        mul+=3;
        div_flo+=2;
        comp++;

        if(EE>2.0e4){
        	ARG=BB1/EE;
        	fioniz=AA1*PR*exp(-ARG)*VD;
            assn+=2;
            mul+=3;
            div_flo++;
            exp_flo++;             
        } 
        else if(EE<5.0e3){
        	ARG=BB3*(1.0/EE-1.0/EE3);
        	fioniz=AA3*PR*(exp(-ARG)-1.0)*VD;
            assn+=2;
            add+=3;
            mul+=4;
            div_flo+=2;
            exp_flo++;
        } 
        else{
        	ARG=BB2/EE;
        	fioniz=AA2*PR*exp(-ARG)*VD;
            assn+=2;
            add++;
            mul+=3;
            div_flo++;
            exp_flo++;
        }
        return fioniz;
    }

    void anim(){
        //printf("In anim\n");
        double e_total;
        char fil[50];
        char* buffer = (char *)malloc(sizeof(int));
        snprintf(buffer, sizeof(buffer) - 1, "%d", ani);
        strcpy(fil,"anim/");
        strcat(fil, buffer);
        strcat(fil,".dat");

        file_xelec=fopen(fil,"w");
        
        for(j=2;j<ny;j++){
            for(i=2;i<nx;i++){
            //      e_total = sqrt(pow(eyt[i][j],2) + pow(ext[i][j],2));
            //e_total = sqrt(pow(ext[i][j],2));
                // printf("animmmmmmm %lf\n",den[i][j]);
            fprintf(file_xelec,"%e ",(den[i][j]));
        }
        fprintf(file_xelec,"\n");
        }
        fclose(file_xelec);
    }

    void ELEC_DENS(){
        //printf("In ELEC_DENS\n");
        /*
        routine to calculate electron DENsity
        free diffusion + ambipolar dIFfusion + ioniZation
        */

        int ie,ied,je,jed,iii;
        double fioniz,aad,taumij;
        double coref,ee,cf,da,dac,fnui,fnua;
        double omgc2,rec1;
        double d0,dimax,ecm;
        double frqij,dtacmax,tcycle;

        //!=============================================================
        dnma=0.0;   //! max density at grids
        //FNUIMAX=0.0;
        dimax=0.0;

        ACCEL=naccel; 
        //      ACCEL=FLOAT(naccel)*(1.-(AMIN1(0.2*PARC,1.))**0.2*(1.-1./ACCEL))  
        assn+=3;
        comp++;

        if(ACCEL<=1){
            ACCEL=1.0;
            assn++;
        }

        // printf("%e 1\n", nu2omeg2);
        coref=FNUM/sqrt(nu2omeg2);
        amu=EMOB;
        assn+=2;
        div_flo++;
        sqrt_flo++;
        //c---------------------------------------------------------
        // #pragma omp parallel for 
        int i;
        #pragma omp parallel for private(j,EDIF,ETEM,ee,cf,frqij,ecm,taumij,aad) //for private(EDIF,ETEM,ee,cf,frqij,ecm,taumij,aad,i) schedule(dynamic)
        
        
        for(i=1;i<=nx;i++){
    //        #pragma omp parallel for private(EDIF,ETEM,ee,cf,frqij,ecm,taumij,aad,i) schedule(dynamic)
            #pragma simd//private(EDIF,ETEM,ee,cf,frqij,ecm,taumij,aad) //private(EDIF,ETEM,ee,cf,frqij,ecm,taumij,aad)
            for(j=2;j<ny;j++){ 
                denp[i][j]=den[i][j];
                ee=ERMS[i][j];
                ee=ee/PRESSURE*coref;
                assn+=3;
                mul++;
                div_flo++;
                
                //c------------ calculation of ionization frequency---------
                frqio[i][j]=FIONIZ(ee,PRESSURE);
                //c---------------------------------------------------------
                
                cf=1.0;
                assn++;
                comp++;
                if(den[i][j]<=0.1){ 
                    cf=0.0;
                    assn++;
                }
                frqij=frqio[i][j];
                assn++;
                //c----- calculation of diffusion coefficient--------------

                EDIF=dife;
                assn++;
                comp++;
                if(frqij<0){
                    frqij=0.0;
                    assn++;        
                }
                comp++;
                if(cene==0.0){
                    ecm=ee*0.01;    
                    ETEM=0.1+2.0/3.0*pow((0.0021*ecm*(91.0+ecm)),0.3333);
                    EDIF=EMOB*ETEM;
                    assn+=3;
                    add+=2;
                    mul+=5;
                    div_flo++;
                    pow_flo++;
                }
                comp++;
                if(cene<0.0)
                {
                    taumij=eps0/(qe*(den[i][j]+1.0)*amu);
                    aad=frqij*taumij;
                    EDIF=(aad*dife+difa)/(aad+1.0);
                    assn+=3;
                    add+=3;
                    mul+=4;
                    div_flo+=2;
                }
            
                //c---------------------------------------------------------              
                DIFFUSION[i][j]=EDIF;
             //EDIT   #pragma omp critical
                    dimax=Max(dimax,DIFFUSION[i][j]);
                // SIO[i][j]=frqio[i][j]*den[i][j];
                assn++;
                max_min_flo++;
            }
        
    }
        //c---------------------------------------------------------              

        //c------- time step for fluid equation------------
        dtacmax=0.20*ds*ds/dimax;
        tcycle=1.0/FREQ;
        DTAC=Min(dtacmax, ndifmax*tcycle*ACCEL);
        assn+=3;
        mul+=4;
        div_flo+=2;
        max_min_flo++;
        //c       DTAC=AMIN1(DTAC,0.2*ds*ds/dimax/FLOAT(NDEN)/FLOAT(NDEN));
        //c-----------------------------------------------------      
        //c------- storing the previous density values------
       //  #pragma omp parallel for private(i) schedule(dynamic)   //This loop is merged with the upar wala loop
       //  for(j=2;j<ny;j++)
       // //     #pragma omp parallel for private(i) schedule(dynamic)
       //      #pragma simd //private(i)
       //      for(i=1;i<nx;i++)
       //          denp[i][j]=den[i][j];

        //c---------- density calculation  --------- 
        #pragma omp parallel for private(j,da,dac,rec1,fnua,fnui) //for private(da,dac,rec1,fnua,fnui,i) schedule(dynamic)     
        for(i=1;i<nx;i++){
       //     #pragma omp parallel for private(da,dac,rec1,fnua,fnui,i) schedule(dynamic)
            #pragma simd //private(da,dac,rec1,fnua,fnui)//private(da,dac,rec1,fnua,fnui)
            for(j=2;j<ny;j++){
                da=DIFFUSION[i][j];
                dac=da*DTAC/(ds*ds);
                
                //c          IF(IDIEL[i][j]!=0) EXIT
                rec1=RECOMB*den[i][j]*DTAC;
                
                //c------------ ionization and attachment frequency-------------
                fnua=0.0;
                fnui=frqio[i][j];
                assn+=4;
                mul+=4;
                div_flo++;
                comp++;
                if(fnui<0.0){
                    fnua=-fnui;
                    fnui=0.0;
                    assn+=2;
                }

                //c---------- Density equation updates------------------------

                den[i][j]=denp[i][j]*exp(fnui*DTAC)+dac*(denp[i+1][j]+denp[i-1][j]+denp[i][j+1]+denp[i][j-1]-4.0*denp[i][j]);

                den[i][j]=den[i][j]/(1.0+rec1+fnua*DTAC);
                assn+=2;
                add+=7;
                mul+=5;
                div_flo++;
                exp_flo++;
                comp++;

                if(den[i][j]<=1e-15){
                    den[i][j]=den[i][j]*0.0;
                    mul++;
                    assn++;
                }
		
                dnma=Max(dnma,den[i][j]);
                dent[i][j]=dent[i][j]+den[i][j]*DTAC;
                SIO[i][j]=frqio[i][j]*denp[i][j];
                SIOT[i][j]=SIOT[i][j]+SIO[i][j]*DTAC;
                max_min_flo++;
                assn+=4;
                add+=2;
                mul+=3;
                //c         FNUIMAX=AMAX1(FNUIMAX,fnui) 
            }
        
    }

        //c---------------Actual time calculation -------------------------
        TIMD=TIMD+DTAC;
        //c------------------------------------------------------------
        
        omgc2=qe*qe/cmasse/eps0*dnma;
        PARC=omgc2/(pow(OMEG,2)+pow(FNUM,2));
        assn+=3;
        add+=2;
        mul+=3;
        div_flo+=2;
        pow_flo+=2;
        
        // int j;
     //    #pragma omp parallel for private(i) schedule(dynamic)
     //    for(j=2;j<ny;j++){
     // //       #pragma omp parallel for private(i) schedule(dynamic)
     //        #pragma simd //private(i) //collapse(2)
     //        for(i=1;i<nx;i++){
     //            dent[i][j]=dent[i][j]+den[i][j]*DTAC;
     //            SIOT[i][j]=SIOT[i][j]+SIO[i][j]*DTAC;
     //        }
     //    }
    }

    void SETUP(){
        //printf("In SETUP\n");
        static int k;
        static double ardiy,ardix,dinig;
        static double xd0,yd0,xxi,yyj;
        static double aaa[2012],bbb[2012],ccc[2012],ddd[2012];

        istart = 1;
        iend = nx-1;

        dt=1.0/(float)(nperdt)/FREQ;
        ds=c/(float)(nlamb)/FREQ;

        OMEG=2.0*pi*FREQ;
        QSM=qe/cmasse;
        e2_epsm=qe*qe/eps0/cmasse;
        
        /*
        !  AIR DAT
        !  FNUM=electron-neutral coll frequency
        !  RECOMB=electron-ion recombination coefficient
        !  EMOB=electron mobility
        !  ETEM=electron temperature
        !  DIFE=electron diffusion coeff
        */

        FNUM=5.3*pow(10,9)*PRESSURE;
        RECOMB=crec*1.0*pow(10,-13);
        EMOB=QSM/FNUM;
        ETEM=2.0*abs(cene);
        dife=EMOB*ETEM;
        difa=dife/100.0;
        nu2omeg2=pow(OMEG,2)+pow(FNUM,2);
        nu_omeg=FNUM/OMEG;
        inv_nperdt = 1.0/nperdt;
        //!================Initial density location =======================

        for(K=1;K<=nini;K++)
        {
        if(xmid[K]>0) 
            imid[K] = xmid[K]*nx;
        if(ymid[K]>0) 
            jmid[K] = ymid[K]*ny;
        if(xmid[K]==0) 
            imid[K] = nx*0.5;
        if(ymid[K]==0) 
            jmid[K] =ny/2;
        }

        //!=============================================================
        TEMP0=300.0;
        DENG0=PRESSURE/760.0*101300.0/akb/TEMP0;

        radius = nx/5*ds;
        printf("radius: %f\n", radius);
        printf("ds: %f\n", ds);

        int i;
        #pragma omp parallel for private(j)//for private(i) schedule(dynamic)
        
        
        for(i=1;i<=nx;i++){
      //      #pragma omp parallel for private(i) schedule(dynamic)
            #pragma simd //private(i)
            for(j=1;j<=ny;j++)     
            {
                den[i][j]=0.0;
                ERMS[i][j]=E0/sqrt(2.0);
            }
        }
    

        // for(j=1;j<=ny;j++)
        // {
        //   for(i=1;i<=nx;i++)      
        //   {
        //     ERMS[i][j]=E0/sqrt(2.0);
        //   }
        // }      

        for(K=1;K<=nini;K++){
        xd0=ds*imid[K];
        yd0=ds*jmid[K];

        if(xmid[K]<0) 
            xd0=-xmid[K];
        if(ymid[K]<0) 
            yd0=-ymid[K];

        //!================Initial density, Gaussian, defined =======================
        //!make DEN and DENP =0
        int j;
        //PRAGMA j - private
        // #pragma simd private(yyj,ardiy,dinig,j,xxi,ardix)
     // #pragma omp parallel for private(yyj,ardiy,dinig,j,xxi,ardix) schedule(dynamic)  //change
        for(i=1;i<=nx;i++){
            xxi=ds*i;
            ardix=0.0;
            if(sgdx0[K]>0)
                ardix=(-pow((xxi-xd0),2))/2.0/sgdx0[K]/sgdx0[K];
            // #pragma omp parallel for private(yyj,ardiy,dinig,j) schedule(dynamic)
                // #pragma simd //private(yyj,ardiy,dinig,/*j,*/xxi,ardix)
                for(j=1;j<=ny;j++){    
                    yyj=ds*j;
                    ardiy=0.0;
                    if(sgdy0[K]>0) 
                       ardiy=-pow((yyj-yd0),2)/2.0/sgdy0[K]/sgdy0[K];
                    dinig=DINI[K]*exp(ardix+ardiy);
                    if(dinig<=DINI[K]*1.0*exp(-2)) 
                        dinig=0;
                    den[i][j]=den[i][j]+dinig;
                    //denp[i][j]=den[i][j];
                }
        }
        }

      //c     find relative delay for x, y, Z cell displacement

        dte=dt/eps0;
        dtm=dt/xmu0;
        dteds=dte/ds;
        dtmds=dtm/ds;
      
        //mur constants
        c1=(c*dt-ds)/(c*dt+ds);
        c2=2.0*ds/(c*dt+ds);
        c3=(c*dt*c*dt)/(2.0*ds*(c*dt+ds));
    }

    void ZERO(){
        //printf("In ZERO\n");
        int i;
       #pragma omp parallel for private(j) //for private(j) schedule(dynamic)
        for(i=1;i<=nx;i++){   
            #pragma simd //private(i)
            for(j=1;j<=ny;j++) {        //nyy changed to ny
                ERMS[i][j]=0.0;
                erms2[i][j]=0.0;
                exrms[i][j]=0.0;
                exrms2[i][j]=0.0;
                eyrms[i][j]=0.0;
                eyrms2[i][j]=0.0;
                den[i][j]=0.0;
                dent[i][j]=0.0;
                denp[i][j]=0.0;
                exi[i][j]=0.0;
                eyi[i][j]=0.0;
                exi1[i][j]=0.0;
                eyi1[i][j]=0.0;
                exs[i][j]=0.0;
                eys[i][j]=0.0;
                hzi[i][j]=0.0;
                hzs[i][j]=0.0;
                vx[i][j]=0.0;
                vy[i][j]=0.0;
                ext[i][j]=0.0;
                eyt[i][j]=0.0;
                exs1[i][j]=0.0;
                eys1[i][j]=0.0;
                exs2[i][j]=0.0;
                eys2[i][j]=0.0;
                frqio[i][j]=0.0;
                DIFFUSION[i][j]=0.0;
                SIO[i][j]=0.0;
                SIOT[i][j]=0.0;
            }
        }
    
        PARC=0.0;
    }

void free_all(){
        free(ERMS);
        free(erms2);
        free(exrms);
        free(exrms2);
        free(eyrms);
        free(eyrms2);
        free(erms);
        free(den);
        free(dent);
        free(exi);
        free(eyi);
        free(exi1);
        free(eyi1);
        free(exs);
        free(eys);
        free(hzi);
        free(hzs);
        free(vx);
        free(vy);
        free(ext);
        free(eyt);
        free(exs1);
        free(exs2);
        free(eys1);
        free(eys2);
        free(xmid);
        free(ymid);
        free(sgdx0);
        free(sgdy0);
        free(DINI);
        free(DIFFUSION);
        free(SIO);
        free(SIOT);
        free(pow1);
        free(frqio);
        free(puis);
        free(denp);
        free(frq);
        free(done);
        free(imid);
        free(jmid);
}