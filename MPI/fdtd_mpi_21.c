#include "xyfdtd.h"

int KELEC,KRMS,nmod,NECRIR,icpling,K;
double output[421][2012];
double slambx,slamby,tstop, radius = 0;
time_t itime;
int start, end;
double sine, sine1, x, x0;
double sineb,sine1b,x0b;
int i0, section_size, total_threads;
double aa,betax,betay,alpha,gamma1,extk,eytk,const7,const8;
double omp2x,omp2y,qmdt,const3,const4,const5x,const6x,const5y,const6y;
      
FILE *fptr, *file_a, *file_xelec, *file_yelec;

void SETUP();
void RMS(int k);
void MR_MUR(int row);
double FIONIZ(double EE,double PR);
void ELEC_DENS();
void ZERO();
void anim();
void initiate_array();
void free_all();
void EFIELD_HFIELD();
double t1, t2, t_inc_efield, t_adv_efield_vel, t_mr_mur, t_hfield,t_efield_hfield,t_elec_dens,t_anim,t_rms,t_zero;
int rank, ierr, num_procs, tag = 1;
long long llim[1024], ulim[1024], size, length;
double *temp_send_arr, *temp_recv_arr;
MPI_Status status;
int main(int argc, char **argv) {
    ierr = MPI_Init(&argc, &argv);
    if(ierr!=0){
        printf("Error in Initialization\n");
        exit(0);
    }
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    printf("%d %d\n", rank, num_procs);
    
    initiate_array();
    struct timeval begin, end, total_start, total_end;
    gettimeofday(&begin, NULL);
    t1 = begin.tv_usec;
    gettimeofday(&total_start, NULL);
    fptr=fopen("xstart_1024.dat","r");
    x0b=(nx-1)*ds;
    AA1=15.0*1.0e2;
    BB1=365.0*1.0e2;
    AA2=8.8050*1.0e2;
    BB2=258.450*1.0e2;
    AA3=0.0050*1.0e2;
    BB3=200.0*1.0e2;
    EE3=3.125e3;

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

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    xmid[1]=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    ymid[1]=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    sgdx0[1]=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    sgdy0[1]=atof(line);
    len=0;

    getline(&line, &len, fptr);
    getline(&line, &len, fptr);
    len=10;
    getline(&line, &len, fptr);
    DINI[1]=atof(line);
    len=0; 
      
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
    inv_nperdt=1.0f/(double)nperdt;
    inv_c=1.0f/c;
    inv_cmasse=1.0f/cmasse;
    nx1=nx-1;
    ny1=ny-1;
    t=0;
    nec=0;
    ani=0;
    nstep=0;
    //=============================================================

    //=========================MPI related=========================
    size = nx/num_procs;
    length = size + nx%num_procs;
    llim[0] = 0;
    ulim[0] = length - 1;

    for(i=1;i<num_procs;i++){
        llim[i] = ulim[i-1] + 1;
        ulim[i] = llim[i] + size - 1;
    }

    if(ulim[num_procs-1] != nx - 1){
        printf("Value of ulim(num_procs) is not NX\n");
        exit(0);
    }

    fptr=fopen("process_limits","w");
    fclose(fptr);
    MPI_Barrier(MPI_COMM_WORLD);
    fptr=fopen("process_limits","a+");
    fprintf(fptr, "llim: %d, ulim: %d\n",llim[rank], ulim[rank]);
    fclose(fptr);
    //=============================================================

    gettimeofday(&begin, NULL);
    ZERO();
    SETUP();
    gettimeofday(&end, NULL);
    t_zero += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));

    // total_threads=atoi(getenv("OMP_NUM_THREADS"));
    // done[total_threads]=1;
    // printf("total threads %d\n", total_threads);
    // section_size=nx/total_threads;
    gettimeofday(&end, NULL);

    printf("Initialization time: %f s\n", ((end.tv_sec - begin.tv_sec) + 
      ((end.tv_usec - begin.tv_usec)/1000000.0)));
    ELEC_DENS();

    if(rank==0){
        fptr=fopen("output.txt","w");      

        //C--------------- Initial parameters ouput-----
        fprintf(fptr,"nx=%d   ny=%d\n",nx,ny);
        fprintf(fptr,"DS=%e   DT=%e\n",ds,dt*1e15);
        fprintf(fptr,"Freq=%e   Omega=%e\n",FREQ,OMEG);
        fprintf(fptr,"Time period=%e\n",1.0/FREQ);
        fprintf(fptr,"Lambda=%e\n",3.0e8/FREQ);
        fprintf(fptr,"Collision Freq=%f\n",FNUM);
        fprintf(fptr,"Recombination Coef=%f\n",RECOMB);
        fprintf(fptr,"Mobility=%f\n",EMOB);
        fprintf(fptr,"Electron Temp= (eV)%f\n",ETEM);
        fprintf(fptr,"DIFFUSION Coef=%f\n",EDIF);
        fprintf(fptr,"Initial Gas/Neutral density=%f\n",DENG0);
        fprintf(fptr,"Cutoff-density=%f\n",(eps0*cmasse/pow(qe,2))*(pow(OMEG,2)+pow(FNUM,2)));// check formula
        fprintf(fptr, "Tstop %e\n",tstop); 
        fclose(fptr);
    }
    KELEC=0;
    n=0;
    int f=0;

    qmdt=qe*inv_cmasse*dt;  
    aa=FNUM*dt/2.0; 
    alpha=(1.0-aa)/(1.0+aa);
    gamma1=1+aa;
    const3=.50*qe*qe*inv_cmasse/eps0;
    const4=dt*dt/4.0/gamma1;
    const7=.25*dte*(1.0+alpha);
    const8=(1.0/2.0/gamma1)*qmdt;
    i0=2;
    x0=(i0-1)*ds;

    do{
        n=n+1;
        nstep += 1;
       
        if(TIMD>tstop){
            free_all();
            ierr = MPI_Finalize();
            if(rank==0){
                printf("%e %d %d\n", dt,nx,ny); 
                printf("OUTPUT 1\n");
                gettimeofday(&total_end, NULL);
                printf("EFIELD_HFIELD time: %f s\n", t_efield_hfield);
                printf("ELEC_DENS time: %f s\n", t_elec_dens);
                printf("RMS time: %f s\n", t_rms);
                printf("Zero time: %f s\n", t_zero);
                printf("ANIM time: %f s\n",t_anim);
                printf("Total time: %f s\n", ((total_end.tv_sec - total_start.tv_sec) + 
                        ((total_end.tv_usec - total_start.tv_usec)/1000000.0)));
            }
            exit(0);
        }
        gettimeofday(&begin, NULL);
        EFIELD_HFIELD();

        gettimeofday(&end, NULL);
        t_efield_hfield += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
        t = t + dt;

        KRMS=1;
        if(fmod(n,nperdt)==0) {
            KRMS=2;
        }
        gettimeofday(&begin, NULL);
        RMS(KRMS);
        gettimeofday(&end, NULL);
        t_rms += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
        
        if(fmod(n,nperdt*nmaxwell)==0){
            if(icpling!=0){
                gettimeofday(&begin, NULL);
                ELEC_DENS();
                gettimeofday(&end, NULL);
                t_elec_dens += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
            }
            else{
                DTAC=1.0/FREQ;
                TIMD=(double)(n)*inv_nperdt*nmaxwell*DTAC;
            }
            KELEC=KELEC+1;
           
            if(fmod(KELEC,NECRIR)==0){
                nec=nec+1;
                ani++;
                gettimeofday(&begin, NULL);
                anim();
                gettimeofday(&end, NULL);
                t_anim += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));
            }

            if(rank==0){
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
        }
        l100:
        nmod=2*nperdt;
        
    } while(1);

    free_all();
    ierr = MPI_Finalize();

    if(rank==0){
        printf("%e %d %d\n", dt,nx,ny); 
        printf("OUTPUT 1\n");
        gettimeofday(&total_end, NULL);
        printf("EFIELD_HFIELD time: %f s\n", t_efield_hfield);
        printf("ELEC_DENS time: %f s\n", t_elec_dens);
        printf("RMS time: %f s\n", t_rms);
        printf("Zero time: %f s\n", t_zero);
        printf("ANIM time: %f s\n",t_anim);
        printf("Total time: %f s\n", ((total_end.tv_sec - total_start.tv_sec) + 
                ((total_end.tv_usec - total_start.tv_usec)/1000000.0)));
    }
    return 0;
}

void EFIELD_HFIELD() {    
    int istart, iend, icount;
    istart = llim[rank];
    iend = ulim[rank] + 1;
    icount = ny;

    int itile=2;
    int ii=0;
    int i,j;
    #pragma omp parallel for private(i,x,sine,sine1,j,omp2x,betax,const5x,const6x,extk,omp2y,betay,const6y,const5y,eytk) schedule(guided)
    for(ii=istart;ii<iend;ii+=itile){
        for(i=ii;i<ii+itile;i++){
            x=i*ds;          
            sine=0.0;
            sine1=0.0;
            if(x<=(x0+c*t))
            { 
                sine = sinf(OMEG*(t-(x-x0)*inv_c));
           
            }
            if(x<=(x0+c*(t+dt))) 
            {   
                sine1 = sinf(OMEG*(t+dt-(x-x0)*inv_c));

            }

            #pragma simd
            for(j=0;j<ny;j++){
                omp2x=(den[i][j]+den[i+1][j])*const3;
                betax=omp2x*const4;
                const5x=1.0/(1.0+betax);
                const6x=1.0-betax;
                if(j>0 && j<ny-1 && i<nx-1){
                    extk=ext[i][j];

                    exs[i][j]=const5x*( exs[i][j]*(const6x)+qe*(den[i][j]+den[i+1][j])*vx[i][j]*const7
                                         -(exi[i][j]+exi1[i][j])*betax+(hzs[i][j]-hzs[i][j-1])*dteds); 

                    ext[i][j]=exs[i][j]+exi1[i][j];
                    vx[i][j]=vx[i][j]*alpha - (ext[i][j]+extk)*const8;
                   
                }
                if(i>0 && i<nx-1 && j<ny-1){
                    eyi[i][j] =   E0*(sine);
                    eyi1[i][j] =  E0*(sine1);
                    omp2y=(den[i][j]+den[i][j+1])*const3;
                    betay=omp2y*const4;
                    const5y=1.0/(1.0+betay);
                    const6y=1.0-betay;

                    eytk=eyt[i][j];

                    eys[i][j]=const5y*(eys[i][j]*(const6y)+qe*(den[i][j]+den[i][j+1])*vy[i][j]*const7
                                     -(eyi[i][j]+eyi1[i][j])*betay-(hzs[i][j]-hzs[i-1][j])*dteds);

                   
                    eyt[i][j]=eys[i][j]+eyi1[i][j];
                    vy[i][j]=vy[i][j]*alpha - (eyt[i][j]+eytk)*const8;
                }
            }
            MR_MUR(i);
        }
        for(i=ii;i<ii+itile-1;i++){
            if(i != ulim[rank]){
                #pragma simd
                for(j=0;j<ny-1;j++){
                    double hzsp=hzs[i][j];
                    hzs[i][j]=hzs[i][j]+(-(eys[i+1][j]-eys[i][j])+(exs[i][j+1]-exs[i][j]))*dtmds;
                }
            }
        }
    }

    #pragma omp parallel for private(j)
    for(i=istart+itile-1;i<iend-1;i+=itile){
        if(i != ulim[rank]){
            #pragma simd
            for(j=0;j<ny-1;j++)
            {
                double hzsp=hzs[i][j];
                hzs[i][j]=hzs[i][j]+(-(eys[i+1][j]-eys[i][j])+(exs[i][j+1]-exs[i][j]))*dtmds;
            }
        }
    }

    if(rank==0){
        // ================================================
        // for(j=0;j<ny;j++){
        //     temp_send_arr[j] = exs[ulim[rank]][j];
        // }
        // ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        // if(ierr!=0){
        //     printf("Error in Sending, e_h_field\n");
        //     exit(0);
        // }
        // ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        // if(ierr!=0){
        //     printf("Error in Receiving, e_h_field\n");
        //     exit(0);
        // }
        // for(j=0;j<ny;j++){
        //     exs[ulim[rank]+1][j] = temp_recv_arr[j];
        // }
        // // ================================================
        // for(j=0;j<ny;j++){
        //     temp_send_arr[j] = eys[ulim[rank]][j];
        // }
        // ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        // if(ierr!=0){
        //     printf("Error in Sending, e_h_field\n");
        //     exit(0);
        // }
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            eys[ulim[rank]+1][j] = temp_recv_arr[j];
        }
        // ================================================
        for(j=0;j<ny;j++){
            temp_send_arr[j] = ext[ulim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        // ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        // if(ierr!=0){
        //     printf("Error in Receiving, e_h_field\n");
        //     exit(0);
        // }
        // for(j=0;j<ny;j++){
        //     ext[ulim[rank]+1][j] = temp_recv_arr[j];
        // }
        // ================================================
    } else if (rank < num_procs-1){
        // ================================================
        // ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        // if(ierr!=0){
        //     printf("Error in Receiving, e_h_field\n");
        //     exit(0);
        // }
        // for(j=0;j<ny;j++){
        //     exs[llim[rank]-1][j] = temp_recv_arr[j];
        // }
        // for(j=0;j<ny;j++){
        //     temp_send_arr[j] = exs[llim[rank]][j];
        // }
        // ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        // if(ierr!=0){
        //     printf("Error in Sending, e_h_field\n");
        //     exit(0);
        // }
        // ================================================
        // ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        // if(ierr!=0){
        //     printf("Error in Receiving, e_h_field\n");
        //     exit(0);
        // }
        // for(j=0;j<ny;j++){
        //     eys[llim[rank]-1][j] = temp_recv_arr[j];
        // }
        for(j=0;j<ny;j++){
            temp_send_arr[j] = eys[llim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        // ================================================
        // ================================================
        // for(j=0;j<ny;j++){
        //     temp_send_arr[j] = exs[ulim[rank]][j];
        // }
        // ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        // if(ierr!=0){
        //     printf("Error in Sending, e_h_field\n");
        //     exit(0);
        // }
        // ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        // if(ierr!=0){
        //     printf("Error in Receiving, e_h_field\n");
        //     exit(0);
        // }
        // for(j=0;j<ny;j++){
        //     exs[ulim[rank]+1][j] = temp_recv_arr[j];
        // }
        // ================================================
        // for(j=0;j<ny;j++){
        //     temp_send_arr[j] = eys[ulim[rank]][j];
        // }
        // ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        // if(ierr!=0){
        //     printf("Error in Sending, e_h_field\n");
        //     exit(0);
        // }
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            eys[ulim[rank]+1][j] = temp_recv_arr[j];
        }
        // ================================================
        // ================================================
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            ext[llim[rank]-1][j] = temp_recv_arr[j];
        }
        // for(j=0;j<ny;j++){
        //     temp_send_arr[j] = ext[llim[rank]][j];
        // }
        // ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        // if(ierr!=0){
        //     printf("Error in Sending, e_h_field\n");
        //     exit(0);
        // }
        // ================================================
        for(j=0;j<ny;j++){
            temp_send_arr[j] = ext[ulim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        // ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        // if(ierr!=0){
        //     printf("Error in Receiving, e_h_field\n");
        //     exit(0);
        // }
        // for(j=0;j<ny;j++){
        //     ext[ulim[rank]+1][j] = temp_recv_arr[j];
        // }

        // ================================================
    } else{
        // ================================================
        // ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        // if(ierr!=0){
        //     printf("Error in Receiving, e_h_field\n");
        //     exit(0);
        // }
        // for(j=0;j<ny;j++){
        //     exs[llim[rank]-1][j] = temp_recv_arr[j];
        // }
        // for(j=0;j<ny;j++){
        //     temp_send_arr[j] = exs[llim[rank]][j];
        // }
        // ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        // if(ierr!=0){
        //     printf("Error in Sending, e_h_field\n");
        //     exit(0);
        // }
        // ================================================
        // ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        // if(ierr!=0){
        //     printf("Error in Receiving, e_h_field\n");
        //     exit(0);
        // }
        // for(j=0;j<ny;j++){
        //     eys[llim[rank]-1][j] = temp_recv_arr[j];
        // }
        for(j=0;j<ny;j++){
            temp_send_arr[j] = eys[llim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        // ================================================
        // ================================================
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            ext[llim[rank]-1][j] = temp_recv_arr[j];
        }
        // for(j=0;j<ny;j++){
        //     temp_send_arr[j] = ext[llim[rank]][j];
        // }
        // ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        // if(ierr!=0){
        //     printf("Error in Sending, e_h_field\n");
        //     exit(0);
        // }
        // ================================================
    }

    // #pragma omp parallel for private(j)
    // for(i=istart;i<iend;i++){
    //     #pragma simd
    //     for(j=0;j<ny-1;j++){
    //         double hzsp=hzs[i][j];
    //         hzs[i][j]=hzs[i][j]+(-(eys[i+1][j]-eys[i][j])+(exs[i][j+1]-exs[i][j]))*dtmds;
    //     }
    // }
    i = ulim[rank];
    #pragma simd 
    for(j=0;j<ny-1;j++){
        double hzsp=hzs[i][j];
        hzs[i][j]=hzs[i][j]+(-(eys[i+1][j]-eys[i][j])+(exs[i][j+1]-exs[i][j]))*dtmds;                
    }

    if(rank==0){
        // ================================================
        for(j=0;j<ny;j++){
            temp_send_arr[j] = hzs[ulim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        // ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        // if(ierr!=0){
        //     printf("Error in Receiving, e_h_field\n");
        //     exit(0);
        // }
        // for(j=0;j<ny;j++){
        //     hzs[ulim[rank]+1][j] = temp_recv_arr[j];
        // }
        // ================================================
    } else if (rank < num_procs-1){
        // ================================================
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            hzs[llim[rank]-1][j] = temp_recv_arr[j];
        }
        // for(j=0;j<ny;j++){
        //     temp_send_arr[j] = hzs[llim[rank]][j];
        // }
        // ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        // if(ierr!=0){
        //     printf("Error in Sending, e_h_field\n");
        //     exit(0);
        // }
        // ================================================
        // ================================================
        for(j=0;j<ny;j++){
            temp_send_arr[j] = hzs[ulim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        // ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        // if(ierr!=0){
        //     printf("Error in Receiving, e_h_field\n");
        //     exit(0);
        // }
        // for(j=0;j<ny;j++){
        //     hzs[ulim[rank]+1][j] = temp_recv_arr[j];
        // }        
        // ================================================
    } else{
        // ================================================
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            hzs[llim[rank]-1][j] = temp_recv_arr[j];
        }
        // for(j=0;j<ny;j++){
        //     temp_send_arr[j] = hzs[llim[rank]][j];
        // }
        // ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        // if(ierr!=0){
        //     printf("Error in Sending, e_h_field\n");
        //     exit(0);
        // }
        // ================================================
    }
}

void RMS(int k){
    double z1,z2;
    int istart, iend, icount;
    istart = llim[rank];
    iend = ulim[rank] + 1;
    icount = ny;
    if(rank==0)
        istart += 1;

    int i,j;
    #pragma omp parallel for private(z1,z2,j)
    for(i=istart;i<iend;i++) {        
        #pragma simd
        for(j=1;j<ny;j++){
            z1=(ext[i][j]*ext[i][j]+ext[i-1][j]*ext[i-1][j])*.5f;
            z2=(eyt[i][j]*eyt[i][j]+eyt[i][j-1]*eyt[i][j-1])*.5f;
            erms2[i][j]=erms2[i][j]+(z1+z2)*inv_nperdt;
        }
    }

    if(k==2){
        if(rank==0)
            istart -= 1;
        int j;
        #pragma omp parallel for private(j)
        for(i=istart;i<iend;i++){
            #pragma simd
            for(j=0;j<ny;j++){
                ERMS[i][j] = sqrt(erms2[i][j]);  
                erms2[i][j]=0.0f;
            } 
        }
    }
}

void MR_MUR(int row){
    double csym;
    csym=1.0;
    if(IABSOR==2)
        csym=0.0;

    int column=0;
    if(row==1)
    {   
        for(column=0;column<ny;column++){
            eys[0][column]=eys1[1][column]+c1*(eys[1][column]-eys1[0][column]);
            eys1[0][column]=eys[0][column];  
            eys1[1][column]=eys[1][column];      
        } 
    }
    if(row==nx-2)
    {
        for(column=0;column<ny;column++){
            eys[nx-1][column]=eys1[nx-2][column]+c1*(eys[nx-2][column]-eys1[nx-1][column]);
            eys1[nx-1][column]=eys[nx-1][column];
            eys1[nx-2][column]=eys[nx-2][column];
        }   
    }
    
    exs[row][0]=exs1[row][1] +csym*c1*(exs[row][1]-exs1[row][0]);
    exs1[row][0]=exs[row][0];
    exs1[row][1]=exs[row][1];
    exs[row][ny-1]=exs1[row][ny-2] +csym*c1*(exs[row][ny-2]-exs1[row][ny-1]);
    exs1[row][ny-1]=exs[row][ny-1];
    exs1[row][ny-2]=exs[row][ny-2];
}

double FIONIZ(double EE,double PR){
    double fioniz,ARG,VD;

    amu=QSM/FNUM*PRESSURE/PR;
    VD=amu*EE*PR;
    if(EE>2.0e4){
        ARG=BB1/EE;
        fioniz=AA1*PR*exp(-ARG)*VD; 
    } 
    else if(EE<5.0e3){
        ARG=BB3*(1.0/EE-1.0/EE3);
        fioniz=AA3*PR*(exp(-ARG)-1.0)*VD;
    } 
    else{
        ARG=BB2/EE;
        fioniz=AA2*PR*exp(-ARG)*VD;
    }
    return fioniz;
}

void anim(){
    int i,j;
    int istart, iend, icount;
    istart = llim[rank];
    iend = ulim[rank] + 1;
    icount = ny;
    if(rank>0){
        for(i=istart;i<iend;i++){
            for(j=0;j<ny;j++){
                temp_send_arr[j] = den[i][j];
            }
            ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
            if(ierr!=0){
                printf("Error in Sending, e_h_field\n");
                exit(0);
            }
        }
        for(i=istart;i<iend;i++){
            for(j=0;j<ny;j++){
                temp_send_arr[j] = ERMS[i][j];
            }
            ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
            if(ierr!=0){
                printf("Error in Sending, e_h_field\n");
                exit(0);
            }
        }
    }
    else if(rank==0){
        for(i=1;i<num_procs;i++){
            istart = llim[i];
            iend = ulim[i] + 1;
            for(j=istart;j<iend;j++){
                ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
                if(ierr!=0){
                    printf("Error in Receiving, e_h_field\n");
                    exit(0);
                }
                long k;
                for(k=0;k<ny;k++){
                    den[j][k] = temp_recv_arr[k];
                }
            }
        }
        for(i=1;i<num_procs;i++){
            istart = llim[i];
            iend = ulim[i] + 1;
            for(j=istart;j<iend;j++){
                ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
                if(ierr!=0){
                    printf("Error in Receiving, e_h_field\n");
                    exit(0);
                }
                long k;
                for(k=0;k<ny;k++){
                    ERMS[j][k] = temp_recv_arr[k];
                }
            }
        }
        double e_total;
        char fil[50];
        char* buffer = (char *)malloc(sizeof(int));
        snprintf(buffer, sizeof(buffer) - 1, "%d", ani);
        strcpy(fil,"anim_300/elec_dens_data/");
        strcat(fil, buffer);
        strcat(fil,".dat");

        file_xelec=fopen(fil,"w");
        for(j=2;j<ny-1;j++){
            for(i=2;i<nx-1;i++){
                fprintf(file_xelec,"%e ",(den[i][j]));
            }
            fprintf(file_xelec,"\n");
        }
        fclose(file_xelec);

        strcpy(fil,"anim_300/efield_data/");
        strcat(fil, buffer);
        strcat(fil,".dat");

        file_xelec=fopen(fil,"w");
        for(j=2;j<ny-1;j++){
            for(i=2;i<nx-1;i++){
                fprintf(file_xelec,"%e ",(ERMS[i][j]));
            }
            fprintf(file_xelec,"\n");
        }
        fclose(file_xelec);
    }
}

void ELEC_DENS(){
    int i,j;
    int istart, iend, icount;
    istart = llim[rank];
    iend = ulim[rank] + 1;
    icount = ny;
    if(rank==0)
        istart += 1;
    else if(rank==num_procs-1)
        iend -= 1;

    int ie,ied,je,jed,iii;
    double fioniz,aad,taumij;
    double coref,ee,cf,da,dac,fnui,fnua;
    double omgc2,rec1;
    double d0,dimax,ecm;
    double frqij,dtacmax,tcycle;

    dnma=0.0;   //! max density at grids
    dimax=0.0;

    ACCEL=naccel; 
    if(ACCEL<=1) 
        ACCEL=1.0;
    
    coref=FNUM/sqrt(nu2omeg2);
    amu=EMOB;

    #pragma omp parallel for private(j,EDIF,ETEM,ee,cf,frqij,ecm,taumij,aad)
    for(i=istart;i<iend;i++){
        #pragma simd
        for(j=1;j<ny-1;j++){
            denp[i][j]=den[i][j];
            ee=ERMS[i][j];
            ee=ee/PRESSURE*coref;
            
            //c------------ calculation of ionization frequency---------
            frqio[i][j]=FIONIZ(ee,PRESSURE);
            //c---------------------------------------------------------
            
            cf=1.0;
            if(den[i][j]<=0.1) 
                cf=0.0;
            frqij=frqio[i][j];

            //c----- calculation of diffusion coefficient--------------
            EDIF=dife;
            if(frqij<0)
                frqij=0.0;        

            if(cene==0.0){
                ecm=ee*0.01;    
                ETEM=0.1+2.0/3.0*pow((0.0021*ecm*(91.0+ecm)),0.3333);
                EDIF=EMOB*ETEM;
            }

            if(cene<0.0)
            {
                taumij=eps0/(qe*(den[i][j]+1.0)*amu);
                aad=frqij*taumij;
                EDIF=(aad*dife+difa)/(aad+1.0);
            }
        
            DIFFUSION[i][j]=EDIF;
            dimax=Max(dimax,DIFFUSION[i][j]);
        }
    }

    int temp_start = istart - 1;
    int temp_end = iend + 1;
    // if(rank==0){
    //     temp_start++;
    // } else if(rank==num_procs-1){
    //     temp_end--;
    // }
    for(i=temp_start;i<temp_end;i++){
        for(j=1;j<ny-1;j++){
            denp[i][j] = den[i][j];
        }
    }

    if(rank>0){
        double temp_send = dimax;
        ierr = MPI_Send( &temp_send, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }

        double temp_recv;
        ierr = MPI_Recv( &temp_recv, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        dimax = temp_recv;
    }
    else if(rank==0){
        for(i=1;i<num_procs;i++){
            double temp_recv;
            ierr = MPI_Recv( &temp_recv, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
            if(ierr!=0){
                printf("Error in Receiving, e_h_field\n");
                exit(0);
            }
            dimax = Max(dimax, temp_recv);
        }
        for(i=1;i<num_procs;i++){
            double temp_send = dimax;
            ierr = MPI_Send( &temp_send, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
            if(ierr!=0){
                printf("Error in Sending, e_h_field\n");
                exit(0);
            }
        }
    }
    
    //c------- time step for fluid equation------------
    dtacmax=0.20*ds*ds/dimax;
    tcycle=1.0/FREQ;
    DTAC=Min(dtacmax, ndifmax*tcycle*ACCEL);
    
    #pragma omp parallel for private(j,da,dac,rec1,fnua,fnui)
    for(i=istart;i<iend;i++){
        #pragma simd
        for(j=1;j<ny-1;j++){
            da=DIFFUSION[i][j];
            dac=da*DTAC/(ds*ds);
            
            rec1=RECOMB*den[i][j]*DTAC;
            
            fnua=0.0;
            fnui=frqio[i][j];
            if(fnui<0.0){
                fnua=-fnui;
                fnui=0.0;
            }

            //c---------- Density equation updates------------------------
            den[i][j]=denp[i][j]*exp(fnui*DTAC)+dac*(denp[i+1][j]+denp[i-1][j]+denp[i][j+1]+denp[i][j-1]-4.0*denp[i][j]);

            den[i][j]=den[i][j]/(1.0+rec1+fnua*DTAC);

            if(den[i][j]<=1e-15)
                den[i][j]=den[i][j]*0.0;

            dnma=Max(dnma,den[i][j]);
        }
    }

    if(rank>0){
        double temp_send = dnma;
        ierr = MPI_Send( &temp_send, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        
        double temp_recv;
        ierr = MPI_Recv( &temp_recv, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        dnma = temp_recv;
    }
    else if(rank==0){
        for(i=1;i<num_procs;i++){
            double temp_recv;
            ierr = MPI_Recv( &temp_recv, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
            if(ierr!=0){
                printf("Error in Receiving, e_h_field\n");
                exit(0);
            }
            dnma = Max(dnma, temp_recv);
        }
        for(i=1;i<num_procs;i++){
            double temp_send = dnma;
            ierr = MPI_Send( &temp_send, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
            if(ierr!=0){
                printf("Error in Sending, e_h_field\n");
                exit(0);
            }
        }
    }

    //c---------------Actual time calculation -------------------------
    TIMD=TIMD+DTAC;
    //c------------------------------------------------------------
    
    omgc2=qe*qe*inv_cmasse/eps0*dnma;
    PARC=omgc2/(pow(OMEG,2)+pow(FNUM,2));

    if(rank==0){
        // ================================================
        for(j=0;j<ny;j++){
            temp_send_arr[j] = den[ulim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            den[ulim[rank]+1][j] = temp_recv_arr[j];
        }
        // ================================================
    } else if(rank<num_procs-1){
        // ================================================
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            den[llim[rank]-1][j] = temp_recv_arr[j];
        }        
        for(j=0;j<ny;j++){
            temp_send_arr[j] = den[llim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        // ================================================
        // ================================================
        for(j=0;j<ny;j++){
            temp_send_arr[j] = den[ulim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            den[ulim[rank]+1][j] = temp_recv_arr[j];
        }
        // ================================================
    } else if(rank==num_procs-1){
        // ================================================
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            den[llim[rank]-1][j] = temp_recv_arr[j];
        }        
        for(j=0;j<ny;j++){
            temp_send_arr[j] = den[llim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        // ================================================
    }
}

void SETUP(){
    int istart, iend, icount;
    istart = llim[rank];
    iend = ulim[rank] + 1;
    icount = ny;
    
    static int k;
    static double ardiy,ardix,dinig;
    static double xd0,yd0,xxi,yyj;
    static double aaa[2012],bbb[2012],ccc[2012],ddd[2012];

    dt=1.0/(double)(nperdt)/FREQ;
    ds=c/(double)(nlamb)/FREQ;

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

    //!================Initial density location =======================

    if(xmid[1]>0) 
        imid[1] = xmid[1]*nx;
    if(ymid[1]>0) 
        jmid[1] = ymid[1]*ny;
    if(xmid[1]==0) 
        imid[1] = nx*0.5;
    if(ymid[1]==0) 
        jmid[1] =ny/2;

    //!=============================================================
    TEMP0=300.0;
    DENG0=PRESSURE/760.0*101300.0/akb/TEMP0;

    radius = nx/5*ds;
    printf("radius: %f\n", radius);
    printf("ds: %f\n", ds);

    int i,j;
    for(i=istart;i<iend;i++){
        for(j=0;j<ny;j++)     
        {
            den[i][j]=0.0;
            ERMS[i][j]=E0/sqrt(2.0);
        }
    }

    xd0=ds*imid[1];
    yd0=ds*jmid[1];

    if(xmid[1]<0) 
        xd0=-xmid[1];
    if(ymid[1]<0) 
        yd0=-ymid[1];

    //!================Initial density, Gaussian, defined =======================
    //!make DEN and DENP =0
    for(i=istart;i<iend;i++){
        xxi=ds*i;
        ardix=0.0;
        if(sgdx0[1]>0)
            ardix=(-pow((xxi-xd0),2))/2.0/sgdx0[1]/sgdx0[1];
        for(j=0;j<ny;j++){    
            yyj=ds*j;
            ardiy=0.0;
            if(sgdy0[1]>0) 
               ardiy=-pow((yyj-yd0),2)/2.0/sgdy0[1]/sgdy0[1];
            dinig=DINI[1]*exp(ardix+ardiy);
            if(dinig<=DINI[1]*1.0*exp(-2)) 
                dinig=0;
            den[i][j]=den[i][j]+dinig;
        }
    }

    if(rank==0){
        // ================================================
        for(j=0;j<ny;j++){
            temp_send_arr[j] = den[ulim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            den[ulim[rank]+1][j] = temp_recv_arr[j];
        }        
        // ================================================
    } else if(rank<num_procs-1){
        // ================================================
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            den[llim[rank]-1][j] = temp_recv_arr[j];
        }        
        for(j=0;j<ny;j++){
            temp_send_arr[j] = den[llim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        // ================================================
        // ================================================
        for(j=0;j<ny;j++){
            temp_send_arr[j] = den[ulim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            den[ulim[rank]+1][j] = temp_recv_arr[j];
        }        
        // ================================================
    } else if(rank==num_procs-1){
        // ================================================
        ierr = MPI_Recv( temp_recv_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        if(ierr!=0){
            printf("Error in Receiving, e_h_field\n");
            exit(0);
        }
        for(j=0;j<ny;j++){
            den[llim[rank]-1][j] = temp_recv_arr[j];
        }        
        for(j=0;j<ny;j++){
            temp_send_arr[j] = den[llim[rank]][j];
        }
        ierr = MPI_Send( temp_send_arr, icount, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
        if(ierr!=0){
            printf("Error in Sending, e_h_field\n");
            exit(0);
        }
        // ================================================
    }

  //c     find relative delay for x, y, Z cell displacement
    dte=dt/eps0;
    dtm=dt/xmu0;
    dteds=dte/ds;
    dtmds=dtm/ds;
  
    c1=(c*dt-ds)/(c*dt+ds);
    c2=2.0*ds/(c*dt+ds);
    c3=(c*dt*c*dt)/(2.0*ds*(c*dt+ds));
}

void ZERO(){
    int i,j;
    int istart, iend, icount;
    istart = llim[rank] - 1;
    iend = ulim[rank] + 2;
    icount = ny;
    if(rank==0)
        istart++;
    else if(rank==num_procs-1)
        iend--;
    
    #pragma omp parallel for private(j)
    for(i=istart;i<iend;i++){   
        #pragma simd
        for(j=0;j<ny;j++) {
            ERMS[i][j]=0.0;
            erms2[i][j]=0.0;
            den[i][j]=0.0;
            denp[i][j]=0.0;
            exi[i][j]=0.0;
            eyi[i][j]=0.0;
            exi1[i][j]=0.0;
            eyi1[i][j]=0.0;
            exs[i][j]=0.0;
            eys[i][j]=0.0;
            hzs[i][j]=0.0;
            vx[i][j]=0.0;
            vy[i][j]=0.0;
            ext[i][j]=0.0;
            eyt[i][j]=0.0;
            exs1[i][j]=0.0;
            eys1[i][j]=0.0;
            frqio[i][j]=0.0;
            DIFFUSION[i][j]=0.0;
        }
    }
    PARC=0.0;
}

void initiate_array(){
    ERMS = malloc(sizeof(double *) * 6200);
    if (ERMS){
        for (i = 0; i < 6200; i++){
            ERMS[i] = malloc(sizeof(double) * 6200);
        }
    }
    erms2 = malloc(sizeof(double *) * 6200);
    if (erms2){
        for (i = 0; i < 6200; i++){
            erms2[i] = malloc(sizeof(double) * 6200);
        }
    }
    exrms = malloc(sizeof(double *) * 6200);
    if (exrms){
        for (i = 0; i < 6200; i++){
            exrms[i] = malloc(sizeof(double) * 6200);
        }
    }
    exrms2 = malloc(sizeof(double *) * 6200);
    if (exrms2){
        for (i = 0; i < 6200; i++){
            exrms2[i] = malloc(sizeof(double) * 6200);
        }
    }
    eyrms = malloc(sizeof(double *) * 6200);
    if (eyrms){
        for (i = 0; i < 6200; i++){
            eyrms[i] = malloc(sizeof(double) * 6200);
        }
    }
    eyrms2 = malloc(sizeof(double *) * 6200);
    if (eyrms2){
        for (i = 0; i < 6200; i++){
            eyrms2[i] = malloc(sizeof(double) * 6200);
        }
    }
    erms = malloc(sizeof(double *) * 6200);
    if (erms){
        for (i = 0; i < 6200; i++){
            erms[i] = malloc(sizeof(double) * 6200);
        }
    }
    den = malloc(sizeof(double *) * 6200);
    if (den){
        for (i = 0; i < 6200; i++){
            den[i] = malloc(sizeof(double) * 6200);
        }
    }
    dent = malloc(sizeof(double *) * 6200);
    if (dent){
        for (i = 0; i < 6200; i++){
            dent[i] = malloc(sizeof(double) * 6200);
        }
    }
    exi = malloc(sizeof(double *) * 6200);
    if (exi){
        for (i = 0; i < 6200; i++){
            exi[i] = malloc(sizeof(double) * 6200);
        }
    }
    eyi = malloc(sizeof(double *) * 6200);
    if (eyi){
        for (i = 0; i < 6200; i++){
            eyi[i] = malloc(sizeof(double) * 6200);
        }
    }
    exi1 = malloc(sizeof(double *) * 6200);
    if (exi1){
        for (i = 0; i < 6200; i++){
            exi1[i] = malloc(sizeof(double) * 6200);
        }
    }
    eyi1 = malloc(sizeof(double *) * 6200);
    if (eyi1){
        for (i = 0; i < 6200; i++){
            eyi1[i] = malloc(sizeof(double) * 6200);
        }
    }
    exs = malloc(sizeof(double *) * 6200);
    if (exs){
        for (i = 0; i < 6200; i++){
            exs[i] = malloc(sizeof(double) * 6200);
        }
    }
    eys = malloc(sizeof(double *) * 6200);
    if (eys){
        for (i = 0; i < 6200; i++){
            eys[i] = malloc(sizeof(double) * 6200);
        }
    }
    hzi = malloc(sizeof(double *) * 6200);
    if (hzi){
        for (i = 0; i < 6200; i++){
            hzi[i] = malloc(sizeof(double) * 6200);
        }
    }
    hzs = malloc(sizeof(double *) * 6200);
    if (hzs){
        for (i = 0; i < 6200; i++){
            hzs[i] = malloc(sizeof(double) * 6200);
        }
    }
    vx = malloc(sizeof(double *) * 6200);
    if (vx){
        for (i = 0; i < 6200; i++){
            vx[i] = malloc(sizeof(double) * 6200);
        }
    }
    vy = malloc(sizeof(double *) * 6200);
    if (vy){
        for (i = 0; i < 6200; i++){
            vy[i] = malloc(sizeof(double) * 6200);
        }
    }
    ext = malloc(sizeof(double *) * 6200);
    if (ext){
        for (i = 0; i < 6200; i++){
            ext[i] = malloc(sizeof(double) * 6200);
        }
    }
    eyt = malloc(sizeof(double *) * 6200);
    if (eyt){
        for (i = 0; i < 6200; i++){
            eyt[i] = malloc(sizeof(double) * 6200);
        }
    }
    exs1 = malloc(sizeof(double *) * 6200);
    if (exs1){
        for (i = 0; i < 6200; i++){
            exs1[i] = malloc(sizeof(double) * 6200);
        }
    }
    exs2 = malloc(sizeof(double *) * 6200);
    if (exs2){
        for (i = 0; i < 6200; i++){
            exs2[i] = malloc(sizeof(double) * 6200);
        }
    }
    eys1 = malloc(sizeof(double *) * 6200);
    if (eys1){
        for (i = 0; i < 6200; i++){
            eys1[i] = malloc(sizeof(double) * 6200);
        }
    }
    eys2 = malloc(sizeof(double *) * 6200);
    if (eys2){
        for (i = 0; i < 6200; i++){
            eys2[i] = malloc(sizeof(double) * 6200);
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
    temp_send_arr = malloc(sizeof(double) * 6200);
    temp_recv_arr = malloc(sizeof(double) * 6200);
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
    free(temp_send_arr);
    free(temp_recv_arr);
}
