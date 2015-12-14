/* sense_example_mex.c, sense done with a mex file 
 * Copyright 10-01-2013, Matthew Muckley, University of Michigan*/

#include <pthread.h> /* for threading */
#include <unistd.h> /* for sleep */
#include <fftw3.h>
#include "mex.h"
#include "matrix.h"

/* a poorly-written ifftshift function
 * warn: not tested for odd N */
void threedifftshift(double *x, int rank, fftw_iodim *dim)
{
    int id1;
    int id2;
    int id3;
    int demcount;
    int demmax;
    double tmp;
    double tmp2;
    double *xp;

    for (id1 = 0; id1 < rank; ++id1) {
        if (id1 == 0) {
            demmax = dim[1].n*dim[2].n;
        } else if (id1 == 1) {
            demmax = dim[0].n*dim[2].n;
        } else if (id1 == 2) {
            demmax = dim[0].n*dim[1].n;
        }
        for (id2 = 0; id2 < demmax; ++id2) {
            if (id1 == 0) {
                xp = x + (id2 % dim[1].n)*dim[0].n + (id2 / dim[1].n)*dim[1].n*dim[0].n;
            } else if (id1 == 1) {
                xp = x + (id2 % dim[0].n) + (id2 / dim[0].n)*dim[1].n*dim[0].n;
            } else if (id1 == 2) {
                xp = x + (id2 % dim[0].n) + (id2 / dim[0].n)*dim[0].n;
            }
            if (dim[id1].n % 2 == 0) {
                for (id3 = 0; id3 < dim[id1].n/2; ++id3) {
                    tmp = *xp;
                    *xp = *(xp+dim[id1].n/2*dim[id1].is);
                    *(xp+dim[id1].n/2*dim[id1].is) = tmp;
                    xp += dim[id1].is;
                }
            } else {
                tmp2 = *(xp+dim[id1].n/2*dim[id1].is);
                xp += (dim[id1].n/2+1)*dim[id1].is;
                for (id3 = 0; id3 < dim[id1].n/2; ++id3) {
                    tmp = *xp;
                    *xp = *(xp-(dim[id1].n/2+1)*dim[id1].is);
                    *(xp-(dim[id1].n/2)*dim[id1].is) = tmp;
                    xp += dim[id1].is;
                }
            }
        }
    }
}

/* a poorly-written fftshift function
 * warn: not tested for odd N */
void threedfftshift(double *x, int rank, fftw_iodim *dim)
{
    int id1;
    int id2;
    int id3;
    int demmax;
    double tmp;
    double tmp2;
    double *xp;

    for (id1 = 0; id1 < rank; ++id1) {
        if (id1 == 0) {
            demmax = dim[1].n*dim[2].n;
        } else if (id1 == 1) {
            demmax = dim[0].n*dim[2].n;
        } else if (id1 == 2) {
            demmax = dim[0].n*dim[1].n;
        }
        for (id2 = 0; id2 < demmax; ++id2) {
            if (id1 == 0) {
                xp = x + (id2 % dim[1].n)*dim[0].n + (id2 / dim[1].n)*dim[1].n*dim[0].n;
            } else if (id1 == 1) {
                xp = x + (id2 % dim[0].n) + (id2 / dim[0].n)*dim[1].n*dim[0].n;
            } else if (id1 == 2) {
                xp = x + (id2 % dim[0].n) + (id2 / dim[0].n)*dim[0].n;
            }
            if (dim[id1].n % 2 == 0) {
                for (id3 = 0; id3 < dim[id1].n/2; ++id3) {
                    tmp = *xp;
                    *xp = *(xp+dim[id1].n/2*dim[id1].is);
                    *(xp+dim[id1].n/2*dim[id1].is) = tmp;
                    xp += dim[id1].is;
                }
            } else {
                tmp2 = *(xp+dim[id1].n/2*dim[id1].is);
                for (id3 = 0; id3 < dim[id1].n/2; ++id3) {
                    tmp = *xp;
                    *xp = *(xp+(dim[id1].n/2+1)*dim[id1].is);
                    *(xp+(dim[id1].n/2)*dim[id1].is) = tmp;
                    xp += dim[id1].is;
                }
                *xp = tmp2;
            }
        }
    }
}

/* sense_args data struct */
typedef struct
{
    double *xreal;
    double *ximag;
    double *smapreal;
    double *smapimag;
    int Nx;
    int Ny;
    int Nz;
    int Ncoil;
    double *yreal;
    double *yimag;
    int* ystart;
    int* ycount;
    fftw_plan *plan;
    fftw_iodim dims[3];
    fftw_iodim howmany_dims[1];
    int nthread;
    int chat;
} sense_args;

/* tric: build a structure with a pointer to data and a
 * thread id number */
typedef struct
{
    sense_args* structpoint;
    int tid;
} ultimate_pointer;

/* single thread workload */
void *sense_thread(void *p)
{
    /* initialize the variables */
    ultimate_pointer* megapoint = (ultimate_pointer*) p;
    sense_args* args = megapoint->structpoint;
    double *xpointreal;
    double *xpointimag;
    double *ypointreal;
    double *ypointimag;
    double *smappointreal;
    double *smappointimag;
    unsigned long Nx = (unsigned long)args->Nx;
    unsigned long Ny = (unsigned long)args->Ny;
    unsigned long Nz = (unsigned long)args->Nz;
    unsigned long Ncoil = args->Ncoil;
    int tid = megapoint->tid;
    int id1;
    int id2;

    /* get the pointers to the right locations
     * unsigned longs probably unnecessary; mostly paranoia */
    ypointreal = args->yreal + (unsigned long)args->ystart[tid]*Nx*Ny*Nz;
    ypointimag = args->yimag + (unsigned long)args->ystart[tid]*Nx*Ny*Nz;
    xpointreal = args->xreal;
    xpointimag = args->ximag;
    smappointreal = args->smapreal + (unsigned long)args->ystart[tid]*Nx*Ny*Nz;
    smappointimag = args->smapimag + (unsigned long)args->ystart[tid]*Nx*Ny*Nz;

    if (args->chat)
        mexPrintf("thread %d started\n", tid);

    /* basis expansion for each sense coil */
    for (id1 = 0; (unsigned long)id1 < args->ycount[tid]; ++id1) {
        xpointreal = args->xreal;
        xpointimag = args->ximag;
        for (id2 = 0; (unsigned long)id2 < Nx*Ny*Nz; ++id2) {
            *ypointreal++ = (*xpointreal)*(*smappointreal) - 
                (*xpointimag)*(*smappointimag);
            *ypointimag++ = (*xpointimag++)*(*smappointreal++) + 
                (*xpointreal++)*(*smappointimag++);
        }
    }

    /* put the pointers back in the right spot */
    ypointreal = args->yreal + (unsigned long)args->ystart[tid]*Nx*Ny*Nz;
    ypointimag = args->yimag + (unsigned long)args->ystart[tid]*Nx*Ny*Nz;

    if (args->chat)
        printf("thread %d is starting ffts\n", tid);
   
    /* do the ifftshift */
    for (id1=0; id1 < args->ycount[tid]; ++id1) { 
        threedifftshift(ypointreal + (unsigned long)id1*Nx*Ny*Nz, (int)3, args->dims);
        threedifftshift(ypointimag + (unsigned long)id1*Nx*Ny*Nz, (int)3, args->dims);
    }
   
    if (args->chat)
        printf("thread %d done with ifftshifts\n", tid);

    /* run the ffts */
    fftw_execute(args->plan[tid]);

    if (args->chat)
        printf("thread %d done with ffts\n", tid);
    
    /* do the fftshift */
    for (id1=0; id1 < args->ycount[tid]; ++id1) {
        threedfftshift(ypointreal + (unsigned long)id1*Nx*Ny*Nz, (int)3, args->dims);
        threedfftshift(ypointimag + (unsigned long)id1*Nx*Ny*Nz, (int)3, args->dims);
    }

    if (args->chat)
        printf("thread %d is done with fftshifts\n", tid);

    pthread_exit(NULL);
}

/* entry point function */
void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    /* declare all the variables */
    sense_args args;
    int rc;
    int y_per_thread;
    int y_leftover;
    int y_claimed;
    int this_claim;
    int id1;
    double* ypointreal;
    double* ypointimag;
    ultimate_pointer *megapoint;
    
    /* check that we have the right inputs, read in the data */
    if (nrhs != 8) {
        mexErrMsgTxt("not correct input arguments");
    } else {
        if (mxIsDouble(prhs[0]) && mxIsComplex(prhs[0])) {
            args.xreal = mxGetPr(prhs[0]);
            args.ximag = mxGetPi(prhs[0]);
        } else {
            mexErrMsgTxt("first entry not double or complex\n");
        }
        if (mxIsDouble(prhs[1]) && mxIsComplex(prhs[1])) {
            args.smapreal = mxGetPr(prhs[1]);
            args.smapimag = mxGetPi(prhs[1]);
        } else {
            mexErrMsgTxt("second entry not double or complex\n");
        }
        if (mxIsInt32(prhs[2]) && mxIsInt32(prhs[3]) && 
                mxIsInt32(prhs[4]) && mxIsInt32(prhs[5]) && 
                mxIsInt32(prhs[6]) && mxIsInt32(prhs[7])) {
            args.Nx = mxGetScalar(prhs[2]);
            args.Ny = mxGetScalar(prhs[3]);
            args.Nz = mxGetScalar(prhs[4]);
            args.Ncoil = mxGetScalar(prhs[5]);
            args.nthread = mxGetScalar(prhs[6]);
            args.chat = mxGetScalar(prhs[7]);
        } else {
            mexErrMsgTxt("one of last six entries not an int\n");
        }
    }

    /* start arrays for dividing up work */
    args.ystart = (int*)calloc(args.nthread, sizeof(int));
    args.ycount = (int*)calloc(args.nthread, sizeof(int));

    y_per_thread = args.Ncoil / args.nthread;
    y_leftover = args.Ncoil % args.nthread;
    y_claimed = 0;

    for (id1 = 0; id1 < args.nthread; ++id1) {
        args.ystart[id1] = y_claimed;
        this_claim = y_per_thread + (id1 < y_leftover ? 1 : 0);
        args.ycount[id1] = this_claim;
        y_claimed += this_claim;
        if (args.chat)
            printf("sense_example_mex: claiming %d loops for thread %d\n", this_claim, id1);
    }

    /* output data matrix */
    plhs[0] = mxCreateNumericMatrix(args.Ncoil*args.Nx*args.Ny*args.Nz, 1, mxDOUBLE_CLASS, mxCOMPLEX);

    args.yreal = mxGetPr(plhs[0]);
    args.yimag = mxGetPi(plhs[0]);
    
    /* trick: structure to pass data structure pointer
     * and thread id number */
    megapoint = (ultimate_pointer*)calloc(args.nthread, sizeof(ultimate_pointer));

    /* array of fftw3 plans, could possibly switch up */
    args.plan = (fftw_plan*)calloc(args.nthread, sizeof(fftw_plan));

    /* dims for fftw3 */
    args.dims[0].n = (int)args.Nx;
    args.dims[0].is = 1;
    args.dims[0].os = 1;
    args.dims[1].n = (int)args.Ny;
    args.dims[1].is = (int)args.Nx;
    args.dims[1].os = (int)args.Nx;
    args.dims[2].n = (int)args.Nz;
    args.dims[2].is = (int)(args.Ny*args.Nx);
    args.dims[2].os = (int)(args.Ny*args.Nx);

    /* build all necessary thread-specific data structures */
    for (id1 = 0; id1 < args.nthread; ++id1) {
        ypointreal = args.yreal + args.ystart[id1]*args.Nx*args.Ny*args.Nz;
        ypointimag = args.yimag + args.ystart[id1]*args.Nx*args.Ny*args.Nz;
        args.howmany_dims[0].n = args.ycount[id1];
        args.howmany_dims[0].is = (int)(args.Nx*args.Ny*args.Nz);
        args.howmany_dims[0].os = (int)(args.Nx*args.Ny*args.Nz);
        args.plan[id1] = fftw_plan_guru_split_dft((int)3, args.dims, (int)1, args.howmany_dims,
            ypointreal, ypointimag, ypointreal, ypointimag, FFTW_ESTIMATE);
        megapoint[id1].structpoint = &args;
        megapoint[id1].tid = id1;
    }

    /* start the threads up */
    pthread_t threads[args.nthread];
    for (id1 = 0; id1 < args.nthread ; ++id1) {
        if (args.chat)
            mexPrintf("creating thread %d\n", id1);
        rc = pthread_create(&threads[id1], NULL, sense_thread, &(megapoint[id1]));
    }
   
    /* wait for thread termination, destroy fftw plan */
    for (id1 = 0; id1< args.nthread; ++id1) {
        pthread_join(threads[id1], NULL);
        fftw_destroy_plan(args.plan[id1]);
    }
    
    free(args.plan);
    free(megapoint);
    free(args.ystart);
    free(args.ycount);
    
    return;
}
