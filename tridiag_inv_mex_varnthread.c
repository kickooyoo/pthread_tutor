/* tridiag_inv_mex.c
 Matlab mex file for tridiagonal solver 2013-11-09
 T * x = y
 assumes real tridiagonal matrix T, rhs y is formulated as NxM matrix, computes M columns of x in parallel
*/

#include "mex.h"
#include "defs-env.h"
#include "jf,mex,def.h"
//#include "def/def,mexarg.h" // in above
#include "jf,thread1.h"
#include "pthread.h"

#define Usage "usage error. see above"
#define VERBOSE false

//pthread_mutex_t mutexout; // global var for locking

static void tridiag_inv_varnthread_mex_help(void)
{
	printf("\n\
	Usage for tridiag_inv_mex: \n\
	output = tridiag_inv_mex(subdiag,diagvals,supdiag,argum) \n\
	\n\
	subdiag has length n-1	\n\
	expect tridiag matrix to be real \n\
	\n");
}


struct thread_data
{
	int thread_id;
	int block_size;
    int num_blocks_for_me;
	float *subdiag_ptr;
	float *diagvals_ptr;
	float *supdiag_ptr;
	float *rhsr_ptr;
    float *rhsi_ptr;
	float *outr_ptr;
    float *outi_ptr;
};

// tridiag_inv()
// tridiag solver over one block
static sof tridiag_inv(float *a, float *b, float *c, float *d, int N, float *x)
{
    float *new_c;
    float *new_d;
    
	new_c = (float *) calloc (N - 1, sizeof(float));
	new_d = (float *) calloc (N, sizeof(float));

	*new_c = *c / *b; //new_c[0] = c[0]/b[0];
	*new_d = *d / *b; //new_d[0] = d[0]/b[0];

	int ii;
	float a_prev, new_c_prev;
	for (ii = 1; ii <= N-2; ii++) {
		a_prev = a[ii - 1];
		new_c_prev = new_c[ii - 1];
		new_c[ii] = c[ii] / (b[ii] - new_c_prev * a_prev);
		new_d[ii] = (d[ii] - new_d[ii - 1] * a_prev) / (b[ii] - new_c_prev * a_prev);
	}
	new_d[N - 1] = (d[N - 1] - new_d[N - 2] * (a[N - 2])) / (b[N - 1] - new_c[N - 2] * (a[N - 2]));

	x[N - 1] = new_d[N - 1];
	for (ii = N-2; ii >= 0; ii--) {
        x[ii] = new_d[ii] - new_c[ii] * x[ii + 1];
	}
#if VERBOSE
	printf("task id: %d \n", taskid);
#endif

//	pthread_mutex_lock(&mutexout);
	
//	pthread_mutex_unlock(&mutexout);
	free(new_c);
	free(new_d);

	Ok
}

// tridiag_inv_loop_thr()
// each thread loops over tridiag_inv
static sof tridiag_inv_loop_thr(void *threadarg, cint tid, cint nthread)
{
    
    // todo: pass in just one set of data and fastforward according to cum_blocks and tid
    struct thread_data *my_data;
    int taskid;
    int num_runs;
    int N;
    float *a;
    float *b;
    float *c;
    float *dr;
    float *xr;
    float *di;
    float *xi;
    my_data = (struct thread_data *) threadarg;
    taskid = my_data -> thread_id;
    N = my_data -> block_size;
    a = my_data -> subdiag_ptr;
    b = my_data -> diagvals_ptr;
    c = my_data -> supdiag_ptr;
    dr = my_data -> rhsr_ptr;
    xr = my_data -> outr_ptr;
    di = my_data -> rhsi_ptr;
    xi = my_data -> outi_ptr;
    num_runs = my_data -> num_blocks_for_me;
#if VERBOSE
    printf("inside looper func with num runs %d \n", num_runs);
#endif
    for (int ii = 0; ii < num_runs; ii++) {
#if VERBOSE
        printf(" looping over index %d /%d \n", ii, num_runs);
        
        printf("address of rhs: %d, incr by %d,  new addr: %d, float size: %d \n", dr, N, dr+N, sizeof(float));
        printf("rhs vals: \n");
        for (int jj = 0; jj < N; jj++) {
            printf("rhs[%d] = %d \n", jj, dr[jj]);
        }
#endif
        tridiag_inv(a, b, c, dr, N, xr);
        dr += N;
        xr += N;
        if (xi != NULL) {
            tridiag_inv(a, b, c, di, N, xi);
            di += N;
            xi += N;
        }
    }
    
    pthread_exit(NULL);
    Ok
}

/*
 static sof check_types_and_sizes(
Const mxArray *subdiag, 
Const mxArray *diagvals, 
Const mxArray *supdiag, 
Const mxArray *rhs, 
Const mxArray *block_size_ptr)
{
	int pass = 1;
	if (!mxIsDouble(subdiag)) {
		printf("subdiag is not double type \n");
		pass = 0;
	}
	if (!mxIsDouble(supdiag)) {
		printf("supdiag is not double type \n");
		pass = 0;
	}
	if (!mxIsDouble(diagvals)) {
		printf("diagvals is not double type \n");
		pass = 0;
	}
	if (!mxIsDouble(rhs)) {
		printf("rhs is not double type \n");
		pass = 0;
	}
	// todo: check lengths using mxGetM
	if (pass) {
		printf(" all are double types \n");
	}
	printf("\n");
    Ok
} */


// currently assume all inputs real
// wrapper function for thread
static sof tridiag_inv_varnthread_mex_thr(
float *subdiag_ptr, float *diagvals_ptr, float *supdiag_ptr, float *rhs_real_ptr, float *rhs_imag_ptr, mwSize block_size, mwSize nblocks, float *out_real_ptr, float *out_imag_ptr, int ncores)
{
	//tridiag_inv_worker_args args;

	//int ii;
	int big_N; // total number of entries, N*nblocks
	int rc;
	long t;
    int block_ndx;
    int *blocks_per_thread;
    int *cum_blocks;
    int remainder_blocks;
    struct thread_data *thread_data_array;
    int nthread;
    
    nthread = ncores; //jf_thread1_ncore(-1); // trouble!
    
    printf(" ncores detected : %d, mandated: %d \n", nthread, ncores);
    
    Mem0(blocks_per_thread, nthread);
    Mem0(cum_blocks, nthread + 1);
    Mem0(thread_data_array, nthread);

    pthread_attr_t attr;
	pthread_t threads[nthread];
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    remainder_blocks = nblocks - (nblocks/nthread)*nthread;
#if VERBOSE
    printf("NUM_THREADS: %d, nblocks: %d \n", nthread, nblocks);
    printf("remainder_blocks: %d \n", remainder_blocks);
#endif
    
    cum_blocks[0] = 0;
    for (int th_ndx = 0; th_ndx < nthread; th_ndx++) {
        blocks_per_thread[th_ndx] = nblocks/nthread;
        if (th_ndx < remainder_blocks) {
            blocks_per_thread[th_ndx]++;
        }
        cum_blocks[th_ndx + 1] = cum_blocks[th_ndx] + blocks_per_thread[th_ndx];
#if VERBOSE
        printf("cum_blocks_per_thread for %d : %d \n", th_ndx+1, cum_blocks[th_ndx+1]);
#endif
    }
#if VERBOSE
    
    printf("rhs vals: \n");
    for (int jj = 0; jj < block_size; jj++) {
        printf("rhs[%d] = %d \n", jj, *(rhs_real_ptr + jj) );
    }
#endif
    
    
    
    for (int th_id = 0; th_id < nthread; th_id++) {
        thread_data_array[th_id].thread_id = th_id;
        thread_data_array[th_id].block_size = block_size;
        thread_data_array[th_id].subdiag_ptr = subdiag_ptr;
        thread_data_array[th_id].diagvals_ptr = diagvals_ptr;
        thread_data_array[th_id].supdiag_ptr = supdiag_ptr;
        thread_data_array[th_id].rhsr_ptr = rhs_real_ptr + cum_blocks[th_id] * block_size;
        thread_data_array[th_id].outr_ptr = out_real_ptr + cum_blocks[th_id] * block_size;
        thread_data_array[th_id].rhsi_ptr = rhs_imag_ptr + cum_blocks[th_id] * block_size;
        thread_data_array[th_id].outi_ptr = out_imag_ptr + cum_blocks[th_id] * block_size;
        thread_data_array[th_id].num_blocks_for_me = blocks_per_thread[th_id];
#if VERBOSE
        printf(" th_id: %d, rhs ptr: %d, curr ptr: %d \n \n", th_id, rhs_real_ptr, thread_data_array[th_id].rhsr_ptr);
#endif
        //rc = pthread_create(&threads[th_id], &attr, (void *) &tridiag_inv_loop_thr, (void *) &(thread_data_array[th_id]));
    }
    
    Free0(blocks_per_thread);
    Free0(cum_blocks);
    Free0(thread_data_array);
    

    Ok
}


// intermediate GateWay routine 
static sof tridiag_inv_varnthread_mex_gw(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
    float *sub;              /* input subdiagonal */
    float *diag;               /* 1xN input diagonal */
    float *sup;
    float *rhs;
    float *rhs_imag;
    int *ncores;

    size_t N;                   /* size of tridiag matrix */
    size_t M;                   /* numcols of rhs matrix */
    float *x_real;                  /* output */
    float *x_imag;                  /* output */
    
	if (nrhs != 5 ) { // hard coding :(
        printf("Incorrect number of inputs. \n");
		tridiag_inv_varnthread_mex_help();
		//Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	// todo: check lengths of vectors
	
    // todo: check type of inputs
    if ( (mxIsSingle(prhs[0]) == 0) || (mxIsSingle(prhs[1]) == 0) || (mxIsSingle(prhs[2]) == 0) ) {
        tridiag_inv_varnthread_mex_help();
        printf("Inputs need to be single precision. \n ");
        //Call(mxu_arg, (nrhs, prhs))
        Fail(Usage)
    }
    
    
    // todo: check complexity of inputs
    
    
    sub = (float *) mxGetData(prhs[0]);
    diag = (float *)  mxGetData(prhs[1]);
    sup =  (float *) mxGetData(prhs[2]);
    rhs = (float *) mxGetData(prhs[3]);
    ncores = (int *) mxGetData(prhs[4]);
    N = mxGetM(prhs[1]); //remember col vec!
    M = mxGetN(prhs[3]);

#if VERBOSE
    printf("M: %d, N: %d \n", M, N);

    printf(" address of rhs in gw: %d \n", rhs);
    printf("rhs vals: \n");
    for (int jj = 0; jj < N; jj++) {
        printf("rhs[%d] = %f \n", jj, *(rhs + jj) );
    }
#endif
    
    if (mxIsComplex(prhs[3])) {
        //printf("rhs is complex \n");
        rhs_imag = (float *) mxGetImagData(prhs[3]);
        //plhs[0] = mxCreateDoubleMatrix((mwSize)N, (mwSize) M, mxCOMPLEX);
        plhs[0] = mxCreateNumericMatrix((int) N * (int) M, 1, mxSINGLE_CLASS, mxCOMPLEX);
    } else {
        rhs_imag = NULL;
        //plhs[0] = mxCreateDoubleMatrix((mwSize)N, (mwSize) M, mxREAL);
        plhs[0] = mxCreateNumericMatrix((int) N *  (int) M, 1, mxSINGLE_CLASS, mxREAL);
    }
    x_real = (float *)mxGetData(plhs[0]);
    x_imag = (float *)mxGetImagData(plhs[0]); // NULL when rhs is real

    
   // check_types_and_sizes(sub, diag, sup, rhs, N, M);
    
    
    if ((rhs_imag == NULL) ^ (x_imag == NULL)) {
        printf("problem: only one NULL vec, should be both or neither");
    }
    
	Call(tridiag_inv_varnthread_mex_thr, (sub, diag, sup, rhs, rhs_imag, N, M, x_real, x_imag, *ncores));
    //tridiag_inv_varnthread_mex_thr(sub, diag, sup, rhs, rhs_imag, N, M, x_real, x_imag, *ncores);
	
	Ok
}

// gateway routine 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		tridiag_inv_varnthread_mex_help();
		return;
	}
	if (!tridiag_inv_varnthread_mex_gw(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("tridiag_inv_mex");
}



