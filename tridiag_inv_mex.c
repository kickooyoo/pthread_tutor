/* mex file for tridiagonal solver 2013-11-09 
 assuming block diagonal with identical, real tridiagonal blocks
 so rhs is formulated as NxM matrix */


#include "mex.h"
#include "defs-env.h"
#include "jf,mex,def.h"
//#include "def/def,mexarg.h" // in above
#include "pthread.h"

//#if !defined(Need_tridiag_inv_mex_gateway)
//#include "tridiag_inv,def.h"
//#endif

#define Usage "usage error. see above"
#define NUM_THREADS 2 // number of cores // 4 for iv1, 2 for vega

//pthread_mutex_t mutexout; // global var for locking

//#if defined(Mmex)

static void tridiag_inv_mex_help(void)
{
	printf("\n\
	Usage for tridiag_inv_mex: \n\
	output = tridiag_inv_mex(subdiag,diagvals,supdiag,argum) \n\
	\n\
	subdiag has length n-1	\n\
	expect tridiag matrix to be real \n\
	\n");
}

/*typedef struct 
{
	float *yy;
	int Ny; // etc??
} tridiag_inv_worker_args; */


struct thread_data
{
	int thread_id;
	int block_size;
    int num_blocks_for_me;
	float *subdiag_ptr;
	float *diagvals_ptr;
	float *supdiag_ptr;
	float *rhs_ptr;
	float *out_ptr;
};

struct thread_data thread_data_array[NUM_THREADS];

static sof tridiag_inv(float *a, float *b, float *c, float *d, int N, float *x)
{
	//tridiag_inv_worker_args args;
/*	int taskid;
	int N;
	float *a;
	float *b;
	float *c;
    float *d;
    float *x;

	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
	taskid = my_data -> thread_id;
	N = my_data -> block_size;
	a = my_data -> subdiag_ptr;
	b = my_data -> diagvals_ptr;
	c = my_data -> supdiag_ptr;
	d = my_data -> rhs_ptr;
	x = my_data -> out_ptr; */
    
    float *new_c;
    float *new_d;
    
	// need to malloc for x_real, x_imag, and then free them ??
	new_c = (float *) calloc (N - 1, sizeof(float));
	new_d = (float *) calloc (N, sizeof(float));

	*new_c = *c / *b; //new_c[0] = c[0]/b[0];
	*new_d = *d / *b; //new_d[0] = d[0]/b[0];

	int ii;
	float a_prev, new_c_prev;
	for (ii = 1; ii <= N-2; ii++) {
		/*
		a_prev = *(a + ii - 1);
		new_c_prev = *(new_c + ii - 1);
		*(new_c + ii) = *(c + ii) / (*(b + ii) - new_c_prev * a_prev);
		*(new_d + ii) = (*(d + ii) - *(new_d + ii - 1) * a_prev) / (*(b + ii) - new_c_prev * a_prev);
		*/
		a_prev = a[ii - 1];
		new_c_prev = new_c[ii - 1];
		new_c[ii] = c[ii] / (b[ii] - new_c_prev * a_prev);
		new_d[ii] = (d[ii] - new_d[ii - 1] * a_prev) / (b[ii] - new_c_prev * a_prev);
	}
	//*(new_d + N - 1) = (*(d + N - 1) - *(new_d + N - 2) * (*(a + N - 1))) / (*(b + N - 1) - *(new_c + N - 2) * (*(a + N - 2)));
	new_d[N - 1] = (d[N - 1] - new_d[N - 2] * (a[N - 2])) / (b[N - 1] - new_c[N - 2] * (a[N - 2]));

	//*(x + N - 1) = *(new_d + N - 1);
	x[N - 1] = new_d[N - 1];
	for (ii = N-2; ii >= 0; ii--) {
		//*(x + ii) = *(new_d + ii) - *(new_c + ii) * (*(x + ii + 1));
		x[ii] = new_d[ii] - new_c[ii] * x[ii + 1];
	}

//	printf("task id: %d \n", taskid);

//	pthread_mutex_lock(&mutexout);
	
//	pthread_mutex_unlock(&mutexout);
	free(new_c);
	free(new_d);

	Ok
}


static sof tridiag_inv_loop_thr(void *threadarg) // for loop over tridiag_inv_thr
{
    struct thread_data *my_data;
    int taskid;
    int num_runs;
    int N;
    float *a;
    float *b;
    float *c;
    float *d;
    float *x;
//    float *new_c;
//    float *new_d;
    my_data = (struct thread_data *) threadarg;
    taskid = my_data -> thread_id;
    N = my_data -> block_size;
    a = my_data -> subdiag_ptr;
    b = my_data -> diagvals_ptr;
    c = my_data -> supdiag_ptr;
    d = my_data -> rhs_ptr;
    x = my_data -> out_ptr;
    num_runs = my_data -> num_blocks_for_me;
    
    printf("inside looper func with num runs %d \n", num_runs);
    for (int ii = 0; ii < num_runs; ii++) {
        printf(" looping over index %d /%d \n", ii, num_runs);
        
        printf("address of rhs: %d, incr by %d,  new addr: %d, float size: %d \n", d, N, d+N, sizeof(float));
        tridiag_inv(a, b, c, d, N, x);
        
                // increment to block, same tridiag for all so no inc a, b, c
        d += N;
        x += N;
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
static sof tridiag_inv_mex_thr(
float *subdiag_ptr, float *diagvals_ptr, float *supdiag_ptr, float *rhs_real_ptr, float *rhs_imag_ptr, mwSize block_size, mwSize nblocks, float *out_real_ptr, float *out_imag_ptr)
{
	//tridiag_inv_worker_args args;

	//int ii;
	int big_N; // total number of entries, N*nblocks
	int rc;
	long t;
    int block_ndx;
    int blocks_per_thread;
    
    pthread_attr_t attr;
	
//	printf("block size : %d \n", block_size);
//    printf("nblocks : %d, NUM_THREADS: %d, nblocks/NUM_THREADS: %d \n", nblocks, NUM_THREADS, (nblocks/ NUM_THREADS));
	pthread_t threads[NUM_THREADS];

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    blocks_per_thread = nblocks/NUM_THREADS;
    
    printf(" blocks per thread: %d \n", blocks_per_thread);
    
    
    

    
    //    printf("\n");
    // do all real values first
    //for (int th_rep = 0; th_rep <= nblocks/NUM_THREADS; th_rep++) {
    for (int th_id = 0; th_id < NUM_THREADS; th_id++) {
        //block_ndx = th_rep * NUM_THREADS + th_id;
        //            printf("th_rep: %d, th_id: %d, block_ndx: %d \n", th_rep, th_id, block_ndx);
        //if (block_ndx <= nblocks - 1) {
        thread_data_array[th_id].thread_id = th_id;
        thread_data_array[th_id].block_size = block_size;
        thread_data_array[th_id].subdiag_ptr = subdiag_ptr;// + th_id * block_size;
        thread_data_array[th_id].diagvals_ptr = diagvals_ptr;// + th_id * block_size;
        thread_data_array[th_id].supdiag_ptr = supdiag_ptr;// + th_id * block_size;
        thread_data_array[th_id].rhs_ptr = rhs_real_ptr + th_id * blocks_per_thread * block_size;
        thread_data_array[th_id].out_ptr = out_real_ptr + th_id * blocks_per_thread * block_size;
        thread_data_array[th_id].num_blocks_for_me = blocks_per_thread;
        
        printf(" th_id: %d, block size: %d, incr: %d , nblocks: %d \n", th_id, block_size, blocks_per_thread * block_size, nblocks);
        rc = pthread_create(&threads[th_id], &attr, (void *) &tridiag_inv_loop_thr, (void *) &(thread_data_array[th_id]));
        //}
    }
    Ok
        /*
    //        printf("\n");
    for (int th_id = 0; th_id < NUM_THREADS; th_id++) {
        //  block_ndx = th_rep * NUM_THREADS + th_id;
        //  if (block_ndx <= nblocks - 1) {
        rc = pthread_join(threads[th_id], NULL);
        if (rc) {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            exit(-1);
        }
        //                printf("done with thread ndx %d ", block_ndx);
        // }
    }
    //        printf("\n");
    //        printf("done with thread blocks %d - %d \n", th_rep*(nblocks/NUM_THREADS), block_ndx);
//}

    if ((rhs_imag_ptr != NULL) && (out_imag_ptr != NULL)) {
     
    // do all complex values next
    for (int th_rep = 0; th_rep <= nblocks/NUM_THREADS; th_rep++) {
        for (int th_id = 0; th_id < NUM_THREADS; th_id++) {
            block_ndx = th_rep * NUM_THREADS + th_id;
            //            printf("th_rep: %d, th_id: %d, block_ndx: %d \n", th_rep, th_id, block_ndx);
            if (block_ndx <= nblocks - 1) {
                thread_data_array[th_id].thread_id = th_id;
                thread_data_array[th_id].block_size = block_size;
                thread_data_array[th_id].subdiag_ptr = subdiag_ptr;// + th_id * block_size;
                thread_data_array[th_id].diagvals_ptr = diagvals_ptr;// + th_id * block_size;
                thread_data_array[th_id].supdiag_ptr = supdiag_ptr;// + th_id * block_size;
                thread_data_array[th_id].rhs_ptr = rhs_imag_ptr + block_ndx * block_size;
                thread_data_array[th_id].out_ptr = out_imag_ptr + block_ndx * block_size;
                rc = pthread_create(&threads[th_id], &attr, (void *) &tridiag_inv_loop_thr, (void *) &(thread_data_array[th_id]));
            }
        }
        //        printf("\n");
        for (int th_id = 0; th_id < NUM_THREADS; th_id++) {
            block_ndx = th_rep * NUM_THREADS + th_id;
            if (block_ndx <= nblocks - 1) {
                rc = pthread_join(threads[th_id], NULL);
                if (rc) {
                    printf("ERROR; return code from pthread_join() is %d\n", rc);
                    exit(-1);
                }
                //                printf("done with thread ndx %d ", block_ndx);
            }
        }
        //        printf("\n");
        //        printf("done with thread blocks %d - %d \n", th_rep*(nblocks/NUM_THREADS), block_ndx);
    }
    }
    
    
	Ok */
}


// intermediate GateWay routine 
static sof tridiag_inv_mex_gw(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
    float *sub;              /* input subdiagonal */
    float *diag;               /* 1xN input diagonal */
    float *sup;
    float *rhs;
    float *rhs_imag;

    size_t N;                   /* size of tridiag matrix */
    size_t M;                   /* numcols of rhs matrix */
    float *x_real;                  /* output */
    float *x_imag;                  /* output */
    
	if (nrhs != 4 ) { // hard coding :(
		tridiag_inv_mex_help();
		//Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	// todo: check lengths of vectors
	
    // todo: check type of inputs
    if ( (mxIsSingle(prhs[0]) == 0) || (mxIsSingle(prhs[1]) == 0) || (mxIsSingle(prhs[2]) == 0) ) {
        tridiag_inv_mex_help();
        printf("inputs need to be single precision \n ");
        //Call(mxu_arg, (nrhs, prhs))
        Fail(Usage)
    }
    
    
    
    // todo: check complexity of inputs
    
    
    sub = (float *) mxGetData(prhs[0]);
    diag = (float *)  mxGetData(prhs[1]);
    sup =  (float *) mxGetData(prhs[2]);
    rhs = (float *) mxGetData(prhs[3]);
    N = mxGetM(prhs[1]); //remember col vec!
    M = mxGetN(prhs[3]);
    

    
    
    
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
    
/*     int out_size;
     out_size = (int) N *  (int) M;
     //plhs[0] = mxCreateNumericMatrix(out_size, 1, mxSINGLE_CLASS, mxREAL);
     //x_real = (float *)mxGetData(plhs[0]);
     for (int ii = 0; ii < out_size; ii++) {
         x_real[ii] = rhs[ii];
     } */
    
    if ((rhs_imag == NULL) ^ (x_imag == NULL)) {
        printf("problem: only one NULL vec, should be both or neither");
    }
    
	Call(tridiag_inv_mex_thr, (sub, diag, sup, rhs, rhs_imag, N, M, x_real, x_imag));
    //tridiag_inv_mex_thr(sub, diag, sup, rhs, rhs_imag, N, M, x_real, x_imag);
	
	Ok
}


//#if defined (Need_tridiag_inv_mex_gateway)
// gateway routine 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		tridiag_inv_mex_help();
		return;
	}
	if (!tridiag_inv_mex_gw(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("tridiag_inv_mex");
}

//#endif


