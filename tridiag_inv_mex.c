/* mex file for tridiagonal solver 2013-11-09 
 assuming block diagonal with identical, real tridiagonal blocks
 so rhs is formulated as NxM matrix */


#include "mex.h"
#include "matrix.h"
//#include "def,type.h"
#include "defs-env.h"
#include "jf,mex,def.h" 
#include "jf,thread1.h"
#include "pthread.h"

//#if !defined(Need_tridiag_inv_mex_gateway)
//#include "tridiag_inv,def.h"
//#endif

#define Usage "usage error. see above"
#define NUM_THREADS 4 // number of cores // 4 for iv1, 2 for vega

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
	double *subdiag_ptr;
	double *diagvals_ptr;
	double *supdiag_ptr;
	double *rhs_ptr;
	double *out_ptr;
};

struct thread_data thread_data_array[NUM_THREADS];

//static sof tridiag_inv_thr(void *threadarg)
void *tridiag_inv_thr(void *threadarg)
{
	//tridiag_inv_worker_args args;
	int taskid;
	int N;
	double *a;
	double *b;
	double *c;
       	double *d;
       	double *x;
	double *new_c;
	double *new_d;
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
	taskid = my_data -> thread_id;
	N = my_data -> block_size;
	a = my_data -> subdiag_ptr;
	b = my_data -> diagvals_ptr;
	c = my_data -> supdiag_ptr;
	d = my_data -> rhs_ptr;
	x = my_data -> out_ptr;
/*
	// print out args to check
	printf("in thread... \n");
	printf("subdiag: \n");
	for (int ii = 0; ii < N-1; ii++) {
		printf("%f ", a[ii]);
	} 
	printf("\n");
	printf("diag: \n");
	for (int ii = 0; ii < N; ii++) {
		printf("%f ", b[ii]);
	} 
	printf("\n");
	printf("supdiag: \n");
	for (int ii = 0; ii < N-1; ii++) {
		printf("%f ", c[ii]);
	} 
	printf("\n");
	printf("rhs: \n");
	for (int ii = 0; ii < N; ii++) {
		printf("%f ", d[ii]);
	} */
	// need to malloc for x_real, x_imag, and then free them ??
	new_c = (double *) calloc (N - 1, sizeof(double));
	new_d = (double *) calloc (N, sizeof(double));

	*new_c = *c / *b; //new_c[0] = c[0]/b[0];
	*new_d = *d / *b; //new_d[0] = d[0]/b[0];

	int ii;
	double a_prev, new_c_prev;
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

	printf("task id: %d \n", taskid);
    /*
	printf("result: \n");
	for (ii = 0; ii < N; ii++) {
		printf("%f ", x[ii]);
	} 
	printf("\n"); */

//	pthread_mutex_lock(&mutexout);
	
//	pthread_mutex_unlock(&mutexout);
	free(new_c);
	free(new_d);

	pthread_exit(NULL);
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
double *subdiag_ptr, double *diagvals_ptr, double *supdiag_ptr, double *rhs_real_ptr, mwSize block_size, mwSize nblocks, double *out_real_ptr)
{
	//tridiag_inv_worker_args args;

	int ii;
	int big_N; // total number of entries, N*nblocks
	int rc;
	long t;
	pthread_attr_t attr;
	
	printf("block size : %d \n", block_size);
    printf("nblocks : %d \n", nblocks);
    printf("wubba lubba dub dub \n");
	pthread_t threads[NUM_THREADS];
	
    /*
	// print out args to check
	printf("subdiag: \n");
	for (ii = 0; ii < block_size - 1; ii++) {
		//printf("%d ", subdiag_ptr[ii]);
		printf("%f", *(subdiag_ptr + ii));
	} 
	printf("\n");
	printf("diag: \n");
	for (ii = 0; ii < block_size; ii++) {
		printf("%f ", diagvals_ptr[ii]);
	} 
	printf("\n");
	for (ii = 0; ii < block_size-1; ii++) {
		printf("%f ", supdiag_ptr[ii]);
	} 
	printf("\n");
	printf("rhs: \n");
	for (ii = 0; ii < block_size; ii++) {
		printf("%f ", rhs_real_ptr[ii]);
	} 
	printf("\n"); */

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
	//for (int th_id = 0; th_id < NUM_THREADS; th_id++)
	for (int th_id = 0; th_id < nblocks; th_id++) // for now bc nblocks < NUM_THREADS
	{
		//Call(tridiag_inv_thr, ());
		thread_data_array[th_id].thread_id = th_id;
		thread_data_array[th_id].block_size = block_size;
        thread_data_array[th_id].subdiag_ptr = subdiag_ptr;// + th_id * block_size;
        thread_data_array[th_id].diagvals_ptr = diagvals_ptr;// + th_id * block_size;
        thread_data_array[th_id].supdiag_ptr = supdiag_ptr;// + th_id * block_size;
		thread_data_array[th_id].rhs_ptr = rhs_real_ptr + th_id * block_size;
		thread_data_array[th_id].out_ptr = out_real_ptr + th_id * block_size;
		rc = pthread_create(&threads[th_id], &attr, tridiag_inv_thr, (void *) &(thread_data_array[th_id]));
        printf("done with thread %d \n", th_id);
		//rc = pthread_create(&threads[th_id], &attr, tridiag_inv_thr, (void *) th_id);
	}

	/*
	//for (int th_id = 0; th_id < NUM_THREADS; th_id++)
	for (int th_id = 0; th_id < nblocks; th_id++) // for now bc nblocks < NUM_THREADS
	{
		rc = pthread_join(threads[th_id], NULL); // can check rc to see if success      
		if (rc) {
         		printf("ERROR; return code from pthread_join() is %d\n", rc);
         		exit(-1);
        	}
	} */
	Ok
}


// intermediate GateWay routine 
static sof tridiag_inv_mex_gw(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
    double *sub;              /* input subdiagonal */
    double *diag;               /* 1xN input diagonal */
    double *sup;
    double *rhs;
    size_t N;                   /* size of tridiag matrix */
    size_t M;                   /* numcols of rhs matrix */
    double *x;                  /* output */
    
	if (nrhs != 4 ) { // hard coding :(
		tridiag_inv_mex_help();
		//Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	// todo: check lengths of vectors
	//
    
    
    sub = mxGetPr(prhs[0]);
    diag = mxGetPr(prhs[1]);
    sup = mxGetPr(prhs[2]);
    rhs = mxGetPr(prhs[3]);
    N = mxGetM(prhs[1]); //remember col vec!
    M = mxGetN(prhs[3]);
    
/*    if (mxIsComplex(prhs[3])) {
        printf("rhs is complex \n");
        double *rhs_imag;
        rhs_imag = mxGetPi(prhs[3]);
        //out_imag_ptr = mxGetPi(plhs[0]);
        plhs[0] = mxCreateDoubleMatrix((mwSize)N, (mwSize) M, mxCOMPLEX);
    } else { */
        plhs[0] = mxCreateDoubleMatrix((mwSize)N, (mwSize) M, mxREAL);
        
//    }
    x = mxGetPr(plhs[0]);
    
   // check_types_and_sizes(sub, diag, sup, rhs, N, M);
    

	
	//Call(tridiag_inv_mex_thr, (nlhs, plhs, nrhs, prhs));
    tridiag_inv_mex_thr(sub, diag, sup, rhs, N, M, x);
	
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


