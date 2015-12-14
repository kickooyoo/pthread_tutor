/* mex file for tridiagonal solver 2013-11-09 */

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
#define NUM_THREADS 4 // number of cores

pthread_mutex_t mutexout; // global var for locking

//#if defined(Mmex)

static void tridiag_inv_mex_help(void)
{
	printf("Usage for tridiag_inv_mex: \n\
	output = tridiag_inv_mex(subdiag,diagvals,supdiag,argum) \n\
	subdiag has length n-1	");
}

/*typedef struct 
{
	float *yy;
	int Ny; // etc??
} tridiag_inv_worker_args; */


struct thread_data
{
	int thread_id;
	int N; 
	double *subdiag_ptr;
	double *diagvals_ptr;
	double *supdiag_ptr;
	double *rhs_ptr;
	double *out_ptr;
};

struct thread_data thread_data_array[NUM_THREADS];

static sof tridiag_inv_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])	
{
	#define x plhs[0]
	//Const mxArray *a; //subdiag
	//Const mxArray *b; //diagvals
	//Const mxArray *c; //supdiag
	double *a; //subdiag
	double *b; //diagvals
	double *c; //supdiag
	double *d; //argum
	double a_prev, new_c_prev;

	mxArray *new_c; 
	mxArray *new_d;
	//mxArray *x;
	double *new_c_db;
	double *new_d_db;
	double *x_db;

	int N; // size of block
	int ndims;
	int output_dims[ndims];
	
	// todo: check data type, real/complex
	a = mxGetPr(prhs[0]);
	b = mxGetPr(prhs[1]);
	c = mxGetPr(prhs[2]);
	d = mxGetPr(prhs[3]);

	//x = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
	x_db = mxGetPr(x);

	new_c = mxCreateDoubleMatrix(N, 1, mxREAL);// todo: set to zero
	new_c_db = mxGetPr(new_c);
	new_d = mxCreateDoubleMatrix(N, 1, mxREAL);// todo: set to zero
	new_d_db = mxGetPr(new_d);

	N = mxGetM((mxArray *) d);
	printf("M: %i, N: %i", mxGetM((mxArray *) d), mxGetN((mxArray *) d));

	// need to do the following for real and imag components of d

	*new_c_db = (*c) / (*b); //new_c[0] = c[0]/b[0];
	*new_d_db = (*d) / (*b); //new_d[0] = d[0]/b[0];

	int ii;
	for (ii = 1; ii <= N-2; ii++) {
		a_prev =  *(a + ii -1);
		new_c_prev =  *(new_c_db + ii -1);
		*(new_c_db + ii) = *(c + ii) / (*(b + ii) - new_c_prev * a_prev);
		//new_c[ii] = c[ii]/(b[ii]-new_c[ii-1]*a[ii-1]);	
		*(new_d_db + ii) = (*(d + ii) - *(new_d_db + ii -1) * a_prev) / (*(b + ii) - new_c_prev * a_prev);
		//new_d[ii] = (d[ii]-new_d[ii-1]*a[ii-1])/(b[ii]-new_c[ii-1]*a[ii-1]);
	}
	*(new_d_db + N - 1) = (*(d + N - 1) - *(new_d_db + N - 2) * (*(a + N -1))) / (*(b + N - 1) - *(new_c_db + N - 2) * (*(a + N - 2)));
	//new_d[N-1] = (d[N-1]-new_d[N-2]*a[N-1])/(b[N-1]-new_c[N-2]*a[N-2]);

	*(x_db + N - 1) = *(new_d_db + N -1);
	//x[N-1] = new_d[N-1];
	for (ii = N-2; ii >= 1; ii--) {
		*(x_db + ii) = *(new_d_db + ii) - *(new_c_db + ii) * (*(x_db + ii +1));
		//x[ii] = new_d[ii]-new_c[ii]*x[ii+1];
	}
	
	/*ndims = 2;
	output_dims[0] = N; 
	output_dims[1] = mxGetN(d); // should be 1
	*/ //FIGURE THIS SHIT OUT
	//Call(plhs[0] = mxCreateNumericArray, (ndims, output_dims, mxDOUBLE_CLASS,mxCOMPLEX))
	
	//Call(tridiag_inv, (args)) // why???
	Ok
}

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
	N = my_data -> N;
	a = my_data->subdiag_ptr;
	b = my_data -> diagvals_ptr;
	c = my_data -> supdiag_ptr;
	d = my_data -> rhs_ptr;
	x = my_data -> out_ptr;

	// need to malloc for x_real, x_imag, and then free them ??
	new_c = (double *) calloc (N - 1, sizeof(double));
	new_d = (double *) calloc (N, sizeof(double));

	*new_c = *c / *b; //new_c[0] = c[0]/b[0];
	*new_d = *d / *b; //new_d[0] = d[0]/b[0];

	int ii;
	int a_prev, new_c_prev;
	for (ii = 1; ii <= N-2; ii++) {
		a_prev = *(a + ii - 1);
		new_c_prev = *(new_c + ii - 1);
		*(new_c + ii) = *(c + ii) / (*(b + ii) - new_c_prev * a_prev);
		*(new_d + ii) = (*(d + ii) - *(new_d + ii - 1) * a_prev) / (*(b + ii) - new_c_prev * a_prev);
	}
	*(new_d + N - 1) = (*(d + N - 1) - *(new_d + N - 2) * (*(a + N - 1))) / (*(b + N - 1) - *(new_c + N - 2) * (*(a + N - 2)));

	*(x + N - 1) = *(new_d + N - 1);
	for (ii = N-2; ii >= 1; ii--) {
		*(x + ii) = *(new_d + ii) - *(new_c + ii) * (*(x + ii +1));
	}

//	pthread_mutex_lock(&mutexout);
	
//	pthread_mutex_unlock(&mutexout);

	pthread_exit(NULL);
}

static sof tridiag_inv_mex_thr(
int nlhs, 
mxArray *plhs[],
Const mxArray *subdiag,
Const mxArray *diagvals,
Const mxArray *supdiag,
Const mxArray *rhs,
int N)	
{
	//tridiag_inv_worker_args args;
	double *subdiag_ptr;
	double *diagvals_ptr;
	double *supdiag_ptr;
	double *rhs_real_ptr;
	double *rhs_imag_ptr;
	double *out_real_ptr;
	double *out_imag_ptr;
	int ii;
	int big_N; // total number of entries, N*nblocks
	int nblocks;
	int rc;

	long t;
	pthread_attr_t attr;

	subdiag_ptr = mxGetPr(subdiag);
	diagvals_ptr = mxGetPr(diagvals);
	supdiag_ptr = mxGetPr(supdiag);
	rhs_real_ptr = mxGetPr(rhs);
	rhs_imag_ptr = mxGetPi(rhs);
	big_N = mxGetM(diagvals);
	nblocks = big_N / N;

	pthread_t threads[NUM_THREADS];
	

	plhs[0] = mxCreateNumericMatrix(big_N, 1, mxDOUBLE_CLASS, mxCOMPLEX);
	out_real_ptr = mxGetPr(plhs[0]);
	out_imag_ptr = mxGetPi(plhs[0]);

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	for (int th_id = 0; th_id < NUM_THREADS; th_id++)
	{
		//Call(tridiag_inv_thr, ());
		thread_data_array[th_id].thread_id = th_id;
		thread_data_array[th_id].N = N;
		thread_data_array[th_id].subdiag_ptr = subdiag_ptr + th_id * N;
		thread_data_array[th_id].diagvals_ptr = diagvals_ptr + th_id * N;
		thread_data_array[th_id].supdiag_ptr = supdiag_ptr + th_id * N;
		thread_data_array[th_id].rhs_ptr = rhs_real_ptr + th_id * N;
	//	thread_data_array[th_id].rhs_ptr = rhs_imag + th_id * N;
		thread_data_array[th_id].out_ptr = out_real_ptr + th_id * N;
	//	thread_data_array[th_id].out_imag_ptr = out_imag_ptr + th_id * N;
		//rc = pthread_create(&threads[th_id], &attr, tridiag_inv_thr, (void *) &(thread_data_array[th_id]));
		rc = pthread_create(&threads[th_id], &attr, tridiag_inv_thr, (void *) th_id);
	}

	
	for (int th_id = 0; th_id < NUM_THREADS; th_id++)
	{
		rc = pthread_join(threads[th_id], NULL); // can check rc to see if success
	}
	Ok
}


// intermediate GateWay routine 
static sof tridiag_inv_mex_gw(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	if (nrhs != 4 ) {
		tridiag_inv_mex_help();
		//Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	// todo: check lengths of vectors
	//
	
	//Call(tridiag_inv_mex, (nlhs, plhs, nrhs, prhs));
	//Call(tridiag_inv_mex_thr, (nlhs, plhs, nrhs, prhs));
	tridiag_inv_mex_thr(nlhs, plhs, prhs[0], prhs[1], prhs[2], prhs[3], mxGetScalar(prhs[4]));	
	
	Ok 
}


#if defined (Need_tridiag_inv_mex_gateway)
// gateway routine 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		tridiag_inv_mex_help();
		return;
	}
	if (!tridiag_inv_mex_gw(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("penalty_mex");
}

#endif


