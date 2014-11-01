typedef double *vector;

typedef struct
{
	double R;
	double *T;
} pi_R_T_type;

typedef struct
{
	int cnt;
	int *col;
	pi_R_T_type *val;
	int *lamb;
	int *ck; //record the iteration # of the sigma
} sigma_type;

typedef struct
{
	int *col;
	pi_R_T_type **val;
} delta_type;

typedef struct
{
	int delta;
	int sigma;
} i_type;

void CudaMain(void);
void launch_kernel(sigma_type *sigma, delta_type *delta, i_type *istar, double *argmax, vector pi_Tbar_x, vector Xvect, int rv_cols, int ictr, int count);