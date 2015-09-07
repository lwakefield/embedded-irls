#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/*#include "platform.h"*/

typedef struct {
    int m, n;
    double ** v;
} mat_t, mat;

mat matrix_new(int m, int n)
{
    mat x;
    x.v = malloc(sizeof(double) * m);
    x.v[0] = calloc(sizeof(double), m * n);
    int i;
    for (i = 0; i < m; i++)
     x.v[i] = x.v[0] + n * i;
    x.m = m;
    x.n = n;
    return x;
}

void matrix_delete(mat m)
{
    free(m.v[0]);
    free(m.v);
//    free(m);
}

void matrix_transpose(mat m)
{
    int i,j;
    for (i = 0; i < m.m; i++) {
        for (j = 0; j < i; j++) {
            double t = m.v[i][j];
            m.v[i][j] = m.v[j][i];
            m.v[j][i] = t;
        }
    }
}

mat matrix_copy(int m, int n, double a[][n])
{
    mat x = matrix_new(m, n);
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            x.v[i][j] = a[i][j];
    return x;
}

mat matrix_mul(mat x, mat y)
{
    /*if (x.n != y.m) return 0;*/
    mat r = matrix_new(x.m, y.n);
    int i, j, k;
    for (i = 0; i < x.m; i++)
        for (j = 0; j < y.n; j++)
            for (k = 0; k < x.n; k++)
                r.v[i][j] += x.v[i][k] * y.v[k][j];
    return r;
}

mat matrix_minor(mat x, int d)
{
    mat m = matrix_new(x.m, x.n);
    int i, j;
    for (i = 0; i < d; i++)
        m.v[i][i] = 1;
    for (i = d; i < x.m; i++)
        for (j = d; j < x.n; j++)
            m.v[i][j] = x.v[i][j];
    return m;
}

/* c = a + b * s */
double *vmadd(double a[], double b[], double s, double c[], int n)
{
    int i;
    for (i = 0; i < n; i++)
        c[i] = a[i] + s * b[i];
    return c;
}

/* m = I - v v^T */
mat vmul(double v[], int n)
{
    mat x = matrix_new(n, n);
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            x.v[i][j] = -2 *  v[i] * v[j];
    for (i = 0; i < n; i++)
        x.v[i][i] += 1;

    return x;
}

/* ||x|| */
double vnorm(double x[], int n)
{
    double sum = 0;
    int i;
    for (i = 0; i < n; i++) sum += x[i] * x[i];
    return sqrt(sum);
}

/* y = x / d */
double* vdiv(double x[], double d, double y[], int n)
{
    int i;
    for (i = 0; i < n; i++) y[i] = x[i] / d;
    return y;
}

double* mcol_from(mat m, double *v, int c, int s)
{
    int i;
    for (i = s; i < m.m; i++)
        v[i-s] = m.v[i][c];
    return v;
}

/* take c-th column of m, put in v */
double* mcol(mat m, double *v, int c)
{
    return mcol_from(m, v, c, 0);
}

void matrix_show(mat m)
{
    int i, j;
    for(i = 0; i < m.m; i++) {
        for (j = 0; j < m.n; j++) {
            printf(" %8.3f", m.v[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void vec_show(double *vec, int n)
{
    int i;
    for (i = 0; i < n; i++) {
        printf(" %8.3f\n", vec[i]);
    }
    printf("\n");
}

void householder_vector(double *beta, double *v, int n)
{
	double sigma = 0;
	int i;
	double x_0 = v[0];
	v[0] = 1;
	for (i = 1; i < n; i++) {
		sigma += v[i] * v[i];
	}
	if (sigma == 0 && x_0 >= 0) {
		*beta = 0;
	} else if (sigma == 0 && x_0 < 0) {
		*beta = -2;
	} else {
		double mu = sqrt(x_0 * x_0 + sigma);
		if (x_0 <= 0) {
			v[0] = x_0 - mu;
		} else {
			v[0] = -sigma / (x_0 + mu);
		}
        *beta = (2 * (v[0] * v[0])) / (sigma + v[0] * v[0]);
		for (i = n-1; i > -1; i--) {
			v[i] = v[i] / v[0];
		}
	}
}

void back_accumulate(mat m)
{

}

void householder(mat m, mat *R, mat *Q)
{
    mat z1;
    int i, j, k;
    for (k = 0; k < m.n && k < m.m - 1; k++) {
    	printf("%d loop\r\n", k);
        z1 = matrix_minor(m, k);

        double beta;
        int n = z1.m;
        double v[n-k];
        mcol_from(z1, v, k, k);
        householder_vector(&beta, v, n);

        mat beta_v_vt = matrix_new(n,n);
        for (i = k; i < n; i++) {
        	for (j = k; j < n; j++) {
        		beta_v_vt.v[i][j] = beta * v[i - k] * v[j - k];
        	}
        }
        mat update = matrix_mul(beta_v_vt, z1);
        for (i = k; i < m.m; i++) {
			for (j = k; j < m.n; j++) {
				m.v[i][j] = z1.v[i][j] - update.v[i][j];
			}
        }

        /*Save the householder vectors*/
        for (i = 1; i < n - k; i++) {
        	m.v[i + k][k] = v[i];
        }
        matrix_delete(update);
        matrix_delete(beta_v_vt);
        matrix_delete(z1);
    }
//    matrix_delete(z);
}


//#include "test_data.h"
double in[][3] = {
  { 12, -51,   4},
  {  6, 167, -68},
  { -4,  24, -41},
  { -1, 1, 0},
  { 2, 0, 3},
};

int main()
{
    /*init_platform();*/
//    int n = 402;
//    int p = 120;

	int n = 5;
	int p = 3;
    mat R, Q;
    mat x = matrix_copy(n, p, in);
    printf("Got a new matrix\r\n");
    householder(x, &R, &Q);
    matrix_show(x);

//    puts("Q"); matrix_show(Q);
//    puts("R"); matrix_show(R);

    // to show their product is the input matrix
//    mat m = matrix_mul(Q, R);
//    puts("Q * R"); matrix_show(m);
    puts("Done.\r\n");

    matrix_delete(x);
    /*matrix_delete(R);*/
    /*matrix_delete(Q);*/
//    matrix_delete(m);
    return 0;
}
