#include <stdio.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>

#include "simpletestdata.h"

int
dofit(const gsl_multifit_robust_type *T,
      const gsl_matrix *X, const gsl_vector *y,
      gsl_vector *c, gsl_matrix *cov)
{
  int s;
  gsl_multifit_robust_workspace * work 
    = gsl_multifit_robust_alloc (T, X->size1, X->size2);
  work->maxiter = 100;

  s = gsl_multifit_robust (X, y, c, cov, work);
  gsl_multifit_robust_free (work);

  return s;
}

int main (int argc, char **argv)
{
  int i;
  const size_t n = 10;
  const size_t p = 1;

  gsl_matrix *cov;
  gsl_vector *c, *c_ols;

  gsl_matrix_view y 
    = gsl_matrix_view_array (y_data, n, p);

  gsl_vector_view x
    = gsl_vector_view_array (x_data, n);

  c = gsl_vector_alloc (p);
  c_ols = gsl_vector_alloc (p);
  cov = gsl_matrix_alloc (p, p);

  gsl_set_error_handler_off();
  dofit(gsl_multifit_robust_bisquare, &y.matrix, &x.vector, c, cov);

  FILE *f;
  f = fopen("reconstruct_soln", "w");
  gsl_vector_fprintf (f, c, "%g");
  gsl_vector_fprintf (stdout, c, "%g");
  fclose(f);

  printf("done.\n");
  return 0;
}
