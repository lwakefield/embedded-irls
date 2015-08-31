#include <stdio.h>
#include <gsl/gsl_multifit.h>
/*#include <gsl/gsl_randist.h>*/
#include <gsl/gsl_errno.h>

#include "testdata.h"
#include "hamming320.h"

/*int*/
/*dofit(const gsl_multifit_robust_type *T,*/
      /*const gsl_matrix *X, const gsl_vector *y,*/
      /*gsl_vector *c, gsl_matrix *cov)*/
/*{*/
  /*int s;*/
  /*gsl_multifit_robust_workspace * work */
    /*= gsl_multifit_robust_alloc (T, X->size1, X->size2);*/
  /*work->maxiter = 100;*/

  /*s = gsl_multifit_robust (X, y, c, cov, work);*/
  /*gsl_multifit_robust_free (work);*/

  /*return s;*/
/*}*/

const size_t n = 320;
const size_t p = 1;

gsl_vector_complex* load_spectrum()
{
    gsl_vector_complex *spectrum;
    spectrum = gsl_vector_complex_alloc (n);

    FILE *f;
    f = fopen("spectrum_in", "r");
    gsl_vector_complex_fscanf(f, spectrum);
    fclose(f);

    return spectrum;
}

gsl_vector* load_valid_points()
{
    gsl_vector *points;
    points = gsl_vector_alloc (n);

    FILE *f;
    f = fopen("valid_points", "r");
    gsl_vector_fscanf(f, points);
    fclose(f);

    return points;
}

int main (int argc, char **argv)
{
    int i;

    gsl_vector_complex *spectrum;
    gsl_vector *valid_points;

    spectrum = load_spectrum();
    valid_points = load_valid_points();

    gsl_vector_complex_view hamming320
        = gsl_vector_complex_view_array (hamming320_data, n);

    gsl_vector_complex_mul(spectrum, &hamming320.vector);

    gsl_vector_complex_free(spectrum);
    gsl_vector_free(valid_points);
    printf("done.\n");
    return 0;
}
