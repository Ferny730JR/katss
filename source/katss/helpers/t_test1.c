#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "toms708.h"

struct t_test1_aggregate {
	/* What we want to compute */
	double t_stat;
	double df;
	double pval;

	/* What we need to compute rolling t-test */
	double mean;
	double M2;
	unsigned int count;
};
typedef struct t_test1_aggregate t_test1_aggregate;

t_test1_aggregate *
t_test1_create(void)
{
	t_test1_aggregate *aggregate = malloc(sizeof *aggregate);
	aggregate->t_stat  = 0.0;
	aggregate->df      = 0.0;
	aggregate->pval    = 0.0;
	aggregate->mean    = 0.0;
	aggregate->M2      = 0.0;
	aggregate->count   = 0U;

	return aggregate;
}

void
t_test1_destroy(t_test1_aggregate *aggregate)
{
	free(aggregate);
}

void
t_test1_update(t_test1_aggregate *aggregate, double value)
{
	if(!isnan(value)) {
		aggregate->count++;

		/* Update the running mean */
		double delta = value - aggregate->mean;
		aggregate->mean += delta / aggregate->count;

		/* Update the running variance */
		double delta2 = value - aggregate->mean;
		aggregate->M2 += delta * delta2;
	}
}

/**
 * @brief Regularized incomplete beta function
 * 
 * @param x 
 * @param a 
 * @param b 
 * @param lower_tail 
 * @param log_p 
 * @return double 
 */
static double
reg_incomplete_beta(double x, double a, double b, int lower_tail, int log_p)
{
	double x1 = 0.5 - x + 0.5, w, wc;
	int ierr;
	bratio(a, b, x, x1, &w, &wc, &ierr, log_p);
	return lower_tail ? w : wc;
}

/**
 * @brief T-distribution CDF with `df` degrees of freedom.
 * 
 * @param t   Value of T-test statistic
 * @param df  Number of degrees of freedom
 * @param lower_tail Probabilities from lower tail or upper tail
 * @param log_p Probabilities p are given as log(p)
 * @return double Probability
 */
static double
t_test_cdf(double t, double df, bool lower_tail, bool log_p)
{
		double nx = 1 + (t/df)*t;

	double val = (df > (t * t))
	    ? reg_incomplete_beta(t*t/(df+t*t), 0.5, df/2.0, 0, log_p)
	    : reg_incomplete_beta(1.0 / nx,     df/2.0, 0.5, 1, log_p);

	if(t<=0)
		lower_tail = !lower_tail;
	if(log_p) {
		if(lower_tail) {
			return log1p(-0.5*exp(val));
		} else {
			return val - M_LN2;
		}
	} else {
		val /= 2;
		// return 1-val if lower_tail, otherwise val
		return lower_tail ? 0.5 - val + 0.5 : val;
	}
}

void
t_test1_finalize(t_test1_aggregate *aggregate, double mu0)
{
	if(aggregate->count < 2)
		return;

	/* Compute the sample (unbiased) variance */
	double variance = aggregate->M2 / (aggregate->count - 1);

	/* Compute the student's T-test */
	aggregate->t_stat = (aggregate->mean - mu0) / sqrt(variance / aggregate->count);

	/* Begin computing the number of degrees of freedom */
	aggregate->df = aggregate->count - 1;

	/* Compute the p-value */
	aggregate->pval = 2 * t_test_cdf(-fabs(aggregate->t_stat), aggregate->df, true, false);
}
