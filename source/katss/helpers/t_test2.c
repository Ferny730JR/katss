#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "toms708.h"

struct t_test2_aggregate {
	/* What we want to compute */
	double t_stat;
	double df;
	double pval;

	/* What we need to compute rolling t-test */
	double x_mean;
	double x_M2;
	unsigned int x_count;

	double y_mean;
	double y_M2;
	unsigned int y_count;
};
typedef struct t_test2_aggregate t_test2_aggregate;

t_test2_aggregate *
t_test2_create(void)
{
	t_test2_aggregate *aggregate = malloc(sizeof *aggregate);
	aggregate->t_stat  = 0.0;
	aggregate->df      = 0.0;
	aggregate->pval    = 0.0;
	aggregate->x_mean  = 0.0;
	aggregate->x_M2    = 0.0;
	aggregate->x_count = 0U;
	aggregate->y_mean  = 0.0;
	aggregate->y_M2    = 0.0;
	aggregate->y_count = 0U;

	return aggregate;
}

void
t_test2_destroy(t_test2_aggregate *aggregate)
{
	free(aggregate);
}

void
t_test2_update(t_test2_aggregate *aggregate, double x_value, double y_value)
{
	if(!isnan(x_value)) {
		aggregate->x_count++;

		/* Update the running mean */
		double delta = x_value - aggregate->x_mean;
		aggregate->x_mean += delta / aggregate->x_count;

		/* Update the running variance */
		double delta2 = x_value - aggregate->x_mean;
		aggregate->x_M2 += delta * delta2;
	}
	if(!isnan(y_value)) {
		aggregate->y_count++;

		/* Update the running mean */
		double delta = y_value - aggregate->y_mean;
		aggregate->y_mean += delta / aggregate->y_count;

		/* Update the running variance */
		double delta2 = y_value - aggregate->y_mean;
		aggregate->y_M2 += delta * delta2;
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
t_test2_finalize(t_test2_aggregate *aggregate)
{
	if(aggregate->x_count < 2 || aggregate->y_count < 2)
		return;

	/* Compute the sample (unbiased) variance */
	double x_var = aggregate->x_M2 / (aggregate->x_count - 1);
	double y_var = aggregate->y_M2 / (aggregate->y_count - 1);

	/* Compute the student's T-test */
	aggregate->t_stat = (aggregate->x_mean - aggregate->y_mean) / \
	  sqrt((x_var / aggregate->x_count) + (y_var / aggregate->y_count));

	/* Begin computing the number of degrees of freedom */
	double x_var_avg = x_var / aggregate->x_count;
	double y_var_avg = y_var / aggregate->y_count;

	double num = (x_var_avg + y_var_avg) * (x_var_avg + y_var_avg);
	double denom = (x_var_avg * x_var_avg) / (aggregate->x_count - 1) + \
	               (y_var_avg * y_var_avg) / (aggregate->y_count - 1);
	aggregate->df = num / denom;

	/* Compute the p-value */
	aggregate->pval = 2 * t_test_cdf(-fabs(aggregate->t_stat), aggregate->df, true, false);
}
