#ifndef KATSS_T_TEST_H
#define KATSS_T_TEST_H


/**
 * @brief Two sample T-test aggregate
 */
struct t_test2_aggregate {
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


/**
 * @brief One-sample T-test aggregate
 * 
 */
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

/**
 * @brief Create the student's T-test aggregate.
 * 
 * @return t_test2_aggregate* 
 */
t_test2_aggregate *t_test2_create(void);


/**
 * @brief Free all resources allocated to the aggregate.
 * 
 * @param aggregate T-test aggregate to destroy
 */
void t_test2_destroy(t_test2_aggregate *aggregate);


/**
 * @brief Add new values to update the T-test of an aggregate.
 * 
 * @param aggregate Aggregate to update by adding values
 * @param x_value X value to add, NAN if none
 * @param y_value Y value to add, NAN if none
 */
void t_test2_update(t_test2_aggregate *aggregate, double x_value, double y_value);


/**
 * @brief Compute the two-sample T-test.
 * 
 * Computes t-statistic, degrees of freedom, and p-value, and stores it into
 * the aggregate based on the values that have been added to the aggregate.
 * 
 * @param aggregate Which aggregate to compute from
 */
void t_test2_finalize(t_test2_aggregate *aggregate);


/**
 * @brief Create the student's T-test aggregate.
 * 
 * @return t_test2_aggregate* 
 */
t_test2_aggregate *t_test1_create(void);


/**
 * @brief Free all resources allocated to the aggregate.
 * 
 * @param aggregate T-test aggregate to destroy
 */
void t_test1_destroy(t_test1_aggregate *aggregate);


/**
 * @brief Add new values to update the T-test of an aggregate.
 * 
 * @param aggregate Aggregate to update by adding values
 * @param x_value value to add
 */
void t_test1_update(t_test1_aggregate *aggregate, double value);


/**
 * @brief Compute the two-sample T-test.
 * 
 * Computes t-statistic, degrees of freedom, and p-value, and stores it into
 * the aggregate based on the values that have been added to the aggregate.
 * 
 * @param aggregate Which aggregate to compute from
 * @param mu0 Hypothesized population mean
 */
void t_test1_finalize(t_test1_aggregate *aggregate, double mu0);

#endif // KATSS_T_TEST_H
