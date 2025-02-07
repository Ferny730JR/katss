struct minunit {
	int unit_tests_failed;
	int unit_tests_total;
};

#define UTEST_TYPE \
	struct minunit

#define init_unit_tests(testname) \
	printf("%-50s\n", testname); \
	int unit_tests_failed = 0; \
	int unit_tests_total = 0;

#define unit_tests_end \
	return (struct minunit){.unit_tests_failed=unit_tests_failed, \
	                        .unit_tests_total=unit_tests_total};

#define mu_assert(message, test) \
	do { \
		printf("    %-46s", message); \
		unit_tests_total++; \
		if(test) { \
			printf("OK\n"); \
		} else { \
			unit_tests_failed++; \
			printf("FAIL\n"); \
		} \
	} while(0)

#define init_run_test \
	struct minunit test_run_return; \
	int tests_failed = 0; \
	int tests_total = 0;

#define run_test_end \
	printf("TOTAL PASSED %d/%d\n", tests_total-tests_failed, tests_total)

#define mu_run_test(test) \
	do { \
		test_run_return = test(); \
		tests_failed += test_run_return.unit_tests_failed; \
		tests_total  += test_run_return.unit_tests_total; \
		printf("    PASSED %d/%d\n", test_run_return.unit_tests_total - \
		  test_run_return.unit_tests_failed, \
		  test_run_return.unit_tests_total); \
	} while (0)
