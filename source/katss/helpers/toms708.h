#ifndef TOMS708_H
#define TOMS708_H


/**
 * @brief Compute the incomplete function beta ratio
 * 
 * @param a 
 * @param b 
 * @param x 
 * @param y 
 * @param w 
 * @param w1 
 * @param ierr 
 * @param log_p 
 */
void
bratio(double a, double b, double x, double y, double *w, double *w1,
       int *ierr, int log_p);

#endif // TOMS708_H