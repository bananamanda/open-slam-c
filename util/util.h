/**
 * This file contains various helper functions, mostly matrix multiplication for small matricies of pre-known sizes.
 **/

void pi_to_pi(double *angle, int nang, int start);
void multiply1(double x[3][3], double y[3][3], double (*res)[3][3]);
void multiply2(double x[3][2], double y[2][2], double (*res)[3][2]);
void multiply3(double x[3][2], double y[2][3], double (*res)[3][3]);
double det(double M[2][2]);
void inv(double M[2][2], double (*r)[2][2]);
void multiply2by1(double M[2][2], double N[2], double (*r)[2]);
double dot(double *a, double *b);
