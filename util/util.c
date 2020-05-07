#include "util.h"
#include <math.h>
#include <stdio.h>

void pi_to_pi(double *angs, int nangs, int start) {
  // TODO optimize this with bit logic

  int i;
  for (i = start; i < nangs; i++) {
    if (angs[i] > 6.28 || angs[i] < -6.28)
      angs[i] = fmod(angs[i], 6.28);

    if (angs[i] > 3.14)
      angs[i] -= 6.28;
    else if (angs[i] < -3.14)
      angs[i] += 6.28;
  }
}

void multiply1(double x[3][3], double y[3][3], double (*res)[3][3]) {
  int i, j, k, temp;
  for (i = 0; i < 3; i++) {
    temp = 0;
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        temp = temp + x[i][k] * y[k][j];
      }
    }
    (*res)[i][j] = temp;
  }
}

void multiply2(double x[3][2], double y[2][2], double (*res)[3][2]) {
  int i, j, k, temp;
  for (i = 0; i < 3; i++) {
    temp = 0;
    for (j = 0; j < 2; j++) {
      for (k = 0; k < 2; k++) {
        temp = temp + x[i][k] * y[k][j];
      }
    }
    (*res)[i][j] = temp;
  }
}

void multiply3(double x[3][2], double y[2][3], double (*res)[3][3]) {
  int i, j, k, temp;
  for (i = 0; i < 3; i++) {
    temp = 0;
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 2; k++) {
        temp = temp + x[i][k] * y[k][j];
      }
    }
    (*res)[i][j] = temp;
  }
}

double det(double M[2][2]) { return M[0][0] * M[1][1] - M[0][1] * M[1][0]; }

void inv(double M[2][2], double (*r)[2][2]) {
  double d = det(M);
  (*r)[0][0] = M[1][1] / d;
  (*r)[0][1] = -M[1][0] / d;
  (*r)[1][0] = -M[0][1] / d;
  (*r)[1][1] = M[0][0] / d;
}

void multiply2by1(double M[2][2], double N[2], double (*r)[2]) {
  (*r)[0] = M[0][0] * N[0] + M[0][1] * N[0];
  (*r)[1] = M[1][0] * N[1] + M[1][1] * N[1];
}

double dot(double a[2], double b[2]) { return a[0] * b[0] + a[1] * b[1]; }
