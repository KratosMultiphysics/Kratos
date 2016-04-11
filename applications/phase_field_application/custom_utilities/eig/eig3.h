
/* Eigen-decomposition for symmetric 3x3 real matrices.
   Public domain, copied from the public domain Java library JAMA. */

#ifndef _eig_h
#define _eig_h

#include "containers/array_1d.h"

namespace Kratos
{
    
/* Symmetric matrix A => eigenvectors in columns of V, corresponding
   eigenvalues in d. */
void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);

/*
* Compute principal stresses and direction using eig3
* http://barnesc.blogspot.de/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
* Remarks: sigma_1, sigma_2, sigma_3 is sorted
*/
void spectral_decomposition(double sigma_xx,
                            double sigma_yy,
                            double sigma_zz,
                            double sigma_xy,
                            double sigma_yz,
                            double sigma_zx,
                            double& sigma_1, double& sigma_2, double& sigma_3,
                            array_1d<double, 3>& dir1, array_1d<double, 3>& dir2, array_1d<double, 3>& dir3);

void spectral_decomposition(double sigma_xx,
                            double sigma_yy,
                            double sigma_zz,
                            double sigma_xy,
                            double sigma_yz,
                            double sigma_zx,
                            double& sigma_1, double& sigma_2, double& sigma_3);
    
}

#endif
