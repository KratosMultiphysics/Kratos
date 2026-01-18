/**
 * Test matrix algebra in 3D.
 *
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/** Include the mmg3d library header file */
// if the header file is in the "include" directory
// #include "libmmg3d.h"
// if the header file is in "include/mmg/mmg3d"
#include "mmg/mmg3d/libmmg3d.h"
#include "mmgcommon_private.h"

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;

  fprintf(stdout,"  -- TEST MATRIX ALGEBRA IN 3D \n");

  if ( argc != 1 ) {
    printf(" Usage: %s\n",argv[0]);
    return(1);
  }


  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh:
   * MMG5_ARG_start: we start to give the args of a variadic func
   * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
   * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
   * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
   * &mmgSol: pointer toward your MMG5_pSol (that store your metric) */

  mmgMesh = NULL;
  mmgSol  = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /* Set inoffensive hmin and hmax for metric intersection */
  if( !MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hmin,1.e-6) )
    return EXIT_FAILURE;
  if( !MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hmax,1.e+6) )
    return EXIT_FAILURE;


  /** ------------------------------ STEP  II -------------------------- */

  /* matrix inversion test */
  if( !MMG5_test_invmat33() )
    return(EXIT_FAILURE);

  /* symmetric matrix eigendecomposition test */
  double m_sym[4][6] = {{2.,0.,0.,3.,4.,9.},
                        {1.,0.,0.,50.5,49.5,50.5},
                        {0.,1.,1.,0.,1.,0.},
                        {1.00495,-1.33511e-23,7.43048e-13,1.00495,-5.35377e-13,1.00502}}; /* Test matrices */
  double lambda_sym[4][3] = {{1.,2.,11.},
                             {1.,1.,100.},
                             {-1.,-1.,2.},
                             {1.00495,1.00495,1.00502}}; /* Exact eigenvalues */
  double vp_sym[4][3][3] = {{{0.,-2./sqrt(5.),1./sqrt(5.)},
                             {1.,0.,0.},
                             {0.,1./sqrt(5.),2./sqrt(5.)}},
                            {{1.,0.,0.},
                             {0.,1./sqrt(2.),-1./sqrt(2.)},
                             {0.,1./sqrt(2.),1./sqrt(2.)}},
                            {{0.,-1./sqrt(2.), 1./sqrt(2.)},
                             {2./sqrt(6.),-1./sqrt(6.),-1./sqrt(6.)},
                             {1./sqrt(3.), 1./sqrt(3.), 1./sqrt(3.)}},
                            {{1.,0.,0.},
                             {0.,1.,0.},
                             {0.,0.,1.}}}; /* Exact eigenvectors */
  for( int8_t i = 0; i < 3; i++ )
    if( !MMG5_test_eigenvmatsym3d(mmgMesh,m_sym[i],lambda_sym[i],vp_sym[i]) )
      return(EXIT_FAILURE);

  /* non-symmetric matrix eigendecomposition test */
  double m_nonsym[3][9] = {{500.5,-499.5,499.5,
                            -49.5,  50.5, 49.5,
                            450., -450., 550.},
                           {50.5,-49.5,49.5,
                             0.,   1.,  0.,
                            49.5,-49.5,50.5},
                           {1.00495,-1.33511e-23,7.43048e-13,
                            -1.33511e-23,1.00495,-5.35377e-13,
                            7.43048e-13,-5.35377e-13,1.00502}}; /* Test matrices */
  double lambda_nonsym[3][3] = {{1.,100.,1000.},
                                {1.,  1., 100.},
                                {1.00495,1.00495,1.00502}}; /* Exact eigenvalues */
  double vp_nonsym[3][3][3] = {{{1./sqrt(2.),1./sqrt(2.),0.},
                                {0.,         1./sqrt(2.),1./sqrt(2.)},
                                {1./sqrt(2.),         0.,1./sqrt(2.)}},
                               {{1./sqrt(2.),1./sqrt(2.),0.},
                                {0.,         1./sqrt(2.),1./sqrt(2.)},
                                {1./sqrt(2.),         0.,1./sqrt(2.)}},
                               {{1.,0.,0.},
                                {0.,1.,0.},
                                {0.,0.,1.}}}; /* Exact right eigenvectors */
  double ivp_nonsym[3][3][3] = {{{ 1./sqrt(2.),-1./sqrt(2.), 1./sqrt(2.)},
                                 { 1./sqrt(2.), 1./sqrt(2.),-1./sqrt(2.)},
                                 {-1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.)}},
                                {{ 1./sqrt(2.),-1./sqrt(2.), 1./sqrt(2.)},
                                 { 1./sqrt(2.), 1./sqrt(2.),-1./sqrt(2.)},
                                 {-1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.)}},
                                {{1.,0.,0.},
                                 {0.,1.,0.},
                                 {0.,0.,1.}}}; /* Exact right eigenvectors inverse */
  for( int8_t i = 0; i < 2; i++ )
    if( !MMG5_test_eigenvmatnonsym3d(mmgMesh,m_nonsym[i],lambda_nonsym[i],vp_nonsym[i],ivp_nonsym[i]) )
      return(EXIT_FAILURE);

  /* test 3x3 matrix transposition */
  if( !MMG5_test_transpose3d() )
    return(EXIT_FAILURE);

  /* test vector scalar product */
  if( !MMG5_test_dotprod() )
    return(EXIT_FAILURE);

  /* test vector product */
  if( !MMG5_test_crossprod3d() )
    return(EXIT_FAILURE);

  /* symmetric matrix multiplication test */
  if( !MMG5_test_mn() )
    return(EXIT_FAILURE);

  /* matrix linear transformation test */
  if( !MMG5_test_rmtr() )
    return(EXIT_FAILURE);

  /* rotation matrix test */
  if( !MMG5_test_rotmatrix() )
    return(EXIT_FAILURE);

  /* simultaneous reduction test */
  double m[2][6] = {{111./2.,-109./2.,  89./2.,111./2.,-91./2.,111./2.},
                    {370.15769255207715, 0., 0., 370.15769255207715, 0., 381.10501238659111}}; /* Test matrix 1 */
  double n[2][6] = {{409./2.,-393./2.,-407./2.,409./2.,391./2.,409./2.},
                    {371.98959740363608, -4.9420191999593769E-21, 2.7504505534267797E-10, 371.98959740363608, -1.9817373851183142E-10, 383.0188865492878 }}; /* Test matrix 2 */
  double dm[2][3] = {{1., 10.,100.},
                     {3.701532847548294e+02, 3.701621003493249e+02, 3.811050123865912e+02}}; /* Exact cobasis projection 1 */
  double dn[2][3] = {{8.,400.,  1.},
                     {3.719851677922659e+02, 3.719940270150061e+02, 3.830188865492879e+02}}; /* Exact cobasis projection 2 */
  double vp[2][3][3] = {{{1./sqrt(2.),1./sqrt(2.),0.},
                         {0.,         1./sqrt(2.),1./sqrt(2.)},
                         {1./sqrt(2.),         0.,1./sqrt(2.)}},
                        {{-7.070983648930086e-01, -7.071067772799207e-01,  1.959912696770420e-09},
                         {-7.071151974800862e-01,  7.071067849929998e-01,  1.204549101753515e-08},
                         {-1.018954010056916e-08,  7.341563327822146e-09, -1.000000000000000e+00}}}; /* Exact cobasis vectors */
  for( int8_t i = 0; i < 1; i++ )
    if( !MMG5_test_simred3d(mmgMesh,m[i],n[i],dm[i],dn[i],vp[i]) )
      return(EXIT_FAILURE);

  /* matrix inverse transformation test */
  if( !MMG5_test_updatemet3d_ani() )
    return EXIT_FAILURE;

  /* metrics intersection test */
  if( !MMG5_test_intersecmet33(mmgMesh) )
    return EXIT_FAILURE;

  /** ------------------------------ STEP III -------------------------- */

  /** 3) Free the MMG3D structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  return(0);
}
