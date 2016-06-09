#include <iostream>
#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Mvngen.hpp"
using namespace Eigen;
using std::cout;
using std::endl;


/** \brief Given variance-covariance matrix, sample from zero mean multinormal distribution
 * \param sigma MatrixXd
 * \param r gsl_rng*
 * \return VectorXd
 */
VectorXd mgen(const MatrixXd& sigma, gsl_rng* r)  {
  int sz = sigma.rows();
  MatrixXd A = dcmp(sigma);

  // Get the standard normals to be transformed
  VectorXd randed(sz);
  for (unsigned int i=0;i<randed.size();++i)  {
    randed(i) = gsl_ran_ugaussian(r);
  }
  return(A*randed); // To be added to mu
}

/** \brief Decompose given matrix. Tries Choleksy, if fails does eigendecomp instead
 * \param sigma const MatrixXd&
 * \return MatrixXd
 */
MatrixXd dcmp(const MatrixXd& sigma)  {
  int sz = sigma.rows();  // .size() gives number of elements
  MatrixXd A(sz,sz);
  // Positive-definite can just use Cholesky decomposition
  LLT<MatrixXd> llt;
  llt.compute(sigma);
  if (llt.info()==Success)  {
    //cout << "Choleskied!" << endl;
    A = llt.matrixL();
  }
  // Positive semi-definite have to Eigendecomposition
  else  {
    //cout << "Not choleskied ):" << endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(sigma);
    A = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }
  return(A);
}

/** \brief Generate multinormal from factorised covariance matrix. Zero mean.
 * \param A eigen::matrix decomposition of the covariacne matrix ie cov=AA'
 * \param r pointer to prng
 * \return VectorXd
 */
VectorXd mgen_dcmp(const Eigen::MatrixXd& A, gsl_rng* r)  {
  int sz = A.rows();
  // Get the standard normals to be transformed
  VectorXd randed(sz);
  for (unsigned int i=0;i<randed.size();++i)  {
    randed(i) = gsl_ran_ugaussian(r);
  }
  return(A*randed);
}









