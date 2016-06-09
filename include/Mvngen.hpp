#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace Eigen;

VectorXd mgen(const MatrixXd&, gsl_rng*);
VectorXd mgen_dcmp(const MatrixXd&, gsl_rng*);
MatrixXd dcmp(const MatrixXd&);
