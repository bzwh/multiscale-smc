#ifndef PARAMS_H_INCLUDED
#define PARAMS_H_INCLUDED

#include <vector>
#include <gsl/gsl_rng.h>
#include <Eigen/Dense>
#include "Settings.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
/** \brief Container for all model parameters to be estimated.
 * Each replicate takes a newly generated parameter set.
 * Methods of generating new set are here, but choice belongs to 'Chain' algorithm
 * So any new design of parameter set should follow same naming of existing methods
 */
class Params  {
public:
  // METHODS
  Params();
  void storage_setup();                     //!< After model is spec'd, gets storage space ready
  void pscreendump();               //!< Probably dumps parameters to screen
  void pwrite(std::ofstream&);      //!< Write self to given file
  void setacc();                    //!< Set to Mike's estimates
  void setrand(gsl_rng*);           //!< Set to random values in min/max range
  int prior_check();                //!< With flat priors - just checking within min-max range
  void calc_cov(MatrixXd&,VectorXd&,int);
  void parse(const VectorXd&);
  void pread(std::ifstream&);
  void load_transmission();         //!< Load transmission outputs for within-farm model
  void wfm_samp(gsl_rng*);

  double pri;
  // DATA
  int pnum; //!< Total number of parameters being fitted
  int nreg; //!< Num regional parameter sets (just fitting 1 set!)
  int nspc; //!< Num species [usually 2 - cattle + sheep(UK)/pigs(JP)]
  int plaw; //!< toggle power law or linear scaling
  int pker; //!< Num kernel parameters
  int pper; //!< Num parameters per region
  int pdet; //!< Num parameters for delay to detection (gamma distributed)
  int pdcs; //!< Num parameters for DC cull - F (distributed?)
  int pdcf; //!< toggle for DC cull - f (const?)
  int pwfm; //!< toggle for fitting wfm gamma shape/mean latent+infectious periods, beta
  int pvac; //!< Num parameters for vaccination (sus and trans scaling)
  int rwfm; //!< Running within-farm model?
  // BETWEEN FARM PARAMETERS
  // Susceptibility & transmissibility parameters for each region
  // 0:Cubria, 1:Devon, 2:Wales, 3:Scotland, 4:Rest of England
  std::vector< std::vector<double> > sb;    //!< Susceptibility coefficient
  std::vector< std::vector<double> > tb;    //!< Transmissibility coefficient
  // Powerlaws
  std::vector< std::vector<double> > sp;    //!< Susceptibility exponent
  std::vector< std::vector<double> > tp;    //!< Transmissibility exponent
  // Vaccd pars
  std::vector<double> vc;                   //!< Scaling for sus and trans for vac'd farms
  // Kernel
  std::vector<double> ker;                  //!< Kernel parameters k0,d0,alpha

  // CONTROL PARAMETERS
  // Delays - detection/reporting/cull
  int delaydet;                             //!< Delay from infection to detection
  std::vector<double> ddet;                 //!< Gamma parameterised delay
  int delayipc;                             //!< Delay from detection to cull
  int delaydcp;                             //!< Delay from report to cull
  // DC culls
  double f;                                 //!< DC tracing accuracy
  std::vector<double> F;                    //!< DC intensity
  // Vaccinations??
  std::vector<double> vaccrange;            //!< inner/outer radii of vaccination annulus

  // WITHIN-FARM SEmInR parameters - But these were to be read from file and selected randomly...
  /*std::vector<int> m_cows;                   //!< Number of exposed classes
  std::vector<int> m_shps;                   //!< Number of exposed classes
  std::vector<double> sigma_cows;            //!< Mean latent period
  std::vector<double> sigma_shps;            //!< Mean latent period
  std::vector<int> n_cows;                   //!< Number of infectious classes
  std::vector<int> n_shps;                   //!< Number of infectious classes
  std::vector<double> gamma_cows;            //!< Mean infectious period
  std::vector<double> gamma_shps;            //!< Mean infectious period
  std::vector<double> beta_cows;             //!< Transmission parameter
  std::vector<double> beta_shps;             //!< Transmission parameter  */
  // But not anymore... now read the gamma-fitted marginal posteriors and resample
  std::vector< std::vector<double> > latk;    //!< Latent period shape gam point priors
  std::vector< std::vector<double> > latm;    //!< Latent period mean gam point priors
  std::vector< std::vector<double> > infk;    //!< Infectious period shape gam point priors
  std::vector< std::vector<double> > infm;    //!< Infectious period mean gam point priors
  std::vector< std::vector<double> > beta;    //!< Beta gam point priors

  std::vector<double> kE;     //!< Current sample's latent period shape
  std::vector<double> mE;     //!< Current sample's latent period mean
  std::vector<double> kI;     //!< Current sample's infectious period shape
  std::vector<double> mI;     //!< Current sample's infectious period mean
  std::vector< std::vector<double> > bt;
  Eigen::VectorXd par_vec;                    //!< Next proposed set - parsed in to above details





private:
  // Regional species-specific susceptibility range (flat prior)
  std::vector< std::vector<double> > sus_min; //!< Flat prior lower bound
  std::vector< std::vector<double> > sus_max; //!< Flat prior upper bound
  // Regional species-specific transmissibility range (flat prior)
  std::vector< std::vector<double> > trn_min; //!< Flat prior lower bound
  std::vector< std::vector<double> > trn_max; //!< Flat prior upper bound
  // Regional species-specific powerlaws
  double sp_min;  //!< Flat prior exponent bounds
  double sp_max;  //!< Flat prior exponent bounds
  double tp_min;  //!< Flat prior exponent bounds
  double tp_max;  //!< Flat prior exponent bounds
  double dcfmin;
  double dcfmax;
  std::vector<double> detmin;
  std::vector<double> detmax;
  std::vector<double> dcFmin;
  std::vector<double> dcFmax;
  std::vector<double> ker_min;
  std::vector<double> ker_max;
  std::vector<double> vcc_min;
  std::vector<double> vcc_max;
};

#endif
