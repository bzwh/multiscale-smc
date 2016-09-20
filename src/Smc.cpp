#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <omp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "Smc.hpp"
#include "Model.hpp"
#include "Farm.hpp"
#include "Mvngen.hpp"


#define LFLAG 1
#define CFLAG 0

using namespace Eigen;
using std::cout;
using std::endl;
using std::vector;
using std::flush;
using std::ofstream;
using std::to_string;
using std::pow;
using std::exp;
using std::stringstream;
using std::string;
using std::setfill;
using std::setw;

/** \brief Constructor.
 * \param nc int Number of threads.
 * \param nn int Number of particles.
 * \param ss int Number of rounds of filtering
 * \param frms Farms& Big chunk of const data, just pass through to Model (avoid reading repeatedly)
 */
Smc::Smc(int nc,int nn,int ss,Farms& frms)  {
  nthreads = nc;
  r.resize(nthreads);
  for (int cno=0;cno<nthreads;++cno)  {
    r[cno] = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r[cno],(cno+1.0)*time(0));
  }
  N = nn;
  nround = ss;
  models.reserve(nthreads);
  for (int cno=0;cno<nthreads;++cno)  {
    models.push_back(Model(frms));
    models[cno].setup(r[cno]);
  }
  npar = models[0].params.pnum; // == particles[0].size()
  nmet = frms.enreg;                 // Number of metrics to test against...
  particles_old = MatrixXd::Zero(N,npar);
  particles_new = MatrixXd::Zero(N,npar);
  priors.resize(N);
  weights_old.resize(N);
  weights_new.resize(N);
  means.resize(npar);
  sigma = MatrixXd::Zero(npar,npar);
  tolerance.resize(nmet);
  tolerance.fill(10.0);
  //tolerance << 0.105,0.165,0.180,0.140,0.08;
  epsilons = MatrixXd::Zero(N,nmet);
  step = 0;
}

/*
void Smc::restart()  {
  int old_round = 4;
  ifstream wfile("./outputs/wht.txt");
  ifstream sfile("./outputs/sig.txt");
  ifstream efile("./outputs/err"+to_string(old_round)+".txt");
}
*/

/** \brief Loop rounds of filtering. Opens output files and calls write();
 * \return void
 */
void Smc::run()  {
  ofstream wht("./outputs/wht.txt");
  ofstream sig("./outputs/sig.txt");
  for (step=0;step<nround;++step)  {
    stringstream d_tmp;
    d_tmp  << "./outputs/dat" << setfill('0') << setw(2) << to_string(step) << ".txt";
    ofstream dat(d_tmp.str());
    stringstream e_tmp;
    e_tmp  << "./outputs/err" << setfill('0') << setw(2) << to_string(step) << ".txt";
    ofstream eps(e_tmp.str());
    cout << "ROUND: "<<step << endl;
    iterate();
    update();
    cout << "\n" << endl;
    write(dat,eps,wht,sig);
    dat.close();
    eps.close();
  }
  wht.close();
  sig.close();
  rngfree();
}


/** \brief Dump state after each round
 * \param dat ofstream & For particles
 * \param eps ofstream & For metric values
 * \param wht ofstream & Particle weights
 * \param sig ofstream & Dump particle means and (component-wise) variances
 * \return void
 */
void Smc::write(ofstream& dat, ofstream& eps, ofstream& wht, ofstream& sig)  {
  dat << particles_old << endl;
  wht << ess << "\t" << weights_old.transpose() << "\n"<<flush;
  eps << tolerance.transpose() << "\n" << epsilons << "\n\n"<<flush;
  sig << means <<"\n\n" << sigma << "\n\n" << sginv << endl;
}

/** \brief Looping through all particles in this round. Spawns however many threads to run through.
 * \return void
 */
void Smc::iterate()  {
  #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
  for (int i=0;i<N;++i)  {
    stringstream lstream;
    int cno = omp_get_thread_num();
    lstream << i << "-" << cno << ": ";
    cout << lstream.str() << flush;
    int acc = 0;
    while (acc==0)  {
      gen(i,cno);
      runmodel(i,cno);
      acc = test(i,cno);
    }
    priors(i) = models[cno].params.pri;
    cout << models[cno].Itot << endl;
    lstream.str(string());
  }
  cout << "Tolerance: " << tolerance.transpose() << endl;
}


/** \brief Generate a new particle. Random to start OTW resample+perturb
 * \param i int Particle id
 * \param cno int Which thread/Model/RNG to use
 * \return void
 */
void Smc::gen(int i,int cno)  {
  models[cno].resetsim();
  int pcheck = 1;
  while (pcheck!=0)  {
    int ready = 0;
    if (step==0)  {
      particles_old.row(i) = models[cno].init_samp(ready);
      particles_new.row(i) = particles_old.row(i);
    }
    else  {
      vector<unsigned int> samp(N,0);
      gsl_ran_multinomial(r[cno],N,1,&weights_old[0],&samp[0]);
      int j = find(samp.begin(),samp.end(),1)-samp.begin();
      particles_new.row(i) = particles_old.row(j) + perturbation(cno).transpose();
    }
    /*stringstream lstream;
    lstream << particles_new.row(i) << "\n";
    cout << lstream.str() << flush;
    lstream.str(string());*/
    models[cno].parse(particles_new.row(i));
    pcheck = models[cno].params.prior_check();  // 0 for all good
    //cout << pcheck << flush;
  }
}


/** \brief Generate perturbation using component-wise variances.
 * \param cno int
 * \return VectorXd
 */
VectorXd Smc::perturbation(int cno)  {
  if (CFLAG)  {
    return(mgen(0.68*sigma,r[cno]));
  }
  else  {
    VectorXd v(npar);
    for (int p=0;p<npar;++p)  {
      v(p) = gsl_ran_gaussian(r[cno],0.68*(sginv(p)));
    }
    return(v);
  }
}


/** \brief Pretty self explanatory...
 * \param i int Particle id
 * \param cno int Thread id
 * \return void
 */
void Smc::runmodel(int i,int cno)  {
  models[cno].run(particles_new.row(i));
}


/** \brief Test simulation against all metrics
 * \param i int
 * \param cno int
 * \return int
 */
int Smc::test(int i,int cno)  {
  VectorXd v = models[cno].errcalc();
  epsilons.row(i) = v;
  int accflag = 0;
  for (int j=0;j<nmet;++j)  { // TODO must be a better way/vectorise?
    accflag += (epsilons(i,j)<=tolerance(j));
  }
  return(accflag==nmet); // 1 for all metrics satisfied
}


/** \brief Have all particles accepted, update state.
 * \return void
 */
void Smc::update()  {
  // Recalculate weights
  if (LFLAG) cout << "w" << flush;
  if (step==0)  {
    weights_new.fill(1.0);
  }
  else  {
    for (int i=0;i<N;++i)  {
      double w_sum = 0.0;
      for (int j=0;j<N;++j)  {
        w_sum += weights_old(j)*pert_dens(i,j);
      }
      weights_new[i] = priors(i)/w_sum; // Flat priors already adhered to in gen
      // FIXME Gamma distributed priors for within-farm. Where to specify requirement?
    }
  }
  // Normalise
  weights_old = weights_new/weights_new.sum();

  // Effective sample size
  if (LFLAG) cout << "e" << flush;
  ess = 1.0/weights_old.squaredNorm();

  // Calc mu and sigma for next round of perturbations
  if (LFLAG) cout << "m" << flush;
  means = particles_new.colwise().mean().transpose();
  if (LFLAG) cout << "s" << flush;
  MatrixXd centered = particles_new.rowwise() - means.transpose();
  if (CFLAG)  {
    sigma = (centered.adjoint() * centered) / double(N-1.0);
    sginv = sigma.inverse();
  }
  else  {
    sigma = ((centered.cwiseProduct(centered)).colwise().sum())/double(N-1.0);// Component-wise vars
    sginv = sigma.array().sqrt();
  }

  // Shrink tolerance
  if (LFLAG) cout << "t" << flush;
  for (int i=0;i<nmet;++i)  {
    VectorXd v = epsilons.col(i);
    std::sort(v.data(),v.data()+v.size());
    tolerance(i) = v[round(0.5*N)];         // TODO: median vs nth percentile efficiency?
  }

  // Particles
  particles_old = particles_new;
}


/** \brief Probability old particle j moves to where new i is now
 * \param i int
 * \param j int
 * \return double
 */
double Smc::pert_dens(int i, int j)  {
  if (CFLAG)  {
    VectorXd dist = (particles_new.row(i)-particles_old.row(j)).transpose()-means;
    VectorXd mtmp = dist.transpose()*sginv*dist;
    if (mtmp.size()!=1) cout << "mtmpsizewtf" << endl;
    return(exp(double(-0.5*mtmp(0))));
  }
  else  {
    double ker = 1.0;
    VectorXd dist = particles_new.row(i)-particles_old.row(j);
    for (int p=0;p<npar;++p)  {
      ker *= gsl_ran_gaussian_pdf(dist(p),sginv(p));
    }
    return(ker);
  }
}


void Smc::rngfree()  {
  for (unsigned int i=0;i<r.size();++i)  {
    gsl_rng_free(r[i]);
  }
}




