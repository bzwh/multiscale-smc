#include <iostream>
#include <fstream>
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

using namespace Eigen;
using std::cout;
using std::endl;
using std::vector;
using std::flush;
using std::ofstream;
using std::to_string;
using std::pow;
using std::exp;

/** \brief Constructor.
 * \param nc int Number of threads.
 * \param nn int Number of particles.
 * \param ss int Number of rounds of filtering
 * \param frms Farms& Big chunk of const data, just pass through to Model (avoid reading repeatedly)
 */
Smc::Smc(int nc,int nn,int ss,Farms& frms)  {
  nchains = nc;
  r.resize(nchains);
  for (int cno=0;cno<nchains;++cno)  {
    r[cno] = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r[cno],(cno+1.0)*time(0));
  }
  N = nn;
  nround = ss;
  models.reserve(nchains);
  for (int cno=0;cno<nchains;++cno)  {
    models.push_back(Model(frms));
    models[cno].setup(r[cno]);
  }
  npar = models[0].params.pnum; // == particles[0].size()
  nmet = frms.enreg;                 // Number of metrics to test against...
  particles_old = MatrixXd::Zero(N,npar);
  particles_new = MatrixXd::Zero(N,npar);
  weights_old.resize(N);
  weights_new.resize(N);
  means.resize(npar);
  sigma = MatrixXd::Zero(npar,npar);
  tolerance.resize(nmet);
  tolerance.fill(0.2);
  epsilons = MatrixXd::Zero(N,nmet);
  step = 0;
}


/** \brief Loop rounds of filtering. Opens output files and calls write();
 * \return void
 */
void Smc::run()  {
  ofstream wht("./outputs/wht.txt");
  ofstream sig("./outputs/sig.txt");
  for (step=0;step<nround;++step)  {
    ofstream dat("./outputs/dat"+to_string(step)+".txt");
    ofstream eps("./outputs/err"+to_string(step)+".txt");
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
 * \param dat ofstream& For particles
 * \param eps ofstream& For metric values
 * \param wht ofstream& Particle weights
 * \param sig ofstream& Dump particle means and (component-wise) variances
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
  #pragma omp parallel for num_threads(nchains)
  for (int i=0;i<N;++i)  {
    int cno = omp_get_thread_num();
    cout << i << "-" << cno << ":" << flush;
    int acc = 0;
    while (acc==0)  {
      gen(i,cno);
      runmodel(i,cno);
      acc = test(i,cno);
    }
  }
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
      //VectorXd pb = perturbation(cno);
      //if (pb.size()!=particles_new.row(i).size())  {cout << "!!WTF!!"<<endl;}
      //particles_new.row(i) = particles_old.row(j) + pb.transpose();// CRASHED HERE??
      particles_new.row(i) = particles_old.row(j) + perturbation(cno).transpose();// CRASHED HERE??
    }
    models[cno].parse(particles_new.row(i));
    pcheck = models[cno].params.prior_check();  // 0 for all good
  }
}


/** \brief Generate perturbation using component-wise variances.
 * \param cno int
 * \return VectorXd
 */
VectorXd Smc::perturbation(int cno)  {
  return(mgen(0.68*sigma,r[cno]));
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
      weights_new[i] = 1.0/w_sum; // Flat priors already adhered to in gen
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
  sigma = (centered.adjoint() * centered) / double(N-1.0);
  sginv = sigma.inverse();

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
  //cout << endl << "sinv" << endl << sginv<< endl;
  //cout << endl << "mean" << endl << means<< endl;
  //cout << endl << "pars" << endl << particles_new.row(i)-particles_old.row(j)<< endl;
  VectorXd dist = (particles_new.row(i)-particles_old.row(j)).transpose()-means;
  //cout << endl << "dist" << endl << dist << endl;
  VectorXd mtmp = dist.transpose()*sginv*dist;
  if (mtmp.size()!=1) cout << "mtmpsizewtf" << endl;
  //cout << "mtmp:" << mtmp.size() << endl;
  return(exp(double(-0.5*mtmp(0))));
}


void Smc::rngfree()  {
  for (unsigned int i=0;i<r.size();++i)  {
    gsl_rng_free(r[i]);
  }
}




