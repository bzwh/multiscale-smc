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
  models.resize(nchains,Model());
  for (int cno=0;cno<nchains;++cno)  {
    models[cno].setup(r[cno],frms);
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
  sig << means.transpose() <<"\n" << sigma << endl;
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
  while (pcheck)  {
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
    pcheck = models[cno].params.prior_check();
  }
}


/** \brief Generate perturbation using component-wise variances.
 * \param cno int
 * \return VectorXd
 */
VectorXd Smc::perturbation(int cno)  {
  /*VectorXd v(npar);
  for (int p=0;p<npar;++p)  {
    v(p) = gsl_ran_gaussian(r[cno],0.68*sigma(p));
  }
  return(v);*/
  return(mgen(sigma,r[cno]));
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
      weights_new[i] = models[0].params.prior_check() / w_sum;  // FIXME flat prior, so...
    }
  }
  // Normalise
  weights_old = weights_new/weights_new.sum();

  // Effective sample size
  if (LFLAG) cout << "e" << flush;
  ess = 1.0/weights_old.squaredNorm();
  // Calc sigma and shrink for next round of perturbations
  if (LFLAG) cout << "m" << flush;
  means = particles_new.colwise().mean();
  if (LFLAG) cout << "s" << flush;

/*  sigma.fill(0.0);
  for (int i=0;i<N;++i)  {  // (particles_new-means).sum() ???
    VectorXd diff = particles_new.row(i)-means.transpose();
    sigma += diff.cwiseProduct(diff);
  }
  sigma /= double(N);*/
  MatrixXd centered = particles_new.rowwise() - means;
  sigma = (centered.adjoint() * centered) / double(N-1.0);

  // Shrink tolerance
  if (LFLAG) cout << "t" << flush;
  for (int i=0;i<nmet;++i)  {
    //tolerance(i) = std::max(0.0,0.8*epsilons.col(i).maxCoeff());
    VectorXd v = epsilons.col(i);
    std::sort(v.data(),v.data()+v.size());
    tolerance(i) = v[round(0.8*N)];
  }
  //std::sort(epsilons.data(),epsilons.data()+N);
  //tolerance = epsilons(round(0.5*N));

  // Particles
  particles_old = particles_new;
  if (LFLAG) cout << endl;
}


/** \brief Probability old particle j moves to where new i is now
 * \param i int
 * \param j int
 * \return double
 */
double Smc::pert_dens(int i, int j)  {
  double ker = 1.0;
  VectorXd dist = particles_new.row(i)-particles_old.row(j);
  for (int p=0;p<npar;++p)  {
    ker *= gsl_ran_gaussian_pdf(dist(p),sigma(p));
  }
  return(ker);

  return(pow(6.2831853,-N*0.5) * sigma.determinant() * exp(double(-0.5*(dist-means)*sigma.inverse()*(dist-means))));
}


void Smc::rngfree()  {
  for (unsigned int i=0;i<r.size();++i)  {
    gsl_rng_free(r[i]);
  }
}




