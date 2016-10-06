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

// Logging
#define LFLAG 1
// Full covariance matrix vs component-wise
#define CFLAG 0


using namespace Eigen;
using std::cout;
using std::endl;
using std::vector;
using std::flush;
using std::ofstream;
using std::ifstream;
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
 * \param pf int Flag for continuing a run
 * \param frms Farms& Big chunk of const data, just pass through to Model (avoid reading repeatedly)
 */
Smc::Smc(int nc,int nn,int ss,int pf,Farms& frms)  {
  nthreads = nc;
  r.resize(nthreads);
  for (int cno=0;cno<nthreads;++cno)  {
    r[cno] = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r[cno],(cno+1.0)*time(0));
  }
  N = nn;
  nround = ss;
  pflag = pf;
  models.reserve(nthreads);
  for (int cno=0;cno<nthreads;++cno)  {
    models.push_back(Model(frms));
    models[cno].setup(r[cno]);
  }
  npar = models[0].params.pnum; // == particles[0].size()
  nmet = frms.enreg;                 // Number of metrics to test against...
  particles_old = MatrixXd::Zero(N,npar);
  particles_new = MatrixXd::Zero(N,npar);
  priors = VectorXd::Zero(N);
  weights_old = VectorXd::Zero(N);
  weights_new = VectorXd::Zero(N);
  means = VectorXd::Zero(npar);
  sigma = MatrixXd::Zero(npar,npar);
  tolerance.resize(nmet);
  tolerance.fill(2.2);      // FIXME
  //tolerance << 0.105,0.165,0.180,0.140,0.08;
  epsilons = MatrixXd::Zero(N,nmet);
  step = 0;
  ess = 0.0;
  n_done = 0;
  if (pflag)  {
    restart();
  }
}


void Smc::restart()  {
  if (LFLAG)  {
    cout << "...reloading:" << endl;
  }
  // To read in means sigma and sginv---------------------------------------------------------------
  ifstream sfile("./outputs/sig.txt");
  VectorXd tmp_in = Eigen::VectorXd::Zero(3*npar);
  int numrnds = 0;
  if (sfile.is_open())  {
    string line;
    while (getline(sfile,line))  {
      ++numrnds;
    }
    sfile.clear();
    sfile.seekg(0,sfile.beg);
    step = numrnds/4; // Becuase 3 lines for mu,sigma,sginv and a blank line separator...
    cout << "sig numlines: " << numrnds << " Step number: " << setfill('0') << setw(2) << step << endl;
    int to_skip = npar*(step-1)*3-1;    // Mean, variance, stdev...
    double read_in = 0.0;
    int isig = 0;
    while(sfile >> read_in)  {
      if (isig++==to_skip)  {
        break;           // Discard up to last round
      }
      if (sfile.eof())  {
        cout << "but shouldn't have reached end of sig file..." << endl;
        exit(-1);
      }
    }
    isig = 0;
    while (sfile >> read_in)  {
      if (sfile.eof())  {
          cout << "!" << endl;
        break;
      }
      tmp_in(isig++) = read_in;
    }
    if (isig!=3*npar)  {
      cout << "something went wrong? " << isig << endl;
      exit(-1);
    }
  }
  else  {
    cout << "Want to restart but no sig" << endl;
    exit(-1);
  }
  sfile.close();
  means = tmp_in.head(npar);
  sigma = tmp_in.segment(npar,npar);
  sginv = tmp_in.tail(npar);
  if (LFLAG)  {
    cout << " mu:\t" << means.transpose() << endl
         << " sg:\t" << sigma.transpose() << endl
         << " si:\t" << sginv.transpose() << endl;
  }

  // To read in particles_old from datXY.txt and N particles----------------------------------------
  stringstream pname;
  string line;
  int numlines = 0;
  ifstream pdat;
  pname << "./outputs/dat" << setfill('0') << setw(2) << to_string(step-1) << ".txt";
  pdat.open(pname.str());
  if (pdat.is_open())  {
    while (getline(pdat,line))  {
      ++numlines;
    }
    pdat.clear();
    pdat.seekg(0,pdat.beg);
    for (int iline=0;iline<N;++iline)  {
      Eigen::VectorXd ptmp = Eigen::VectorXd::Zero(npar);
      for (int ipar=0;ipar<npar;++ipar)  {
        pdat >> ptmp(ipar);
      }
      particles_old.row(iline) = ptmp;
      if (pdat.eof())  {
        cout << "trying to read too much from " << pname.str() << endl;
        exit(-1);
      }
    }
  }
  else  {
    cout << "can't open old " << pname.str() << " file" << endl;
    exit(-1);
  }
  N = numlines;     // This better be right!
  // Happily read everything expected, but what if not at end of file? ie more particles/parameters?
  // Check model spec'd? Can't really guess Settings.cpp and Country from outputs?

  // Reading particles done so far------------------------------------------------------------------
  VectorXd tmp_par = VectorXd::Zero(npar);
  VectorXd tmp_err = VectorXd::Zero(nmet);
  ifstream dfile("./outputs/dtmp.txt");
  if (dfile.is_open())  {
    while(true)  {
      for (int ipar=0;ipar<npar;++ipar)  {    // Grab complete particle
        dfile >> tmp_par(ipar);
      }
      for (int ierr=0;ierr<nmet;++ierr)  {    // And all nmet epsilon values
        dfile >> tmp_err(ierr);
      }
      if (dfile.eof())  {                     // FIXME tries to read extra npar+nmet times...
        break;
      }
      particles_new.row(n_done) = tmp_par;
      tmp_par.fill(0);
      epsilons.row(n_done) = tmp_err;
      ++n_done;
      tmp_err.fill(0);
    }
  }
  else  {
    cout << "Want to restart but no data" << endl;
    exit(-1);
  }
  dfile.close();
  #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
  for (int ii=0;ii<n_done;++ii)  {
    int cno = omp_get_thread_num();
    models[cno].parse(particles_new.row(ii));
    models[cno].params.prior_check();  // 0 for all good
    priors(ii) = models[cno].params.pri;
  }

  // To read in particle WEIGHTS and STEP-----------------------------------------------------------
  ifstream wfile("./outputs/wht.txt");
  VectorXd tmp_wht = VectorXd::Zero(N);
  if (wfile.is_open())  {
    int to_skip = (N+1)*(step-1)-1;   // Have ESS + N weights per line. now on step+1'th round (line)
    double read_in = 0.0;
    int iwht = 0;         // discard. seekg would be better...?
    while (wfile>>read_in)  {
      if (iwht==to_skip)  {
        break;
      }
      if(wfile.eof())  {
        cout << "wtf, this should all be previous rounds' weights" << endl;
        exit(-1);
      }
      ++iwht;
    }
    wfile >> ess;         // now at last line
    iwht = 0;
    while(wfile >> read_in)  {
      tmp_wht(iwht++) = read_in;
    }
    if (iwht!=N)  {
      cout << "Not go the right number of weights..." << iwht << ". ESS: " <<  ess << endl;
      cout << "Weights read in so far... " << tmp_wht.transpose() << endl;
      exit(-1);
    }
  }
  else  {
    cout << "Want to restart but no weights" << endl;
    exit(-1);
  }
  weights_old = tmp_wht;
  wfile.close();
  if (LFLAG)  {
    cout << " ESS:\t" << ess << endl;
    cout << " WHT:\t" << weights_old.head(1) << "\t" << weights_old.tail(1) << endl << endl;
  }

  // To read in tolerance to use from errXY.txt
  stringstream ename;
  ifstream edat;
  ename << "./outputs/err" << setfill('0') << setw(2) << to_string(step-1) << ".txt";
  Eigen::VectorXd etmp = Eigen::VectorXd::Zero(nmet);
  edat.open(ename.str());
  if (edat.is_open())  {
    for (int imet=0;imet<nmet;++imet)  {
      edat >> etmp(imet);
    }
    if (edat.eof())  {
      cout << "tring to read too much from errxy.txt" << endl;
      exit(-1);
    }
  }
  else  {
    cout << "can't open old errxy.txt file" << endl;
    exit(-1);
  }
  tolerance = etmp;
  if (LFLAG)  {
    cout << " Tolerance:\t" << tolerance.transpose() << endl;
  }

  cout << " Continuing step: " << setfill('0') << setw(2) << step << endl;
  cout << " from particle " << n_done << "/" << N << endl;
  std::cin.get();
}


/** \brief Loop rounds of filtering. Opens output files and calls write();
 * \return void
 */
void Smc::run()  {
  ofstream wht;
  ofstream sig;
  if (pflag)  {
    wht.open("./outputs/wht.txt",ofstream::app);
    sig.open("./outputs/sig.txt",ofstream::app);
  }
  else  {
    wht.open("./outputs/wht.txt");
    sig.open("./outputs/sig.txt");
  }

  while (step<nround)  {
    // Setting up this round's particle data and metric files
    stringstream d_tmp;
    d_tmp  << "./outputs/dat" << setfill('0') << setw(2) << to_string(step) << ".txt";
    ofstream dat(d_tmp.str());
    stringstream e_tmp;
    e_tmp  << "./outputs/err" << setfill('0') << setw(2) << to_string(step) << ".txt";
    ofstream eps(e_tmp.str());
    p_dat.open("./outputs/dtmp.txt",ofstream::app);
    cout << "ROUND: "<<step << endl;
    iterate();
    update();
    cout << "\n Tolerance: " << tolerance.transpose() << "\n" << endl;
    write(dat,eps,wht,sig);
    dat.close();
    eps.close();
    ++step;
    p_dat.close();
    p_dat.open("./outputs/dtmp.txt"); // Must be a better way to close and clear file...
    p_dat.close();
  }
  wht.close();
  sig.close();
  rngfree();
}


/** \brief Dump state after each round.
 * Note: called after update(), so tolerance at head of errXY.txt will be the new criteria,
 *  not what was used to judge the error values in file.
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
  sig << means.transpose() <<"\n" << sigma << "\n" << sginv << "\n\n" << flush;
}


/** \brief Looping through all particles in this round. Spawns however many threads to run through.
 * \return void
 */
void Smc::iterate()  {
  #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
  for (int i=n_done;i<N;++i)  {
    stringstream lstream;
    stringstream pstream;
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
    pstream << particles_new.row(i) << "\t" << epsilons.row(i) << "\n";
    p_dat << pstream.str() << flush;
    pstream.str(string());
    lstream.str(string());
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
      particles_new.row(i) = particles_old.row(j) + perturbation(cno).transpose();
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
  n_done = 0;
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




