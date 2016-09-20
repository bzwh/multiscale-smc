#include <cmath>
#include <functional>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "Model.hpp"
#include "Farm.hpp"
#include "Grid.hpp"
#include "Params.hpp"
#include "Settings.hpp"

#define LFLAG 0

using namespace std;

// TODO vaccination routine and parameters
// TODO JPN CPs

Model::Model(Farms& frms) : farms(frms) {
  tmax = G_CONST::simtmax; // 242 runs 01 Feb - 01 Oct.
  twmx = G_CONST::withtmx;
}


/** \brief Initialise Model. Holds basically everything...
 * Loads up farm data etc. sets up vectors for states and transitions, blahblahblah
 * \param rr gsl_rng* pointer for all your randoming needs
 */
void Model::setup(gsl_rng* rr)  {
  r = rr;
  N = farms.ncnt;
  S = N;
  I = 0;  II.resize(tmax,0);
  Itot = 0;
  C_d.resize(farms.nreg,0);
  R_d.resize(farms.nreg,0);
  RIP = 0;
  RDC = 0;
  RCP = 0;
  t = 1;  // This is 1st Feb
  states.resize(N,0);
  dstate.resize(N,0);
  time_i.resize(N,-1);
  time_r.resize(N,-1);
  time_c.resize(N,-1);
  i_by.resize(N,-1);
  sus.resize(N,0.0);
  tblty.resize(N,vector<double>(0,0.0));
  // Init the grids - requires farm locations
  grid.setup(farms);
  igrids.resize(grid.ngrids,0);// Total infection counts in each grid over all runs
  // Track number of cases & culls for farms and animals
  ninfd.resize(farms.enreg,0);  // Cases
  nculd.resize(farms.enreg,0);  // Culls
  infd = vector< vector<int> >(farms.enreg,vector<int>(3,0)); // Cases in 5 #s, 3 species
  culd = vector< vector<int> >(farms.enreg,vector<int>(3,0)); // Culls in 5 regions, 3 species
  ertot.resize(farms.enreg, vector<double> (6,0.0));

  // TODO move out of constructor - belongs with farm data?
  norm_fi.resize(farms.enreg,0.0); // Normalising constants for each region
  norm_fc.resize(farms.enreg,0.0); // TODO could store species together
  norm_ci.resize(farms.enreg,0.0);
  norm_si.resize(farms.enreg,0.0);
  norm_cc.resize(farms.enreg,0.0);
  norm_sc.resize(farms.enreg,0.0);
  if (G_CONST::err_d_c)  { // These are daily values, want to normalise to cumulative end
   for (int reg=0;reg<farms.enreg;++reg)  {
      norm_fi[reg] = 1.0/max((double)( accumulate(farms.farmsi[reg].begin(),farms.farmsi[reg].begin()+tmax,0) ),1.0);
      norm_fc[reg] = 1.0/max((double)( accumulate(farms.farmsc[reg].begin(),farms.farmsc[reg].begin()+tmax,0) ),1.0);
      norm_ci[reg] = 1.0/max((double)( accumulate(farms.cowssi[reg].begin(),farms.cowssi[reg].begin()+tmax,0) ),1.0);
      norm_si[reg] = 1.0/max((double)( accumulate(farms.sheepi[reg].begin(),farms.sheepi[reg].begin()+tmax,0) ),1.0);
      norm_cc[reg] = 1.0/max((double)( accumulate(farms.cowssc[reg].begin(),farms.cowssc[reg].begin()+tmax,0) ),1.0);
      norm_sc[reg] = 1.0/max((double)( accumulate(farms.sheepc[reg].begin(),farms.sheepc[reg].begin()+tmax,0) ),1.0);
    }
  }
  else  { // Already cumulative totals, so normalise to last value
    for (int reg=0;reg<farms.enreg;++reg)  {
      norm_fi[reg] = 1.0/max((double)(farms.farmsi[reg][tmax]),1.0);
      norm_fc[reg] = 1.0/max((double)(farms.farmsc[reg][tmax]),1.0);
      norm_ci[reg] = 1.0/max((double)(farms.cowssi[reg][tmax]),1.0);
      norm_si[reg] = 1.0/max((double)(farms.sheepi[reg][tmax]),1.0);
      norm_cc[reg] = 1.0/max((double)(farms.cowssc[reg][tmax]),1.0);
      norm_sc[reg] = 1.0/max((double)(farms.sheepc[reg][tmax]),1.0);
    }
  }

}


void Model::setrdat(int cno,string fpath)  {
  stringstream fn;
  fn << fpath << "rundat_" << cno;
  //rdat.open(fn.str().c_str());
  //if (!rdat.is_open())  {
  //    cout << "can't open rundat file" << endl;
  //    exit(-1);
  //
}


/** \brief Initialise replicate with generated parameter set
 * Has a little thing for the kernel, to avoid sqrt'ing
 * \return void
 */
void Model::initrun()  {
  plaw = params.plaw;  // Using linear scaling or powerlaws on sus, trans?
  pdet = params.pdet;
  // NOTE old ker[0] absorbed
  params.ker[0] = 1.0/(params.ker[0]*params.ker[0]);
  params.ker[1] = params.ker[1]*0.5; // kernel uses d^2 to avoid the sqrt
  sus_calc();
  maxrate_calc();
  seed_run();
}


/** \brief Caclulate susceptibilities of all farms
 * Already have stored farm data and parameter values
 * \return void - values all stored in member vector
 */
void Model::sus_calc()  {
  // FIXME HACKED parameter vs farm region. should be 0 anyway, only eregion[i] is a thing
  if (plaw)  {
    for (int i=0;i<N;++i)  {
      int reg = farms.region[i];
      sus[i] = params.sb[reg][0]*pow(farms.N[i][0],params.sp[reg][0])
             + params.sb[reg][1]*pow(farms.N[i][1],params.sp[reg][1]);
    }
  }
  else  {
    for (int i=0;i<N;++i)  {
      int reg = farms.region[i];
      sus[i] = params.sb[reg][0]*farms.N[i][0]
             + params.sb[reg][1]*farms.N[i][1];
    }
  }
}


/** \brief Calculate max transmission rate between grids = maxsus[h]*dker(gdist(g,h))
 * \return void - stored in Grid.mrate
 */
void Model::maxrate_calc()  {
  // First get max susceptibility of each grid
  vector<double> maxsus(grid.ngrids,0.0);
  for (int g=0;g<grid.ngrids;++g)  {
    for (auto currfarm=grid.f_id[g].begin();currfarm!=grid.f_id[g].end();++currfarm)  {
      int i = *currfarm;
      maxsus[g] = max(maxsus[g],sus[i]);
    }
    /*if (maxsus[g]>0.0)
      cout << g << " " << maxsus[g] << endl;*/
  };
  // Now get maxrates
  for (int g=0;g<grid.ngrids;++g)  {
    auto h_it = grid.potgrids[g].begin();
    auto m_it = grid.mrate[g].begin();
    while (h_it!=grid.potgrids[g].end())  {
      *m_it = maxsus[*h_it]*dker(grid.gdist(g,*h_it));
      ++h_it;
      ++m_it;
    }
  }
}


/** \brief Calculate the DC cull intensity at the current simulation time
 *
 */
void Model::dccullF_calc()  { // too many overheads?
  if (G_CONST::fit_dcs==2)  {
    dcF = params.F[0]*t*exp(-params.F[1]*t);
  }
  else if (G_CONST::fit_dcs==1)  {
    dcF = params.F[0];
  }
  else  {
    dcF = 50.0;
  }
}


/** \brief Seed infection, determined by the data read in earlier (which_farms,cull_farms)
 * \return void
 */
void Model::seed_run()  {
  auto sf_it = farms.seedfarms.begin(); // farm to seed iterator
  auto st_it = farms.seedtimes.begin(); // time to seed iterator
  auto end_seed = farms.seedtimes.end();
  auto cf_it = farms.seedculls.begin(); // farm to cull
  auto ct_it = farms.seedtimec.begin(); // farm to seed
  auto end_cull = farms.seedtimec.end();
  // here tt is temporary date, NOT seedtimes index.
  for (int tt=1;tt<=farms.seedtimes.back();++tt)  {
    t = tt;
    // any more farms to be seeded? any today? stupid way to loop through these
    while ((st_it!=end_seed)&&(*st_it<=tt))  { // check order: check value exists before reading
      dstate[*sf_it] = 1;
      ++sf_it;
      ++st_it;
      ++Itot;
    }
    // any more farms to cull? any today?
    while ((ct_it!=end_cull)&&(*ct_it<=tt))  { // check order: check value exists before reading
      dstate[*cf_it] = 6;
      ++cf_it;
      ++ct_it;
    }
    //  dont wan't detrepcul(), no detections! so just manually check...
    for (int i=0;i<N;++i)  {
      if (t==time_r[i])  {
        if (states[i]==1)  {
          // Detected IP... BUT NOT REALLY. just setting up for imminent cull
          dstate[i] = 2;
        }
      }
      else if (t==time_c[i])  { // To be culled. switch-case prettier?
        if (states[i]==2)  {
          dstate[i] = 5;
        }
        else if (states[i]==3)  {
          dstate[i] = 6;
        }
      }
    }
    error_calc();
    //cout << t<< " " << ertot[0] << " " << ertot[1] << " " << ertot[2] << " " << ertot[3] << " " << ertot[4] << endl;
    update();
  }
}


/** \brief Run simulation with proposed parameters on seeded state.
 * \return 0=> going far too big. 1=> outbreak plausible enough to be tested by fit metric
 */
int Model::runsim() {
  initrun();
  while(t<tmax)  {
    iterate();
    detrepcul();
    error_calc();
    update();
    ++t;
    if (Itot>4000)  { // early termination when obv too big
      for (int ii=0;ii<farms.enreg;++ii)  {
        fill(ertot[ii].begin(),ertot[ii].end(),1.0e10);
      }
      //cout << "B" << flush;
      return(0);  // reject
    }
  }
  //cout << "ITOT(" << Itot << ")" << endl;
  return(1);      // accept
}


/** \brief Calc today's infection events. Stick em as dstate[j]=1 for Model::update() to work on.
 * Loop through infectious farms, testing transmission to nearby grids (grid.potgrids[g]) first.
 * Then test against all susceptibles in that grid.
 */
void Model::iterate() {
  // Each farm
  for (int i=0;i<N;++i)  {
    // Check infectious farms
    if ((states[i]==1)&&trans(i))  {  // FIXME trans(i) potentially not allocated!! tblty[i].size()?
      int g = grid.g_id[i]; // infectious farm's grid number
      // First test transmission to susceptible farms in same grid g
      for (auto sfarm=grid.f_id[g].begin();sfarm!=grid.f_id[g].end();++sfarm)  {
        int j = *sfarm;
        // TODO check dstate as well??
        if (states[j]==0)  { // Only target susceptibles, not changing this day?
          double q = 1.0 - exp(-sus[j]*trans(i)*dker(dist2(i,j)));
          if (gsl_rng_uniform(r)<q)  {
            dstate[j] = 1;  // INFECT
          }
        }
      }
      // Now test transmissions to near-enough non-empty grids (def by potgrids-excludes self)
      //for (auto hh=grid.potgrids[g].begin();hh!=grid.potgrids[g].end();++hh)  {
      auto hh_it = grid.potgrids[g].begin();
      auto mr_it = grid.mrate[g].begin();
      while (hh_it!=grid.potgrids[g].end())  {
        int h = *hh_it; // potgrid number
        double transi = trans(i);
        double epab = exp(-transi* (*mr_it));
        double epa2b = gsl_sf_pow_int(epab,grid.left[h]); // 1-P(A to B)
        // Gets in to grid?
        if (gsl_rng_uniform(r)<(1.0-epa2b))  {
          int s = 1; // none in grid h infected by farm i yet
          double invepab = 1.0/epab;
          for (auto sfarm=grid.f_id[h].begin();sfarm!=grid.f_id[h].end();++sfarm)  {
            // farms j in grid h
            int j = *sfarm;            // susceptible farm
            if (states[j])  { // TODO check dpsread? only susceptible ...
              continue;
            }
            else  {
              // TODO farms infected during this iteration? multiple infecterisers?
              double rnd = gsl_rng_uniform(r);
              if (s)  {
                rnd *= (1.0 - epa2b); // P = 1-s(1-P_ab)**Nb
                epa2b *= invepab;     // for the next iteration, P = 1-(oldepa2b/epab)
              }
              // j infected assuming prob is the maximum value P_ab
              if (rnd<(1.0-epab))  {
              // Rnd < P_ab/P
                s = 0;
                //if (dstate[j]) {continue;} // already infected this iteration, maybe doesn't save much
                // prob jth susc in h infected by i
                double q = 1.0 - exp(-sus[j]*transi*dker(dist2(i,j)));
                if (gsl_rng_uniform(r)<q)  {
                  dstate[j] = 1;
                }
              }
            }
          }
        }
        ++hh_it;
        ++mr_it;
      }
    }
  }
}


/** \brief Today's detections, reports and culls.
 * Detection: dstate=2 and need to identify DCs and CPs
 * Report DC/CP: dstate=3 add some delay and then cull
 * Culls: dstate=5 IPcull and dstate=6 pre-emptive cull
 */
void Model::detrepcul()  {
  // FIXME comparing all time_r/time_c heavy on DLmr and branch misses. What do?
  //    better data structure?
  // 50% of runsim() DLmr
  dccullF_calc();
  for (int i=0;i<N;++i)  {
    if (t==time_r[i])  {  // Detection or report
      switch (states[i])  {
        case 1: // Detected IP
          dstate[i] = 2;
          if (G_CONST::dc_cull)  {
            dccull(i);
          }
          if (G_CONST::cp_cull)  {
            cpcull(i);
          }
        break;

        case 0: // Reporting DC/CP
          dstate[i] = 3;
        break;

        default:
          //cout << t << " " << i << " " << states[i] << " " << dstate[i] << " DETREPCUL" << endl; // FIXME Always borks at t=23
          //cin.get();
        break;
      }
    }
    else if (t==time_c[i])  { // Cull
      switch (states[i])  {
        case 2: // Culling IP
          dstate[i] = 5;
        break;

        case 3: // Culling DC/CP
          dstate[i] = 6;
        break;
      }
    }
  }
}


/** \brief Carry out whatever happened today- updating states according to dstate
 */
void Model::update()  {
  for (int i=0;i<N;++i)  {
    switch (dstate[i])  {
      case 0:
      break;

      case 1: // Infection
        states[i] = 1;
        time_i[i] = t;
        (params.rwfm) ? within(i) : without(i); // std::function  functor? function pointer? OTT?
        ++I; ++II[t]; ++Itot;
        --S;
        ++igrids[grid.g_id[i]];
      break;

      case 2: // Detected IP
        states[i] = 2;
        time_r[i] = t;
        ++R_d[farms.region[i]];
      break;

      case 3: // Reporting DC/CP
        states[i] = 3;
        time_r[i] = t;
        time_c[i] = (time_c[i] > 0) ? min(time_c[i],t+params.delaydcp) : t+params.delaydcp;
      break;

      case 4: // Reporting occult
        states[i] = 3; // TODO Occult infections?
        time_r[i] = t;
        time_c[i] = (time_c[i] > 0) ? min(time_c[i],t+params.delaydcp) : t+params.delaydcp;
      break;

      case 5: // Culling IP
        states[i] = -1;
        time_c[i] = t;
        --I;
        ++RIP;
      break;

      case 6: // Culling DC/CP
        states[i] = -2;
        time_c[i] = t;
        --S;
        //++R2;
        ++C_d[farms.region[i]];
      break;
    }
  }
  fill(dstate.begin(),dstate.end(),0);
  fill(C_d.begin(),C_d.end(),0);
  fill(R_d.begin(),R_d.end(),0); // This is a daily counter...
}


/** \brief Identify DCs of the detected infectious farm i. dstate[DC] = 3
 * \param ifarm int - infectious farm
 */
void Model::dccull(int ifarm)  { // TODO grid this? Probably not worth it
  double prob = 0.0;
  int g = grid.g_id[ifarm];
  for (auto hh=grid.potgrids[g].begin();hh!=grid.potgrids[g].end();++hh)  {
    // Each near enough grid [h]
    int h = *hh;
    for (auto jj_it=grid.f_id[h].begin();jj_it!=grid.f_id[h].end();++jj_it)  {
      // Each farm [j] in grid [h]
      int j = *jj_it;
      if (states[j]==0)  { // TODO check dstate here?
        prob = exp(-dcF*sus[j]*trans(ifarm)*dker(dist2(ifarm,j)));
        if (i_by[j]==ifarm)  {
          prob *= params.f;
        }
        if (gsl_rng_uniform(r)<(1.0-prob))  {
          dstate[j] = 3;
          ++RDC;
        }
      }
    }
  }
}


/** \brief Get detected infectious farm i's CPs reported for culls dstate[CP]=3
 * \param i int - infectious farm
 */
void Model::cpcull(int i)  {
  for (auto jj=farms.cps[i].begin();jj!=farms.cps[i].end();++jj)  {
    int j = *jj;
    if (states[j]==0)  {  // TODO check dstate? susceptibles only?
      dstate[j] = 3;
      ++RCP;
    }
  }
}


/** \brief Distance^2 between farms i and j
 * \param i int - farm id
 * \param j int - farm id
 * \return double - The distance squared between them.
 */
double Model::dist2(int i, int j)  {
  double dx = farms.x[i] - farms.x[j];
  double dy = farms.y[i] - farms.y[j];
  if ((dx>15.0)||(dy>15.0))  {
    // FIXME 15km isn't nec safe, depends on kernel parameters. but probably fine. also see Grid
    return(-1.0);
  }
  else  {
    return( dx*dx+dy*dy );
  }
}


/** \brief Transmission kernel ( Euclidean distance^2 ) K = k0 / [1+(d/d0)^a]
 * \param d2 double
 * \return double
 */
double Model::dker(double d2)  {
  if (d2<0.0)  {
    return(0.0);
  }
  else  {
    // Expecting parameters rescaled to avoid having to use d=sqrt(...)
    return(1.0/(1.0+pow(d2*params.ker[0],params.ker[1])));
  }
}


/** \brief Without farm model. Generating step-function profile for farm i
 * \param ifarm int
 * \return void
 */
void Model::without(int ifarm)  {
  double tt = 0.0; // temp calc of transmissibility
  int delay = 0;    // delay to detection after infectiousness(==
  int latent = 5;
  switch (pdet)  {
    case 1:
      delay = params.ddet[0];
      break;
    case 2:
      delay = ceil(gsl_ran_gamma(r,params.ddet[0]/params.ddet[1],params.ddet[1])); // 0:mean 1:scale
      break;
    default:
      delay = params.delaydet;
      break;
  }
  time_r[ifarm] = t + latent + delay; // reported = infection(epi) + latent(fixed) + delay(fit)
  time_c[ifarm] = t + latent + delay + params.delayipc; // culled = reported+more_delay(fixed)

  int reg = farms.region[ifarm];
  if (G_CONST::fit_plw)  {
    tt = params.tb[reg][0]*pow(farms.N[ifarm][0],params.tp[reg][0])
       + params.tb[reg][1]*pow(farms.N[ifarm][1],params.tp[reg][1])
       + params.tb[reg][2]*pow(farms.N[ifarm][2],params.tp[reg][2]);
  }
  else  {
    tt = params.tb[reg][0]*farms.N[ifarm][0]
       + params.tb[reg][1]*farms.N[ifarm][1];
  }
  // a huge time difference with this unrolled loop. but not as general re #species? jsut leave all in and suppress params

  /*for (int k=0;k<latent;++k)  {
    tblty[ifarm].push_back(0.0);
  }
  for (int k=latent;k<20;++k)  {
    tblty[ifarm].push_back(tt);
  }*/
  tblty[ifarm].resize(1+latent + delay + params.delayipc,0.0);
  fill(tblty[ifarm].begin()+latent,tblty[ifarm].end(),tt);
}

/**  \brief Within-farm SEmInR model. Grab parameters from read-in files.
 */
void Model::within(int i)  {
  int HACKFLAG = 1;
  Eigen::Vector2d ptot = Eigen::Vector2d::Zero();
  Eigen::Vector2d n1pop(1.0,1.0);  // Inverse of population...
  ptot << farms.N[i][0],farms.N[i][1];
  n1pop << farms.N1[i][0],farms.N1[i][1];
  if (ptot.sum()==0)  {  // How the hell did it get infected then???
    time_r[i] = t+1;
    time_c[i] = t+1;
    if (i!=66481)  {// belongs to seed farms for the UK outbreak
      cout << "empty farm " << i << " got infected...?" << endl;
      exit(-1);
    }
  }
  else  {
    while(HACKFLAG)  {
      double tt = 0.0;                    // Time post-infection
      double tau = 0.1;                   // Discrete time steps
      double disp_r = 1.12;
      Eigen::Vector2d seed_mu = Eigen::Vector2d::Zero();
      Eigen::Vector2i seeds = Eigen::Vector2i::Zero();
      Eigen::Vector2i m_spc = Eigen::Vector2i::Zero();
      Eigen::Vector2i n_spc = Eigen::Vector2i::Zero();
      Eigen::Vector2d sigma_spc = Eigen::Vector2d::Zero();
      Eigen::Vector2d gamma_spc = Eigen::Vector2d::Zero();
      Eigen::Vector2d nsus = Eigen::Vector2d::Zero();
      Eigen::Vector2d ninf = Eigen::Vector2d::Zero();
      for (int spc=0;spc<ptot.rows();++spc)  {
        n1pop(spc) = (ptot(spc)) ? 1.0/ptot(spc) : 1.0;
        // Initial infections - see Schley ModelInf09 effect vacc undetected persistence
        // FIXME they're only for cows... but what to do with sheep?
        seed_mu(spc) = 4.21+0.00386*ptot(spc);
        while (seeds(spc)==0)  {
          seeds(spc) = gsl_ran_negative_binomial(r,1.0-seed_mu(spc)/(disp_r+seed_mu(spc)),disp_r);
        }
        m_spc(spc) = floor(params.kE[spc]);
        n_spc(spc) = max(int(floor(params.kI[spc])),1);
        sigma_spc(spc) = m_spc(spc)/params.mE[spc];
        gamma_spc(spc) = n_spc(spc)/params.mI[spc];
      }
      Eigen::Matrix2d beta;
      beta << params.bt[0][0] , 0.06,  // 0.06 from de Rueda
              0.06 , params.bt[1][1];

      if (LFLAG)  {
        cout << "\n";
        cout << " a: " << ptot.transpose() << endl;
        cout << " m: " << m_spc.transpose() << "\n n: " << n_spc.transpose() << "\n";
        cout << " s: " << sigma_spc.transpose() << "\n g: " << gamma_spc.transpose() << "\n";
        cout << " b: " << beta << "\n";
        cout << " I: " << seeds.transpose() << "\n";
      }

      // Populations - full vectors
      vector<Eigen::VectorXi> pops(2);
      for (int spc=0;spc<2;++spc)  {
        pops[spc] = Eigen::VectorXi::Zero(m_spc(spc)+n_spc(spc)+2);
        int ptmp = ptot(spc);
        pops[spc](1) = min(ptmp,seeds(spc));
        pops[spc](0) = ptot(spc)-pops[spc](1);
        nsus(spc) = pops[spc](0);
        ninf(spc) = pops[spc].segment(m_spc(spc)+1,n_spc(spc)).sum();
      }
      // rates of change
      vector<Eigen::VectorXd> rate(2);
      rate[0] = Eigen::VectorXd::Zero(m_spc(0)+n_spc(0)+2);
      rate[1] = Eigen::VectorXd::Zero(m_spc(1)+n_spc(1)+2);
      // population changes per time-step
      vector<Eigen::VectorXi> dpop(2);
      dpop[0] = Eigen::VectorXi::Zero(m_spc(0)+n_spc(0)+2);
      dpop[1] = Eigen::VectorXi::Zero(m_spc(1)+n_spc(1)+2);
      //cout << dpop[0].transpose() << endl << dpop[1].transpose() << endl;
      int ind = 1;        // day post infection - also index for tblty vector
      int t_wend = twmx-1;    // Time end of within-farm simulation
      vector<double> tb_tmp(twmx,0.0);
      int iflag = 1;      // flags detection status - goes to zero when time_r recorded
      int delay = 0;
      switch (pdet)  {
        case 1:
          delay = params.ddet[0];
          break;
        case 2:
          delay = ceil(gsl_ran_gamma(r,params.ddet[0]/params.ddet[1],params.ddet[1])); // 0:mean 1:scale
          break;
        default:
          delay = params.delaydet;
          break;
      }

      // Actual loop here
      while(ind<t_wend)  {
        if ((iflag)&&( ninf.sum()>0 ))  { // FIXME threshold triggers detection delay
          time_r[i] = t+ind+delay;        // So
          time_c[i] = t+ind+delay+params.delayipc;
          t_wend = ind+delay+params.delayipc;   // Shrink simulation end point
          if (t_wend>tb_tmp.size())  {tb_tmp.resize(t_wend,0.0);}
          iflag = 0;
        }
        // Calc rates and pop changes
        Eigen::Vector2d rt = tau*nsus.cwiseProduct(beta*ninf.cwiseProduct(n1pop)); // t*S*beta*I/N
        rate[0](0) = rt(0);
        rate[1](0) = rt(1);
        for (int spc=0;spc<2;++spc)  {
          if (ptot(spc))  {
            // Latent periods
            for (int ii=1;ii<m_spc(spc)+1;++ii)  {
              rate[spc](ii) = pops[spc](ii)*sigma_spc(spc)*tau;
            }
            // Infectious periods
            for (int ii=m_spc(spc)+1;ii<m_spc(spc)+n_spc(spc)+1;++ii)  {
              rate[spc](ii) = pops[spc](ii)*gamma_spc(spc)*tau;
            }
            // Population changes ~Pois(rates)
            for (int ii=0;ii<m_spc(spc)+n_spc(spc)+1;++ii)  {
              dpop[spc](ii) = floor(min(pops[spc](ii),int(gsl_ran_poisson(r,rate[spc](ii)))));
            }
          }
        }
        // Update!
        for (int spc=0;spc<2;++spc)  {
          pops[spc] = pops[spc] - dpop[spc];
          for (int ii=1;ii<m_spc(spc)+n_spc(spc)+2;++ii)  {
            pops[spc](ii) = pops[spc](ii) + dpop[spc](ii-1);
          }
        }
        nsus(0) = pops[0](0);
        nsus(1) = pops[1](0);
        ninf(0) = pops[0].segment(m_spc(0)+1,n_spc(0)).sum();
        ninf(1) = pops[1].segment(m_spc(1)+1,n_spc(1)).sum();
        tt += tau;
        // Re-bin for daily steps
        if (trunc(tt)>ind)  {
          tb_tmp[ind] = params.tb[farms.region[i]][0]*ninf(0) + params.tb[farms.region[i]][1]*ninf(1);
          ++ind;
        }
      }
      tblty[i] = tb_tmp;
      if (iflag)  {
        /*cout << "!";
        cout << t << " " << delay << " " << time_r[i] << " " << time_c[i] << " " << ind << endl;
        cout << ptot(0) << "\t" << seeds(0) << "\t" << m_spc(0) << "\t" << sigma_spc(0) << "\t" << n_spc(0) << "\t" << gamma_spc(0) << endl;
        cout << ptot(1) << "\t" << seeds(1) << "\t" << m_spc(1) << "\t" << sigma_spc(1) << "\t" << n_spc(1) << "\t" << gamma_spc(1) << endl << endl;
        cout << pops[0].transpose() << endl;
        cout << pops[1].transpose() << endl;
        cout << beta << endl;
        cin.get();*/
      }
      else  {
        HACKFLAG = 0;
      }
    }
  }
}


/** \brief Get transmissibility of farm i at time t
 * \param i int
 * \return double
 */
inline
double Model::trans(int i)  {
  if (t>time_c[i])  {
    cout << t << " " << i << " " << time_i[i] << " " << time_r[i] << " " << time_c[i] << " " << states[i] << " " << farms.N[i][0] << " " << farms.N[i][1] << " siadfhaksjm" << endl;
    exit(-1);
  }
  else  {
    return(tblty[i][t-time_i[i]]);
  }
}


/** \brief Dumping state of the epidemic daily- current #S #I culled IPs and culled DC/CPs
 * \return void
 */
void Model::dumpdailyir()  {
// Generates pretty not small files
  for (int tt=0;tt<tmax;++tt)  {
    //rdat << tt << " " << II[tt] << "\n";  // regional breakdown. easy to sum afterwards
  }
  //rdat.flush(); // Avoid partially written runs if terminated early
}


/**
 *  \brief Dump state of every farm with infection and cull times
 */
void Model::dumpfarmids()  {
  // Generates frickin massive files
  for (int i=0;i<N;++i)  {
    //rdat << i << " " << time_i[i] << " " << time_c[i] << " " << states[i] << " "<< grid.g_id[i] << "\n";
  }
  //rdat.flush();
}

void Model::dumpoutbreak(std::ofstream& whichdat,std::ofstream& culldat)  {
  for (int i=0;i<N;++i)  {
    switch(states[i])  {
      case -1:
        whichdat << i << " " << time_i[i] << "\n";
        break;

      case -2:
        culldat << i << " " << time_c[i
        ] << "\n";
        break;

      default:
        break;
    }
  }
  whichdat.flush();
  culldat.flush();
}

//FIXME err which of these dump methods are actually still functional?


/** \brief Reset everything for next replicate
 * \return void
 */
void Model::resetsim()  {
  t = 0;
  S = N;
  I = 0;  fill(II.begin(),II.end(),0);
  Itot = 0;
  RIP = 0;
  RDC = 0;
  RCP = 0;
  fill(states.begin(),states.end(),0);
  fill(dstate.begin(),dstate.end(),0);
  fill(time_i.begin(),time_i.end(),-1);
  fill(time_r.begin(),time_r.end(),-1);
  fill(time_c.begin(),time_c.end(),-1);
  fill(i_by.begin(),i_by.end(),-1);
  fill(ninfd.begin(),ninfd.end(),0);
  fill(nculd.begin(),nculd.end(),0);
  infd = vector< vector<int> >(5,vector<int>(3,0));
  culd = vector< vector<int> >(5,vector<int>(3,0));
  // TODO careful of memory overheads here! - see massif profile. trade against time...
  /*for (int itmp=0;itmp<tblty.size();++itmp)  {
    fill(tblty[itmp].begin(),tblty[itmp].end(),0.0);
  }*/
  tblty.clear();
  tblty.resize(N,vector<double>(0,0.0));
  //tblty = vector< vector<double> >(N,vector<double>(twmx,0.0));
  //tblty.clear();
  //tblty.resize(N,vector<double>(twmx,0.0));
  //tblty.resize(N); for(int i=0;i<N;++i)  {tblty[i].reserve(20);} // FIXME SPACE!

  fill(sus.begin(),sus.end(),0.0);
  grid.reset();
  for (int ii=0;ii<farms.enreg;++ii)  {
    fill(ertot[ii].begin(),ertot[ii].end(),0.0);
  }
}


/** \brief Calc and store today's error (daily cumulative totals of IPs and culled farms/animals)
 *  Updates the total regional errors based on dstate - that is the changes that have occurred today
 * \return void
 */
void Model::error_calc()  {
  // TODO Eigen-ise this shit.
  if (G_CONST::err_d_c)  {  // if daily
    fill(ninfd.begin(),ninfd.end(),0);
    fill(nculd.begin(),nculd.end(),0);
    for (int reg=0;reg<farms.enreg;++reg)  {
      fill(infd[reg].begin(),infd[reg].end(),0);
      fill(culd[reg].begin(),culd[reg].end(),0);
    }
  }
  for (int i=0;i<N;++i)  {
    switch (dstate[i])  {
      case 0:   // Nothing changes (most likely)
      break;

      case 1:   // Infection
        ++ninfd[farms.eregion[i]];                    // Infected Farm
        infd[farms.eregion[i]][0] += farms.N[i][0];   // Infected cattle
        infd[farms.eregion[i]][1] += farms.N[i][1];   // Infected sheep
      break;

      case 2:   // IP reported
      break;

      case 3:   // DC/CP notified. Only want actuals here.
      break;

      case 4:   // Occult?
        // TODO Handle occult infections?
      break;

      case 5:   // Culling IP
      break;

      case 6:   // Culling DC/CP - are required here for total number of culls
        ++nculd[farms.eregion[i]];
        culd[farms.eregion[i]][0] += farms.N[i][0];
        culd[farms.eregion[i]][1] += farms.N[i][1];
      break;
    }
  }
  //cout << t << "\t";
  for (int reg=0;reg<farms.enreg;++reg)  {
    // Today's differences in daily/cumulative numbers of reported cases and culled premises
    // On regional basis and normalised by metric value on day tmax
    double dninfd = (ninfd[reg]  -farms.farmsi[reg][t])*norm_fi[reg];
    double dinfdc = (infd[reg][0]-farms.cowssi[reg][t])*norm_ci[reg];
    double dinfds = (infd[reg][1]-farms.sheepi[reg][t])*norm_si[reg];
    double dnculd = (nculd[reg]  -farms.farmsc[reg][t])*norm_fc[reg];
    double dculdc = (culd[reg][0]-farms.cowssc[reg][t])*norm_cc[reg];
    double dculds = (culd[reg][1]-farms.sheepc[reg][t])*norm_sc[reg];
    //cout << ninfd[reg]-farms.farmsi[reg][t] << " ";// << ninfd[reg] << " ";
    ertot[reg][0] = ertot[reg][0] + dninfd*dninfd;
    ertot[reg][1] = ertot[reg][1] + dinfdc*dinfdc;
    ertot[reg][2] = ertot[reg][2] + dinfds*dinfds;
    ertot[reg][3] = ertot[reg][3] + dnculd*dnculd;
    ertot[reg][4] = ertot[reg][4] + dculdc*dculdc;
    ertot[reg][5] = ertot[reg][5] + dculds*dculds;
   //ertot[reg]    = ertot[reg] +(dnculd*dnculd) + dculdc*dculdc + dculds*dculds);
  }
}


// WIP vacc around target farm i
/** \brief WIP - not fully designed or implemented...?
 * Checks for farms inside predefined annulus range,
 * \param i int farm id
 * \return void
 */
void Model::vaccinate(int i)  {
  int g = grid.g_id[i];
  auto hh_it = grid.potgrids[g].begin();
  // Check farms in nearby grids
  while (hh_it!=grid.potgrids[g].end())  {
    int h = *hh_it;
    for (auto j_it=grid.f_id[h].begin();j_it!=grid.f_id[h].end();++j_it)  {
      double d2 = dist2(i,*j_it);
      if ((d2>params.vaccrange[0])&&(d2<params.vaccrange[1]))  { // annulus
        dstate[*j_it] = 99; // to be vaccinated...
      }
    }
  }
}


/** \brief Generates a starting sample for whatever. Rerun from Chain/Smc until ready
 * \param ready int&
 * \return Eigen::VectorXd
 */
Eigen::VectorXd Model::init_samp(int& ready)  {
  // TODO "ready" was holdover from old code...
  params.setrand(r);
  ready = 1;
  return(params.par_vec);
}


// TODO sort out run/runsim names to work as intended with Chain/Smc
/** \brief Selfexplanatory interface with Smc
 * \param v Eigen::VectorXd
 * \return void
 */
void Model::run(Eigen::VectorXd v)  {
  parse(v);
  runsim();
}

// TODO sort out errcalc/error_calc names to work as intended with Chain/Smc
/** \brief Selfexplanatory interface with Smc
 * \return Eigen::VectorXd
 */
Eigen::VectorXd Model::errcalc()  {
  Eigen::VectorXd errs(farms.enreg);
  for (int reg=0;reg<farms.enreg;++reg)  {
    errs(reg) = 0.0;
    for (int e=0;e<ertot[reg].size();++e)  {
      errs(reg) += ertot[reg][e];
    }
  }
  // TODO are you actually fucking stupid?
  return(errs);
}


/** \brief Selfexplanatory interface with Smc
 * \param v const Eigen::VectorXd&
 * \return void
 */
void Model::parse(const Eigen::VectorXd& v)  {
  params.parse(v);
  if (params.rwfm&&(params.pwfm==0))  { // Running wfm from trn exp, not fitting here
    params.wfm_samp(r);
  }
}





