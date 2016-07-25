#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <Eigen/Dense>
#include "Params.hpp"
#include "Settings.hpp"

using namespace std;
using namespace Eigen;

/// NOTE: Current (bad) model selection: priors min/max, nreg, plaw, pker

/// NOTE: change parameters -> priors, sigmas, itopram, pnum, mu/cov updates, #check, MESSY
/// TODO: externalise priors; read from file, store with outputs
/// TODO: just keep cov and mu as MatrixXd VectorXd, instead of jumping between
//  req some map pnum i -> sb/sp/tb/tp/ker/etc
// But where to put this? mu&cov in Chain as dep on #iterations.
// don't want to be copying around with pold=pnew?

// Checks file to spec model. Call setup() for storage. Init priors. Call setacc() starting values
Params::Params()  {
  ifstream pdat;
  pdat.open("INPUT"); // better to just toggle this and open in pread?
  if (pdat.is_open())  {
    pread(pdat);
    pdat.close();
  }
  else{
    // Determines model and parameters to fit
    nreg = G_CONST::fit_reg;         // number of regions to split parameterisation in to
    plaw = G_CONST::fit_plw;         // using power law scaling or not?
    pker = G_CONST::fit_ker ? 2 : 0; // Fitting the 2 kernel parameters?
    nspc = G_CONST::fit_spc;         // cows/pigs/sheep?
    pdet = G_CONST::fit_det ? 2 : 0; // Fitting ~gamma(u,k)
    pdcs = G_CONST::fit_dcs;         // 2-dist, 1-const, 0-nope
    pdcf = G_CONST::fit_dcf;         // Fitting DC f accuracy
    with = G_CONST::multisc;         // 1: Within farm SEmInR model  0: Step function

    pper = (plaw) ? 7 : 3;   // number of regionalised parameters (sus and trans, no powerlaws)
    pnum = nreg*pper+pker+pdet+pdcs+pdcf;
    storage_setup(); // storage

    // Initialising priors here??? - put in list??
    // Susceptibility and Transmissibility priors
    for (int reg=0;reg<nreg;++reg)  {
      fill(sus_min[reg].begin(),sus_min[reg].end(),0.0);
      fill(sus_max[reg].begin(),sus_max[reg].end(),30.0);
      fill(trn_min[reg].begin(),trn_min[reg].end(),0.0);
      trn_max[reg][0] = 1.0e-4;//1.0e-3;
      trn_max[reg][1] = 1.0; // Who cares right now?
      trn_max[reg][2] = 1.0e-5;//1.0e-4;
      //fill(trn_max[reg].begin(),trn_max[reg].end(),1.0e-4);
    }
    // Kernel priors
    if (pker)  {
      ker_min[0] = 0.0; ker_max[0] = 2.0; // K(r) = k0/2
      ker_min[1] = 1.0; ker_max[1] = 10.0;
    }
    // Power law exponent piors
    if (plaw)  {
      sp_min = 0.0;
      sp_max = 1.0;
      tp_min = 0.0;
      tp_max = 1.0;
    }
    // Delays - detection/reporting/cull
    // TODO Regionalised detection and cull delays? probably not...
    if (pdet)  {
      detmin = 1.0;
      detmax = 12.5;
    }
    delaydcp = 2;
    delayipc = 1;

    // DC cull parameters
    if (pdcs==2)  {
      dcFmin[1] = 1.0/150.0;         // So F peaks t=(50,150)
      dcFmax[1] = 1.0/50.0;
      dcFmin[0] = 100.0*2.7183*dcFmin[1]; // So that F can peak at 10?
      dcFmax[0] = 100.0*2.7183*dcFmax[1];
    }
    else if (pdcs==1)  {// fitting constant F. (Mike had F==6)
      dcFmin[0] = 5.0;
      dcFmax[0] = 10.0;
    }
    if (pdcf)  {
      dcfmin = 0.75;//0.5;
      dcfmax = 1.0;//1.0;
    }
    f = 0.85;
    F[0] = 6.0;
    // Kernel parameters
    ker[0] = 1.0;
    ker[1] = 3.0; // Boender kernel, Marleen's estimates

    // Vaccination strategy - wip
    vaccrange.resize(2,0.0);
    vaccrange[1] = 5.0;
  }
}


/** \brief Loads gamma parameters from SG cattle and BH sheep transmission anaylsis.
 * Don't need to call if not running within-farm model
 * \return void
 */
void Params::load_transmission()  {
  ifstream cowdat;
  #ifdef _WIN32_WINNT
    cowdat.open("D:/Ben/data/uk/GammaParameters_COW.txt");
  #endif // _WIN32_WINNT
  #ifdef __linux__
    cowdat.open("../../data/uk/GammaParameters_COW.txt");
  #endif // __linux__

  if(cowdat.is_open())  { // m muE n muI B
    int num = std::count(std::istreambuf_iterator<char>(cowdat),
                      std::istreambuf_iterator<char>(),'\n');
    cowdat.clear();             // reset eof and..
    cowdat.seekg(0,ios::beg);   //   ... jump back to beginning of file

    m_cows.resize(num,0);
    sigma_cows.resize(num,0.0);
    n_cows.resize(num,0);
    gamma_cows.resize(num,0.0);
    beta_cows.resize(num,0.0);
    double mm,nn,ss,gg,bb;
    for (int i=0;i<num;++i)  {
      cowdat >> mm >> ss >> nn >> gg >> bb; // NB ss and gg here are mean periods - not transition rates
      m_cows[i] = mm;
      sigma_cows[i] = mm/ss;  // COW(:,[2 4]) = mean latent and infectious periods. inverse for rate
      n_cows[i] = nn;
      gamma_cows[i] = nn/gg;
      beta_cows[i] = bb;
    }
  }
  else  {
    cout << "where's my goddamn cow data?!?!?!?!?!?!?!?!?" << endl;
    exit(-1);
  }
  cowdat.close();

  ifstream shpdat;
  #ifdef _WIN32_WINNT
    shpdat.open("D:/Ben/data/uk/GammaParameters_SHEEP.txt");
  #endif // _WIN32_WINNT
  #ifdef __linux__
    shpdat.open("../../data/uk/GammaParameters_SHEEP.txt");
  #endif // __linux__
  if(shpdat.is_open())  { // m 1/muE n 1/muI B
    int num = std::count(std::istreambuf_iterator<char>(shpdat),
                      std::istreambuf_iterator<char>(),'\n');
    shpdat.clear();             // reset eof and..
    shpdat.seekg(0,ios::beg);   //   ... jump back to beginning of file

    m_shps.resize(num,0);
    sigma_shps.resize(num,0.0);
    n_shps.resize(num,0);
    gamma_shps.resize(num,0.0);
    beta_shps.resize(num,0.0);
    double mm,nn,ss,gg,bb;
    for (int i=0;i<num;++i)  {
      shpdat >> mm >> ss >> nn >> gg >> bb;
      m_shps[i] = mm;
      sigma_shps[i] = mm/ss;
      n_shps[i] = nn;
      gamma_shps[i] = nn/gg;
      beta_shps[i] = bb;
//      cout << m_shps[i] << "\t" << sigma_shps[i] << "\t" << n_shps[i] << "\t" << gamma_shps[i] << "\t" << beta_shps[i] << "\n";
    }
  }
  else  {
    cout << "where's my goddamn sheep data?!?!?!?!?!?!?!?!?" << endl;
    exit(-1);
  }
  shpdat.close();
}


/** \brief Dead. Can just look where pwrite has been pointed (probably)
 * \return void
 */
void Params::pscreendump()  {
  cout << endl;
  cout << "sb \t\t\t tb \t\t\t sp \t\t\t tp \t\t\t \n";
  for (int i=0;i<nreg;++i)  {
    // each region
    for (int j=0;j<nspc;++j)  {
      // each species
      cout << sb[i][j] << "\t";
    }
    for (int j=0;j<nspc;++j)  {
      // each species
      cout << tb[i][j] << "\t";
    }
    for (int j=0;j<nspc;++j)  {
      // each species
      cout << sp[i][j] << "\t";
    }
    for (int j=0;j<nspc;++j)  {
      // each species
      cout << tp[i][j] << "\t";
    }
    cout << "\n";
  }
}



/** \brief Write parameter set to file. Dumping each region's parameters in block of 7, end with ker
 * \param pdat  opened ofstream to write to
 * \return void
 */
void Params::pwrite(std::ofstream& pdat)  {
  pdat << par_vec << endl;
}



/** \brief With model specified, initialise storage space for parameters and priors
 * Doesn't care whether read from file or hardcoded.
 */
void Params::storage_setup()  {
  // Setup: Within-farm SEmInR parameters
  /*m.resize(nspc,0); sigma.resize(nspc,0.0);
  n.resize(nspc,0); gamma.resize(nspc,0.0);
  beta.resize(nspc,vector<double>(nspc,0.0));*/
  if (with)  {
    load_transmission();
  }
  // Setup: Sus/Trans regional parameters
  sb.resize(nreg,vector<double>(nspc,0.0));
  tb.resize(nreg,vector<double>(nspc,0.0));
  // Setup: Powerlaws
  if (plaw)  { // default to [0,1]
    sp.resize(nreg,vector<double>(nspc,0.5));
    tp.resize(nreg,vector<double>(nspc,0.5));
  }
  ddet.resize(pdet);
  switch (pdet)  {
    case 1:
      ddet[0] = 4;
      break;
    case 2:
      ddet[0] = 4.0;
      ddet[1] = 1.0;
  }
  // Setup: Kernel parameters
  ker.resize(2,0.0);
  ker_min.resize(pker,0.0);
  ker_max.resize(pker,0.0); // cheap, no real need to check if these are actually going to be used
  // Vaccination WIP
  vaccrange.resize(2,0.0);
  // Priors: Sus/Trans
  sus_min.resize(nreg,vector<double>(nspc,0.0));
  sus_max.resize(nreg,vector<double>(nspc,0.0));
  trn_min.resize(nreg,vector<double>(nspc,0.0));
  trn_max.resize(nreg,vector<double>(nspc,0.0));
  if (pdcs)  {
    F.resize(pdcs,0.0);
  }
  else  {
    F.resize(1,0.0);
  }
  dcFmin.resize(pdcs,0.0);
  dcFmax.resize(pdcs,0.0);
  par_vec.resize(pnum);
}

/// TODO init on file. req's model selection and 2 sets of parameters. (means and cov elsewhere)
/// Basically, externalise the initialisation. easier to read what's going on in different runs
/// Worry about continuing a chain later...
/** \brief WIP:Init from values in a file, determines the model to use..?
 * \param pdat ifstream& where model is spec'd, initial values and priors given
 * \return void
 */
void Params::pread(std::ifstream& pdat)  {
cout << "HELLO!" << endl;
  if (pdat.is_open())  { // only gets here if open!
    // Init model
    pdat >> nreg >> plaw >> pker >> nspc;
    pper = (plaw) ? 7:3; // 4*nspc-1 : 2*nspc-1  or will scale to speciifc region's S_sheep?
    pnum = nreg*pper + pker;  // Total number of parameters being fitted. Delays? Vacc?
    // Get storage space
    storage_setup();
    // Read in sus and trans priors (bounded uniform)
    for (int reg=0;reg<nreg;++reg)  {
      for (int spc=0;spc<nspc;++spc)  {
        pdat >> sus_min[reg][spc] >> sus_max[reg][spc] >> trn_min[reg][spc] >> trn_max[reg][spc];
      }
    }
    // Priors: powerlaws
    if (plaw)  {
      pdat >> sp_min >> sp_max >> tp_min >> tp_max;
    }
    // Priors: Kernel
    if (pker)  {
      for (int kk=0;kk<pker;++kk)  {
        pdat >> ker_min[kk] >> ker_max[kk];
      }
    }
    setacc();
  } // FIXME some sort of catch for end of file? ie, spec'd but not init'd?
  else  {
    cout << "No parameter input... Default or quit?" << endl;
    setacc();
  }
}


/// TODO for init and rejection
// get some range to test on and make sure it gets the chain started. depends on the model...
// probs want overdispersed starting point.
void Params::setrand(gsl_rng* r)  {
  for (int reg=0;reg<nreg;++reg)  {
    sb[reg][0] = sus_min[reg][0] + (sus_max[reg][0]-sus_min[reg][0])*gsl_rng_uniform(r);
    sb[reg][1] = 0.00001;
    sb[reg][2] = 1.0; // Sus_sheep = 1.0;
    tb[reg][0] = trn_min[reg][0] + (trn_max[reg][0]-trn_min[reg][0])*gsl_rng_uniform(r);
    tb[reg][1] = 0.00001;
    tb[reg][2] = trn_min[reg][2] + (trn_max[reg][2]-trn_min[reg][2])*gsl_rng_uniform(r);
  }
  if (pker)  {
    ker[0] = ker_min[0] + (ker_max[0]-ker_min[0])*gsl_rng_uniform(r);
    ker[1] = ker_min[1] + (ker_max[1]-ker_min[1])*gsl_rng_uniform(r);
  }
  else  { // Not nec... (depending on original init!)
    ker[0] = 0.12;
    ker[1] = 3.0;
  }
  if (pdet==1)  {
    //delaydet = floor(detmin+(detmax-detmin+1)*gsl_rng_uniform(r));
    ddet[0] = detmin+(detmax-detmin+1)*gsl_rng_uniform(r);
  }
  else if (pdet==2)  {
    ddet[0] = detmin+(detmax-detmin+1)*gsl_rng_uniform(r);
    ddet[1] = 5.0*gsl_rng_uniform(r);
  }
  else  {
    delaydet = 4.0;
  }
  delaydcp = 2;
  delayipc = 1;
  for (int i=0;i<pdcs;++i)  {
    F[i] = dcFmin[i]+(dcFmax[i]-dcFmin[i])*gsl_rng_uniform(r);
  }
  if (pdcf)  {
    f = dcfmin+(dcfmax-dcfmin)*gsl_rng_uniform(r);
  }
  else  {
    f = 0.85;
    F[0] = 6.0;
  }
  // TODO clean this up - populating par Eigen::VectorXd
  // FIXME HORRIBLE NESTED DUPLICATED CODE SHIT - this looks nicer?
  int iblah = 0;
  par_vec[iblah++] = sb[0][0];
  par_vec[iblah++] = tb[0][0];
  par_vec[iblah++] = tb[0][2];
  for(int i=0;i<pker;++i)  {
    par_vec[iblah++] = ker[i];
  }
  for(int i=0;i<pdet;++i)  {
    par_vec[iblah++] = ddet[i];
  }
  for(int i=0;i<pdcs;++i)  {
    par_vec[iblah++] = F[i];
  }
  if (pdcf)  {
    par_vec[iblah++] = f;
  }
  //(iblah==pnum) ? cout<<"!"<<endl : cout<<"WTF"<<endl; // checking right number of parameters...

  // TODO while(boundscheck)???
}


/// NOTE ANOTHER CHANGE HERE FOR MODEL - trans smaller if no plaw
/** \brief Setting to Mike's estimates from '08 Accuracy of models paper (and suppressing pigs)
 */
void Params::setacc()  {
// Parameters fitted in Accuracy08. sort of.
  if (nreg==1)  {
    sb[0][0]=7.0;    sb[0][1]=0.01;    sb[0][2]=1.0;
    tb[0][0]=1.3e-6; tb[0][1]=0.01e-7; tb[0][2]=3.0e-7;
    if (plaw)  {
      tb[0][0] = 1.0e-5;
      tb[0][2] = 1.0e-5;  // Takes a long time to get away from the much smaller plaw-less estimates
      for (int reg=0;reg<nreg;++reg)  {
        // NOTE tried starting close to powerlaw-less estimate, did not work.
        fill(sp[reg].begin(),sp[reg].end(),0.75);
        fill(tp[reg].begin(),tp[reg].end(),0.75);
      }
      /*sb[0][0]=7.0;    sb[0][1]=0.01;    sb[0][2]=1.0;  // These will end up much larger...
      tb[0][0]=1.3e-5; tb[0][1]=0.01e-7; tb[0][2]=3.0e-5;
      sp[0][0]=0.5;    sp[0][1]=0.5;     sp[0][2]=0.5;  // dont really need these...?
      tp[0][0]=0.5;    tp[0][1]=0.5;     tp[0][2]=0.5;*/
    }
  }
  else  {
    sb[0][0]= 5.7;    sb[0][1]=0.01;    sb[0][2]=1.0;
    sb[1][0]= 4.9;    sb[1][1]=0.01;    sb[1][2]=1.0;
    sb[2][0]= 0.7;    sb[2][1]=0.01;    sb[2][2]=1.0;
    sb[3][0]=10.2;    sb[3][1]=0.01;    sb[3][2]=1.0;
    sb[4][0]= 2.3;    sb[4][1]=0.01;    sb[4][2]=1.0;

    tb[0][0]= 8.2e-7; tb[0][1]=0.01e-7; tb[0][2]= 8.3e-7;
    tb[1][0]= 5.8e-4; tb[1][1]=0.01e-4; tb[1][2]=11.0e-4;
    tb[2][0]=30.1e-4; tb[2][1]=0.01e-4; tb[2][2]=36.3e-4;
    tb[3][0]=23.2e-4; tb[3][1]=0.01e-4; tb[3][2]=28.2e-4;
    tb[4][0]= 8.2e-4; tb[4][1]=0.01e-4; tb[4][2]=23.2e-4;

    if (plaw)  {
      sp[0][0]=0.41;    sp[0][1]=0.1;     sp[0][2]=0.20;
      sp[1][0]=0.37;    sp[1][1]=0.1;     sp[1][2]=0.40;
      sp[2][0]=0.31;    sp[2][1]=0.1;     sp[2][2]=0.43;
      sp[3][0]=0.23;    sp[3][1]=0.1;     sp[3][2]=0.33;
      sp[4][0]=0.42;    sp[4][1]=0.1;     sp[4][2]=0.30;

      tp[0][0]=0.42;    tp[0][1]=0.1;     tp[0][2]=0.49;
      tp[1][0]=0.37;    tp[1][1]=0.1;     tp[1][2]=0.42;
      tp[2][0]=0.25;    tp[2][1]=0.1;     tp[2][2]=0.22;
      tp[3][0]=0.20;    tp[3][1]=0.1;     tp[3][2]=0.40;
      tp[4][0]=0.44;    tp[4][1]=0.1;     tp[4][2]=0.37;
    }
  }
  //ker[0]=0.12;
  ker[0] = 1.0;
  ker[1] = 3.0;
  //ker[1] = 1.0/(ker[1]*ker[1]);
  //ker[2] = ker[2]/2.0;
  if (pdet==1)  {
    //delaydet = floor(detmin+(detmax-detmin+1)*gsl_rng_uniform(r));
    ddet[0] = 4.0;
  }
  else if (pdet==2)  {
    ddet[0] = 4.0;
    ddet[1] = 0.1;
  }
  else  {
    delaydet = 4;
  }
  if (pdcs==2)  {
    F[0] = 10.0*65.0*exp(-1.0);
    F[1] = 1.0/65.0;
  }
  else  {
    F[0] = 6;
  }

  delaydcp = 2;
  delayipc = 1;
  f = 0.85;
  // TODO clean this up - populating par Eigen::VectorXd
  // FIXME HORRIBLE NESTED DUPLICATED CODE SHIT
  int iblah = 0;
  par_vec[iblah++] = sb[0][0];
  par_vec[iblah++] = tb[0][0];
  par_vec[iblah++] = tb[0][2];
  if (pker)  {
    par_vec[iblah++] = ker[0];
    par_vec[iblah++] = ker[1];
  }
  if (pdet)  {
    par_vec[iblah++] = ddet[0];
    par_vec[iblah++] = ddet[1];
  }
  if (pdcs)  {
    par_vec[iblah++] = F[0];
    par_vec[iblah++] = F[1];
  }
  if (pdcf)  {
    par_vec[iblah++] = f;
  }
}




/** \brief Parse the vector of parameters in to named fields
 * \param v const Params&
 * \return void
 */
void Params::parse(const VectorXd& v)  {
  par_vec = v;
  int i = 0;  // iterating through the full vector - check for potential bounds errors
  for (int reg=0;reg<nreg;++reg)  {
    // TODO may well work, but only for 1 region. needs splitting for S_sheep. reg*pper?
    sb[reg][0] = v[i++];
    tb[reg][0] = v[i++];
    tb[reg][2] = v[i++];
    if (plaw)  {
      sp[reg][0] = v[i++];
      sp[reg][2] = v[i++];
      tp[reg][0] = v[i++];
      tp[reg][2] = v[i++];
    }
  }
  for (int k=0;k<pker;++k)  { // if not fitting kernel, pker==0
    ker[k] = v[i++];
  }
  for (int k=0;k<pdet;++k)  {// Delay from farm infectiousness/clinical signs to report
    ddet[k] = v[i++];
  }
  for (int k=0;k<pdcs;++k)  {
    F[k] = v[i++]; // Measure of DC:IP ratio
  }
  if (pdcf)  {
    f = v[i++]; // Accuracy of identifying DCs
  }
  // TODO allow for delays and control parameters
  delaydcp = 2;   // Delay from notification to removal
  delayipc = 1;   // Delay from report to removal

}


// Should these checks be going in Chain? Not nec for rejection sampling
/** \brief Check parameter values are within bounds. Uniform priors.
 * \return int
 */
int Params::prior_check()  {
  // Regional species-specific susceptibilities
  for (int reg=0;reg<nreg;++reg)  {
    for (int spc=0;spc<nspc;++spc)  {
      if ((sb[reg][spc]<sus_min[reg][spc])||(sb[reg][spc]>sus_max[reg][spc]))  {
        return(-1);
      }
    }
  }
  // Regional species-specific transmissibilities
  for (int reg=0;reg<nreg;++reg)  {
    for (int spc=0;spc<nspc;++spc)  {
      if ((tb[reg][spc]<trn_min[reg][spc])||(tb[reg][spc]>trn_max[reg][spc]))  {
        return(-2);
      }
    }
  }
  // Fitting powerlaws?
  if (plaw)  {  // Just a boolean state
    // Susceptibility exponent
    for (int reg=0;reg<nreg;++reg)  {
      for (int spc=0;spc<nspc;++spc)  {
        if ((sp[reg][spc]<sp_min)||(sp[reg][spc]>sp_max))  {
          return(-3);
        }
      }
    }
  // Transmissibility exponent
    for (int reg=0;reg<nreg;++reg)  {
      for (int spc=0;spc<nspc;++spc)  {
        if ((tp[reg][spc]<tp_min)||(tp[reg][spc]>tp_max))  {
          return(-4);
        }
      }
    }
  }
  // Fitting kernel parameters?
  for (int i=0;i<pker;++i)  { // Is actual number of parameters being fitted
    if ((ker[i]<ker_min[i])||(ker[i]>ker_max[i]))  {
      return(-5);
    }
  }
  // Fitting DC cull F - functional form or constant?
  for (int i=0;i<pdcs;++i)  {
    if ((F[i]<dcFmin[i])||(F[i]>dcFmax[i]))  {
      return(-6);
    }
  }
  // Fitting DC cull f
  if (pdcf)  {
    if ((f<dcfmin)||(f>dcfmax))  {
      return(-7);
    }
  }
  return(0);  // All within whatever bounds you've decided....
}




