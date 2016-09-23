#include <iostream>
#include "Settings.hpp"

namespace G_CONST  {
  const int rej_amh = 0;   // 0-Rejection sampler

  const int err_reg = 1;   // Metric on 5 regions
  const int err_day = 1;   // Fit to 1:daily or 0:cumulative errors
  const int err_dcp = 0;   // Metric on culled farms & animals (as well as infecteds)
  const int err_vac = 0;   // Metric on vaccinated farms?

  const int fit_reg = 1;   // How many parameter sets
  const int fit_plw = 0;   // Linear scaling for Sus & Trans
  const int fit_spc = 2;   // Cows, pigs and sheep. Though pigs are quietly ignored... FIXME!
  const int fit_ker = 1;   // Flag - fit kernel parameters
  const int fit_det = 1;   // Fit gamma distribution to detection delay. 0:constant, 1:gamma
  const int fit_dcs = 2;   // Time-invariant or some functional form or not fitting at all? 2:dist,1:const,0:off
  const int fit_dcf = 1;   // Fitting f
  const int fit_msc = 0;   // Fitting WFM parameters (k,mu,beta)
  const int fit_vac = 0;   // Fit vaccination scaling parameters

  const int multisc = 0;   // Within-farm model or step function
  const int cp_cull = 0;   // Contiguous premises to be culled or not? (ie UK vs JPN...)
  const int dc_cull = 0;   // Simulate DC culls or not?
  const int rn_vacc = 0;   // Run ring vaccination

  //const int simtmax = 242; // UK
  const int simtmax = 58;  // JP
  const int withtmx = 20;
  const int seedday = 20;  // t=23 to condition on 23rd Feb movement ban (UK)

  const int dump_epi = 0;
}

/*
cout << G_CONST::err_reg << " " << G_CONST::err_d_c << " " << G_CONST::fit_reg << " "
       << G_CONST::fit_plw << " " << G_CONST::fit_spc << " " << G_CONST::fit_ker << " "
       << G_CONST::fit_det << " " << G_CONST::multisc << " " << G_CONST::simtmax << " "
       << G_CONST::dump_epi << " " << G_CONST::seedday << endl;
*/

