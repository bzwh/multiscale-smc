#include <iostream>
#include "Settings.hpp"

namespace G_CONST  {
  const int rej_amh = 0;   // 0-Rejection sampler

  const int err_reg = 5;   // Metric on 5 regions
  const int err_d_c = 1;   // Fit to 1:daily or 0:cumulative errors

  const int fit_reg = 1;   // Single parameter set
  const int fit_plw = 0;   // Linear scaling for Sus & Trans
  const int fit_spc = 3;   // Cows, pigs and sheep. Though pigs are quietly ignored... FIXME!
  const int fit_ker = 1;   // Fit the 2 kernel parameters
  const int fit_det = 0;   // Fit gamma distribution to detection delay. 0:constant
  const int fit_dcs = 0;   // Time-invariant or some functional form or not fitting at all?
  const int fit_dcf = 0;   // Fitting f
  const int multisc = 0;   // Within-farm model or step function

  const int simtmax = 242;
  const int withtmx = 100;
  const int seedday = 23;  // Condition on 23rd Feb

  const int dump_epi = 0;
}

/*
cout << G_CONST::err_reg << " " << G_CONST::err_d_c << " " << G_CONST::fit_reg << " "
       << G_CONST::fit_plw << " " << G_CONST::fit_spc << " " << G_CONST::fit_ker << " "
       << G_CONST::fit_det << " " << G_CONST::multisc << " " << G_CONST::simtmax << " "
       << G_CONST::dump_epi << " " << G_CONST::seedday << endl;
*/

