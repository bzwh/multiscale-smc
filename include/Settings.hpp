#ifndef SETTINGS_HPP
#define SETTINGS_HPP

namespace G_CONST  {
  extern const int rej_amh;     //!< TODO Rejection sampler or adaptive metropolis

  extern const int err_reg;     //!< Daily or cumulative
  extern const int err_d_c;     //!< Regionalised fit metric - for Farms and Epi

  extern const int fit_reg;     //!< Regionalised parameters
  extern const int fit_plw;     //!< Power law or linear scaling on animal numbers
  extern const int fit_spc;     //!< Number of species
  extern const int fit_ker;     //!< Fitting kernel parameters?
  extern const int fit_det;     //!< Fitting detection delay parameters?
  extern const int fit_dcs;     //!< Fitting culls F parameters (distribution?)
  extern const int fit_dcf;     //!< Fitting cull f parameter - tracing accuracy
  extern const int multisc;     //!< Using multi-scale model or step-function

  extern const int cp_cull;     //!< Contiguous premises culled? (UK vs JPN)

  extern const int simtmax;     //!< How many days to run the simulation
  extern const int withtmx;     //!< Max length of within-farm epidemic run
  extern const int seedday;     //!< When to seed epidemic

  extern const int dump_epi;     //!< Dump simulations to WHICH & CULL files;

}
#endif // SETTINGS_HPP
