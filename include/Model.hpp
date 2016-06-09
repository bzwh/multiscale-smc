#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <gsl/gsl_rng.h>
#include <Eigen/Dense>

#include "Farm.hpp"
#include "Grid.hpp"
#include "Params.hpp"
#include "Settings.hpp"

/** \brief Where the epidemic is realised. loads up all farm data and sets up gridding.
  Give it a parameter set, run the simulation, output the fit metric to update chain.
 */
class Model  {
public:
  // TODO
  Eigen::VectorXd init_samp(int&);
  void run(Eigen::VectorXd);
  Eigen::VectorXd errcalc();



  // Methods
  Model();
  void setup(gsl_rng*,const Farms&);        //!< Just grabs Farm data and stores away the RNG
  void setrdat(int,std::string);  //!< Comprehensive output of simulation runs
  void initrun();         //!< Give parameter set to run simulation
  void parse(const Eigen::VectorXd&);
  void sus_calc();            //!< Calc susceptibility for all farms with this parameter set
  void maxrate_calc();        //!< Calc maximum sus*dker between grids
  void dccullF_calc();        //!< Calc the DC cull intensity at time t
  void seed_run();            //!< Seed the epidemic
  int runsim();               //!< Run the epidemic - return 1 if too big!
  void iterate();             //!< Between-farm transmissions
  void detrepcul();           //!< Detection/Report/Cull transitions
  void update();              //!< Implement transitions
  void cpcull(int);           //!< Pre-calculated CPs to cull
  void dccull(int);           //!< Calculate DCs to cull
  double dist2(int,int);      //!< Distance^2 between two farms
  double dker(double);        //!< Transmission kernel
  double trans(int);          //!< Pull transmissibility of farm at current time
  void vaccinate(int);        //!< Ring vaccination around farm id given
  void error_calc();          //!< Update regional numbers of cases and culls farms/cows/sheep
  void dumpdailyir();         //!< Write individual run data to files, just the daily number of new infections
  void dumpfarmids();         //!< Write all farms' state and infection&cull times
  void dumpoutbreak(std::ofstream&,std::ofstream&);        //!< Write outbreak in same style as input Which&CULL
  void resetsim();            //!< Clear everything ready for next replicate
  void without(int);          //!< Step function transmission profile
  void within(int);           //!< Simulated within-farm epidemic
  std::function<void(int)> within_model;


  // Member data
  gsl_rng* r;
  int N;
  int plaw; // to come from params...
  int pdet; // to come from params...
  Farms farms;
  Grid grid;
  Params params;

  // Stuff for each run
  int tmax;                 //!< Time to end simulation. t=242 is 01/10/2001
  int twmx;                 //!< Max within-farm simulation time
  int t;                    //!< Current simulation time
  //std::ofstream rdat;       //!< Where to dump the output of individual runs
  std::vector<int> states;  //!< 0:Sus 1:Inf  2:Rep  3:Notified -1:R-IP -2:R-DC/CP
  std::vector<int> dstate;  //!< Today's state transitions
  std::vector<int> time_i;  //!< When was i infected
  std::vector<int> time_r;  //!< When was i reported
  std::vector<int> time_c;  //!< When was i culled
  std::vector<int> i_by;    //!< Who infected i
  std::vector<double> sus;  //!< Susceptibility of farm i
  std::vector< std::vector<double> > tblty; //!< Transmissibility of farm i
  double dcF;

  // A few random things tracking size of outbreak. To be dumped?
  int S;                    //!< Today's number of S
  std::vector<int> SS;
  int I;                    //!< Today's number of I
  int Itot;
  std::vector<int> II;
  std::vector<int> C_d;     //!< Today's number of something(?)
  std::vector<int> R_d;     //!< Today's number of reported IPs
  int RIP;                   //!< Total number of IPs culled
  int RDC;                  //!< Total number of DCs culled
  int RCP;                  //!< Total number of CPs culled
  std::vector<int> igrids;  //!< Tracking numbers of infections in each grid(?)

  // Tracking epidemic progress to calc fit metric. Store daily stuff for entire run?
  std::vector<int> ninfd;               //!< Regional number of infected farms
  std::vector<int> nculd;               //!< regional number of culled farms
  std::vector< std::vector<int> > infd; //!< Regional numbers of infected cows/pigs/sheep
  std::vector< std::vector<int> > culd; //!< Regional numbers of culled cows/pigs/sheep
  std::vector< std::vector<double> > ertot;            //!< Regional errors, as def by fit metric

  // Normalising constants for fit metric
  std::vector<double> norm_fi;
  std::vector<double> norm_fc;
  std::vector<double> norm_ci;
  std::vector<double> norm_cc;
  std::vector<double> norm_si;
  std::vector<double> norm_sc;





private:

};

#endif  // MODEL_HPP_INCLUDED



