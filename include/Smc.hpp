#ifndef SMC_HPP_INCLUDED
#define SMC_HPP_INCLUDED

#include <fstream>
#include <Eigen/Dense>
#include "Model.hpp"
#include "Farm.hpp"



class Smc{
public:
  Smc(int,int,int,Farms&);
  void run();
  void gen(int,int);
  Eigen::VectorXd perturbation(int);
  void iterate();
  void runmodel(int,int);
  int test(int,int);
  void update();
  double pert_dens(int,int);
  void write(std::ofstream&,std::ofstream&,std::ofstream&,std::ofstream&);
  void rngfree();

  std::vector<gsl_rng*> r;                    /// RNGs
  std::vector<Model> models;                  /// Model, duh.
  Eigen::MatrixXd particles_old;              /// Vector of particles/parameter vectors
  Eigen::MatrixXd particles_new;              /// Particle proposal
  Eigen::VectorXd priors;                     /// Prior probs of accepted particles
  Eigen::VectorXd weights_old;                /// Old weights for each particle
  Eigen::VectorXd weights_new;                /// Temp storage of new weights, for normalising
  Eigen::VectorXd tolerance;                  /// For each metric, to be shrunk with each step
  Eigen::MatrixXd epsilons;                   /// Calcd error for each particle's metrics
  Eigen::VectorXd means;                      /// To iteratively calc
  Eigen::MatrixXd sigma;                      /// Calc variance, to update kernel
  Eigen::MatrixXd sginv;                      /// Store the inverse too, for pdf
  double ess;                                 /// Effective Sample Size, duh.
  int step;                                   /// How many filtering steps done
  int N;                                      /// Number of particles to accept at each step
  int npar;                                   /// Number of parameters.
  int nround;
  int nmet;
  int nthreads;

private:

};

#endif // SMC_HPP_INCLUDED

