#ifndef FARM_H_INCLUDED
#define FARM_H_INCLUDED

#include <vector>
#include "Settings.hpp"

/** \brief Container for all farm data, all in separate vectors rather than individual Farm classes
 *  Reads demography data, CPs lists, conditions on real outbreak, daily regional fit metrics to compare to.
 *  Nothing within this class changes.
 */
class Farms  {
  public:
    int ncnt;                               //!< Number of farms
    int nreg;                               //!< Number of regions
    std::vector< std::vector<int> > N;      //!<  numbers of cows, pigs, sheep
    std::vector< std::vector<double> > N1;  //!<  inverse, for use in within-farm model
    std::vector<double> x;                  //!< easting  (km)
    std::vector<double> y;                  //!< northing (km)
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    std::vector<int> region;                //!< Cumbria, Devon, Wales, Scotland, Rest of England
    int enreg;// for use in error calc, not parameter and model space
    std::vector<int> eregion;
    std::vector< std::vector<int> > cps;    //!< Contiguous premises
    std::vector<int> seedfarms;             //!< Farms to be seeded - condition on movement ban
    std::vector<int> seedtimes;             //!< Corresponding time of infection
    std::vector<int> seedculls;             //!< Farms that are culled before start date
    std::vector<int> seedtimec;             //!< Corresponding times of culls

    // REAL numbers of infected and culled farms by region, each day. farmsi[day][reg]
    // FIXME hardcoded "cow" and "sheep". Poor piggies =[
    std::vector< std::vector<int> > farmsi;
    std::vector< std::vector<int> > farmsc;
    std::vector< std::vector<int> > cowssi;
    std::vector< std::vector<int> > cowssc;
    std::vector< std::vector<int> > sheepi;
    std::vector< std::vector<int> > sheepc;

    Farms();              //!< Init, just resizes all these vectors for ncnt=188496
    void bigload();       //!< Just calls the rest of these. what's the point? i don't really know
    void loadfarms();     //!< Read demography data from file
    void loadcps();       //!< Read farms' contiguous premises from file
    void loadseeds();     //!< Read first t days of epidemic to seed simulation
    void loaderrdat();    //!< Read fit metric values. Cumulative or Daily?
    void dumpfarms();

  private:
    std::string datapath;
    std::string cpdat;
    std::string fmdat;
    std::string whdat;
    std::string cldat;
    std::string erdat;
    std::string cowdat;
    std::string shpdat;
};

#endif


