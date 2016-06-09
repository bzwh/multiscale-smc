#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include "Farm.hpp"


/** \brief Gridding algorithm, want to test entry to grids before considering individuals within
 */
class Grid  {
  public:
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double xrange;
    double yrange;
    double sz;      //!< grid resolution (10km)
    int xgrids;     //!<
    int ygrids;     //!<
    int ngrids;     //!< total number of grids

    std::vector<int> num;                     //!< Number of farms in each grid
    std::vector< std::vector<int> > f_id;     //!< IDs of farms in grid [i]
    std::vector<int> g_id;                    //!< Grid number of farm [i] = gy*xgrids+xg
    std::vector<int> gx;                      //!< x ordinate of grid for farm [i]
    std::vector<int> gy;                      //!< y ordinate of grid for farm [i]
    std::vector< std::vector<int> > potgrids; //!< the grids that g can infect, ie has farms within 15km

    std::vector<int> left;  //!< Number of susceptible farms left in grid during simulation
    std::vector< std::vector<double> > mrate; //!< maxsus[h]*K(d_g,h)

    int whichgrid(int,double,double);         // Identify gridnumber given (x,y)
    double gdist(int,int);                    // Minimum distance between two grids
    void setup(Farms&);
    void reset();
    Grid();

  private:

};


/*
class GridCell  {
  private:

  public:
    int id;  //position in vector
    int num; //needs to be able to be read easily?
    int left;//needs to be reset to num easily?
    std::vector<int> f_id; // ids of farms in cell

    std::vector<int> potgrids; // entries correspond with mrate
    std::vector<double> mrate; // maxsus*gdist to cell's potgrids

};

class GridNew  {
  private:
  public:
    std::vector<GridCell> cells;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double xrange;
    double yrange;
    double sz;
    int xgrids;
    int ygrids;
    int ngrids;
};
*/

#endif
/*
Grid::reset_left()  {
  for (int g=0;g<ngrids;++g)  {
    cells[g].left=cells[g].num;
  }
} // much slower than left=num?
*/
