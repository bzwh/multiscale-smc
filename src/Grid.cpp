#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Grid.hpp"
#include "Farm.hpp"

using namespace std;


/** \brief Init empty grid, 10x10 cells for UK.
 */
Grid::Grid()  {
  sz = 10.0;    // 10x10 km grids
  xmin = 0.0;
  xmax = 0.0;
  ymin = 0.0;
  ymax = 0.0;
  xrange = 0.0;
  yrange = 0.0;
  xgrids = 0;
  ygrids = 0;
  ngrids = 0;
  /*xmin = 0.0;   // UK
  xmax = 700.0;
  ymin = 0.0;
  ymax = 1250.0;

  xmin = 0.0;   // JPN
  xmax = 140.0;
  ymin = 0.0;
  ymax = 200.0;
*/
}


/** \brief Set up the grid data; where farms are, how many in each grid and which grids can infect which.
 * \param farms Farms&
 * \return void
 */
void Grid::setup(const Farms& farms)  {
  xmin = 0.0; //floor(farms.xmax/sz)*sz;
  xmax = ceil(farms.xmax/sz)*sz;
  ymin = 0.0; //floor(farms.ymax/sz)*sz;
  ymax = ceil(farms.ymax/sz)*sz;

  xrange = xmax-xmin;
  yrange = ymax-ymin;

  xgrids = ceil(xrange/sz);
  ygrids = ceil(yrange/sz);
  ngrids = xgrids*ygrids;
  //cout << endl << xmax << "\t" << ymax << endl << ngrids << endl<<endl;
  num.resize(ngrids,0);     // total number of farms in grid i
  left.resize(ngrids,0);    // number of susceptible farms in grid
  potgrids.resize(ngrids);  // ids of non-empty grids near enough to be infected
  mrate.resize(ngrids);     // maxsus*dker
  f_id.resize(ngrids);      // IDs of farms in each grid

  g_id.resize(farms.ncnt,-1);   // Grid number of each farm = gy*xgrids+gx
  gx.resize(farms.ncnt,-1);     // x ordinate of grid
  gy.resize(farms.ncnt,-1);     // y ordinate of grid
  for (int i=0;i<farms.ncnt;++i)  {
    int grid_id = whichgrid(i,farms.x[i],farms.y[i]);
    g_id[i] = grid_id;
    f_id[grid_id].push_back(i);
    ++num[grid_id]; // num[g] == f_id[g].size
  }
  left = num;
  // Set up potgrids - guarantee nonempty and close enough
  for (int g=0;g<ngrids;++g)  {
    if (num[g]==0)  {
      continue; // g is empty!
    }
    else  {
      for (int h=0;h<ngrids;++h)  {
        if (num[h]==0)  {
          continue; // h is empty!
        }
        else if (g==h)  {
          continue; // Ignore self- dealt specially in Epi::iterate()
        }
        else if (gdist(g,h)>=0.0)  {
          potgrids[g].push_back(h);
        }
        else  {
          continue;
        }
      }
    }
    mrate[g].resize(potgrids[g].size());
  }
}


/** \brief Identify gridnumber from easting-northing coordinates
 *
 * \param i int - farm id number
 * \param x double - easting(km)
 * \param y double - northing(km)
 * \return int - grid number
 *
 */
int Grid::whichgrid(int i, double x, double y)  {
  double xx = x-xmin;
  double yy = y-ymin;
  gx[i] = floor(xx*xgrids/xrange);
  gy[i] = floor(yy*ygrids/yrange);
  return( gx[i] + gy[i]*xgrids );
}


/** \brief Shortest distance between two grids, ie 0 for nearest 8 neighbours
 *
 * \param g int a grid number
 * \param h int another grid number
 * \return Distance. Or -1 if further than 15km...?
 *
 */
double Grid::gdist(int g,int h)  {
  int xg = g%xgrids;
  int xh = h%xgrids;
  int yg = floor(g/xgrids);
  int yh = floor(h/xgrids);
  double dx = sz*max(0,abs(xg-xh)-1);
  double dy = sz*max(0,abs(yg-yh)-1);
  if ((dx>15.0)||(dy>15.0))  {
    // FIXME 15km global fixed maximum? depend on kernel parameters?
    return(-1.0);
  }
  else  {
    return(dx*dx+dy*dy);
  }
}


/** \brief Reset whatever has changed during last sim ready for next.
 *
 * \return void
 *
 */
void Grid::reset()  {
  left = num; // reset grids' numbers of susceptibles
  for (int g=0;g<ngrids;++g)  {
    // not really nec... should get over-written anyway in Epi::maxrate_calc
    fill(mrate[g].begin(),mrate[g].end(),-1.0);
  }
  //grid.mrate = vector< vector<double> >(grid.ngrids,vector<double>(grid.ngrids,0.0));
  /*mrate.clear();
  mrate.resize(ngrids);
  for (int g=0;g<ngrids;++g)  {
    mrate[g].reserve(10);
  }
  potgrids.clear();
  potgrids.resize(ngrids);
  for (int g=0;g<ngrids;++g)  {
    potgrids[g].reserve(10);
  }*/
}

