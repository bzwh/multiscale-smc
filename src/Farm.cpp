#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Farm.hpp"
#include "Settings.hpp"

using namespace std;



/** \brief Initialise space required for UK demography data (N=188496)
 */
Farms::Farms()  {
  // NOTE INPUT DATA PATH
  #ifdef _WIN32_WINNT
    // Make relative???
    datapath = "D:/Ben/data/uk/";
  #endif // _WIN32_WINNT
  #ifdef __linux__
    datapath = "../../data/uk/";
  #endif // __linux__

  // NOTE INPUT DATA FILENAMES: demography, CPs, infection/cull seeds, fit metric
  cpdat = datapath + "CPs_NEW";
  whdat = datapath + "Which_Farms_BH";
  cldat = datapath + "CULL_Farms_BH";
  if (G_CONST::err_reg==1)  {
    fmdat = datapath + "Farm_Data_BH_onereg";
    if (G_CONST::err_d_c)  {
      erdat = datapath + "errors_daily_onereg.txt";
    }
    else  {
      erdat = datapath + "errors_cumulative_onereg.txt";
    }
  }
  else  { // Cumbria, Devon, Wales, Scotland, RoUK
    fmdat = datapath + "Farm_Data_BH";
    if (G_CONST::err_d_c)  {
      erdat = datapath + "errors_daily.txt";
    }
    else  {
      erdat = datapath + "errors_cumulative.txt";
    }
  }
  ncnt = 0;
  nreg = G_CONST::fit_reg;
  enreg = 1;
  xmin = 0.0;
  xmax = 0.0;
  ymin = 0.0;
  ymax = 0.0;
}


void Farms::bigload()  { // hahaha.
  // why not just put in constructor? cos then part of Epi::Epi()...
  cout << "Loading farms: " + fmdat << endl;
  loadfarms();
  if (G_CONST::cp_cull)  {
    cout << "Loading CPs: " + cpdat << endl;
    loadcps();
  }
  cout << "Loading seeds: " + whdat << "\t" << flush;
  loadseeds(G_CONST::seedday); // (t=1) == (Feb 1st)
  cout <<   seedfarms.size() << endl;
  cout << "Loading metric: " + erdat << flush;
  loaderrdat();
  cout <<"!"<<endl;
}


/** Read in UK '01 demography data from my re-formatted file
 *  0: Region 0-Cumbria, 1-Devon, 2-Wales, 3-Scotland, 4-Rest of England
 *  1: X km
 *  2: Y km
 *  3: Cows
 *  4: Pigs
 *  5: Sheep
 *  then the inverses of these for within farm.... weird? could just calc when needed hahaha
 *
 * Line number (1:188496) used to cross ref the Demography to CPs IPs and DCs
 * Zero indexed here
 * Must ensure those 4 files match up. +vacc for Japan.
 */
void Farms::loadfarms()  {
  ifstream fdat;
  fdat.open(fmdat.data());
  if (fdat.is_open())  {
    ncnt = std::count(std::istreambuf_iterator<char>(fdat),
                      std::istreambuf_iterator<char>(),'\n'); // This many newlines==farms (I hope)
    N.resize(ncnt,vector<int>(G_CONST::fit_spc,0));
    N1.resize(ncnt,vector<double>(G_CONST::fit_spc,1.0));  // inverses... don't want accidental zero-div
    x.resize(ncnt,-1.0);      // default to off map?
    y.resize(ncnt,-1.0);
    region.resize(ncnt,0);    // default to all region 0, if fitting more, then read.
    eregion.resize(ncnt,-1);  // default to not a region, catch reg<0 to ignore state
    cps.resize(ncnt);         // farm ids are push_back'd in to this
    fdat.clear();             // reset eof and..
    fdat.seekg(0,ios::beg);   //   ... jump back to beginning of file
    double reg,n0,n1;         // temp storage vars to put in to the vectors
    double xx,yy,n10,n11;     // temp storage for location and inverse livestock numbers
    for (int i=0;i<ncnt;++i)  { // Takes first ncnt lines from file
      //fdat >> farms.region[i] >> farms.x[i] >> farms.y[i] >> farms.N[i][0] >> farms.N[i][1] >> farms.N[i][2] >> farms.N1[i][0] >> farms.N1[i][1] >> farms.N[i][2];
      fdat >> reg >> xx >> yy >> n0 >> n1 >> n10 >> n11 ;
      eregion[i]=int(round(reg)); x[i] = xx; y[i]=yy;
      N[i][0] = int(round(n0));
      N[i][1] = int(round(n1));
      N1[i][0] = n10;
      N1[i][1] = n11;
    }
    fdat.close();
    if (G_CONST::fit_reg==5)  {
      region = eregion;
    }
  }
  else {
    cout << "\n\nFarms.loadfarms() data error!\n\n";
    exit(-1);
  }
  // Farm location bounding box. floor/ceil'd in grid for nice grid layout
  xmin = 0.0;
  xmax = *max_element(x.begin(),x.end());
  ymin = 0.0;
  ymax = *max_element(y.begin(),y.end());

  vector<int> rtmp = eregion;
  sort(rtmp.begin(),rtmp.end());  // sort region and count changes for unique regions
  auto r_it = rtmp.begin();
  while (r_it!=rtmp.end()-1)  {
    if (*(r_it++)!=*r_it)  { // changes
      ++enreg;
    }
  }
}



/** \brief Read in CPs from sorted file
 * \return void
 */
void Farms::loadcps()  {
  ifstream cpsdat;
  cpsdat.open(cpdat.data());
  if (cpsdat.is_open())  {
    int junk = 0;
    for (int i=0;i<ncnt;++i)  { // TODO could clean up files, remove -1s and >188496s
      cps[i].reserve(20);
      cpsdat >> junk; // Discard: first input == line number i
      for (int j=0;j<20;++j)  {
        cpsdat >> junk;
        if ((junk<0)||(junk>ncnt))  {// ignore '-1' and >188496 imaginary farms
          continue;
        }
        else  {
          // Input file indexed from 1 (line numbers), not zero. Hence -1
          cps[i].push_back(junk-1);
        }
      }
    }
    cpsdat.close();
  }
  else  {
    cout << "\n\n NO CONTIGUITY DATA \n\n";
    exit(-1);
  }
  /*ofstream cptest("cptest");  // dump CPs as read/stored.
  for (int i=0;i<188496;++i)  {
    cptest << i;
    for (auto jj=cps[i].begin();jj!=cps[i].end();++jj)  {
      cptest << " " << *jj;
    }
    cptest << "\n";
  }
  cptest.close();*/
}

/** \brief Reads in real outbreak data for seeding simulations. Assumes data are time-sorted
 *  \param tmax int - date to seed epidemic
 *  \return void
 */
void Farms::loadseeds(int tmax)  {
  ifstream ipdat;
  ifstream dcdat;

  // TODO check JPN complies to format. UK dates are fine.
  ipdat.open(whdat.data()); // confirmed infected
  dcdat.open(cldat.data()); // pre-emptively culled
  // Here be the infected farms
  int seedf = -1;
  int seedt = 0;
  seedfarms.reserve(177); // Only reserving... not a big deal that this size is set
  seedtimes.reserve(177); //
  if (ipdat.is_open())  {
    while (1)  {
      ipdat >> seedf >> seedt;
      seedt += 4; // NOTE Rescale to t=1 == Feb 1st.  File has: t=2 == 6th Feb.
      if (seedt>tmax)  {
        break;
      }
      //cout << seedf << "\t" << seedt << endl;
      seedfarms.push_back(seedf-1); // File refers to line number in farm dat, not 0 indexed as here
      seedtimes.push_back(seedt);
    }
  }
  else  {
    cout << "\n\nFarms.loadfarms() data error - Which_Farms\n\n";
    exit(-1);
  }
  ipdat.close();
  // Here be the culled farms
  seedf = -1;
  seedt = 0;
  seedculls.reserve(2);
  seedtimec.reserve(2);
  if (dcdat.is_open())  {
    while (1)  {
      dcdat >> seedf >> seedt;
      seedt += 9; // NOTE Rescale to t=1 == Feb 1st. File: t=13 == 22nd Feb
      if (seedt>tmax)  {
        break;
      }
      //cout << seedf << "\t" << seedt << endl;
      seedculls.push_back(seedf-1); // File refers to line number in farm dat, not 0 indexed.
      seedtimec.push_back(seedt);
    }
  }
  else  {
    cout << "\n\nFarms.loadfarms() data error - CULLS_Farms \n\n";
    exit(-1);
  }
  dcdat.close();
  /*cout << endl; for (int i=0;i<seedfarms.size();++i)  {
    cout << seedtimes[i] << "\t" << seedfarms[i] << "\n";
  } cout << endl<<endl;*/
}


/** \brief Load daily cumulative errors from file
 *  Parameter acceptance depends on this (daily/regional?) fit metric
 */
void Farms::loaderrdat()  {
  // Storage for cumulative numbers of infections and culls. 01/02-23/12
  // reads in entirety of outbreak.
  // 326 obv for uk outbreak. lots of extra(empty) space when jpn
  int real_epi_t_end = 326; // FIXME hardcoded end of file ish
  farmsi = vector< vector<int> >(enreg,vector<int>(real_epi_t_end,0));
  farmsc = vector< vector<int> >(enreg,vector<int>(real_epi_t_end,0));
  cowssi = vector< vector<int> >(enreg,vector<int>(real_epi_t_end,0));
  cowssc = vector< vector<int> >(enreg,vector<int>(real_epi_t_end,0));
  sheepi = vector< vector<int> >(enreg,vector<int>(real_epi_t_end,0));
  sheepc = vector< vector<int> >(enreg,vector<int>(real_epi_t_end,0));

  ifstream rd;
  // erdat already defined at Farms::Farms()
  rd.open(erdat.data());
// set simulation length off data given? But then have to pre-gen new file...
// for tt=1 is 1st feb, this is where the file starts
  if (rd.is_open())  {
    for (unsigned int tt=1;tt<farmsi[0].size();++tt)  {
      for (int reg=0;reg<enreg;++reg)  {
        rd >> farmsi[reg][tt];  // Infected premises
      }
      for (int reg=0;reg<enreg;++reg)  {
        rd >> farmsc[reg][tt];  // Culled premises
      }
      for (int reg=0;reg<enreg;++reg)  {
        rd >> cowssi[reg][tt];  // Infected cows
      }
      for (int reg=0;reg<enreg;++reg)  {
        rd >> cowssc[reg][tt];  // Culled cows
      }
      for (int reg=0;reg<enreg;++reg)  {
        rd >> sheepi[reg][tt];  // Infected sheep
      }
      for (int reg=0;reg<enreg;++reg)  {
        rd >> sheepc[reg][tt];  // Culled sheep
      }
    }
    rd.close();
  }
  else  {
    cout << "Failed reading in cumulative outbreak data!" << endl;
    exit(-1);
  }
}



void Farms::dumpfarms()  {
  ofstream fdump("fdump");
  for (int i=0;i<ncnt;++i)  {
    fdump << region[i] << "\t" << x[i] << "\t" << y[i] << "\t" << N[i][0] << "\t" << N[i][1] << "\t" << N[i][2] << "\n";
  }
}





















