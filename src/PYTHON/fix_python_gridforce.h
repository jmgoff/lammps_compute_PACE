/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(python/gridforce,FixPythonGridForce);
FixStyle(python,FixPythonGridForce);
// clang-format on
#else

#ifndef LMP_FIX_PYTHON_GRIDFORCE_H
#define LMP_FIX_PYTHON_GRIDFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPythonGridForce : public Fix {
 public:
  FixPythonGridForce(class LAMMPS *, int, char **);
  ~FixPythonGridForce() override;
  int setmask() override;
  void setup();
  
  int quadraticflag;
  int twojmax, switchflag, bzeroflag, bnormflag, bikflag, dgradflag;
  int chemflag, wselfallflag;
  int switchinnerflag;

  //void allocate();                     // allocate pairstyle arrays
  //void allocate_grid();                // create grid arrays
  //void deallocate_grid();              // free grid arrays
  //void grid2x(int, int, int, double*); // convert global indices to coordinates
  //void set_grid_global();              // set global grid
  //void set_grid_local();               // set bounds for local grid
  //void assign_coords();                // assign coords for grid
  //void copy_gridlocal_to_local_array();// copy 4d gridlocal array to 2d local array
  //void settings(int, char **);
  //void grid_settings(int, char **);
  //void compute(int, int);
  //void init_list(int, class NeighList *);
  //void init_sna();
  int ntypes;

  void end_of_step() override;
  void post_force(int) override;
  void pre_force(int) override;

  double rcutfac;
  double *radelem;
  double *wjelem;
  double **cutsq;
  int nxlo, nxhi, nylo, nyhi, nzlo, nzhi; // local grid bounds, inclusive
  int nx, ny, nz;                      // global grid dimensions
  double **alocal;                     // pointer to local array
  //double **aall;                     // pointer to global array
  double ****gridlocal;                // local grid, redundant w.r.t. alocal
  //int ngridlocal;                      // number of local grid points
  //int ndesc;                           // number of descriptors
  //int ndesc_base;                      // number of columns used for coords, etc.

 private:
  void allocate();                     // allocate pairstyle arrays
  void allocate_grid();                // create grid arrays
  void deallocate_grid();              // free grid arrays
  void grid2x(int, int, int, double*); // convert global indices to coordinates
  void set_grid_global();              // set global grid
  void set_grid_local();               // set bounds for local grid
  void assign_coords();                // assign coords for grid
  void copy_gridlocal_to_local_array();// copy 4d gridlocal array to 2d local array
  void settings(int, char **);
  void grid_settings(int, char **);
  void compute(int, int);
  void init_list(int, class NeighList *);
  void init_sna();
  void *lmpPtr;
  void *pFunc;
  bigint lasttime;
  int selected_callback;
  //int nx, ny, nz;                      // global grid dimensions
  //int nxlo, nxhi, nylo, nyhi, nzlo, nzhi; // local grid bounds, inclusive
  int ngridlocal;                      // number of local grid points
  int nvalues;                         // number of values per grid point
  //double ****gridlocal;                // local grid, redundant w.r.t. alocal
  //double **alocal;                     // pointer to Compute::array_local
  int triclinic;                       // triclinic flag
  double *boxlo, *prd;                 // box info (units real/ortho or reduced/tri)
  double *sublo, *subhi;               // subdomain info (units real/ortho or reduced/tri)
  double delxinv,delyinv,delzinv;      // inverse grid spacing
  double delx,dely,delz;               // grid spacing
  int nargbase;                        // number of base class args
  double cutmax;                       // largest cutoff distance
  int ndesc;                           // number of descriptors
  int ndesc_base;                      // number of columns used for coords, etc.
  int gridlocal_allocated;             // shows if gridlocal allocated
  double **beta;                       // betas for all local grid points in list
  int beta_max;                        // length of beta
  /*
  void allocate();                     // allocate pairstyle arrays
  void allocate_grid();                // create grid arrays
  void deallocate_grid();              // free grid arrays
  void grid2x(int, int, int, double*); // convert global indices to coordinates
  void set_grid_global();              // set global grid
  void set_grid_local();               // set bounds for local grid
  void assign_coords();                // assign coords for grid
  void copy_gridlocal_to_local_array();// copy 4d gridlocal array to 2d local array
  //void compute_beta();                 // get betas from someplace
  */
  int ncoeff;
  int nelements;
  class NeighList *list;
  class SNA *snaptr;
  double *sinnerelem;
  double *dinnerelem;
  double rfac0, rmin0;
  int *map;    // map types to [0,nelements)
  Fix *test_fix;

};

}    // namespace LAMMPS_NS

#endif
#endif
