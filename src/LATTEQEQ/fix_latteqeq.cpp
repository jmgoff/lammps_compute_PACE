// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Christian Negre (LANL)
------------------------------------------------------------------------- */

#include "fix_latteqeq.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "compute_pace.h"

#include <cstring>
#include <ctime>

using namespace LAMMPS_NS;
using namespace FixConst;

extern "C" {
  void latteqeq(int *, int *, double *, int *, int *,
             double *, double *, double *, double *,
             double *, double *, double *, int *,
             double *, double *, double *, double *, 
	     int *, //double *,  double *, 
	     double *,
	     int * , bool *);
  int latteqeq_abiversion();
}

// the ABIVERSION number here must be kept consistent
// with its counterpart in the LATTE library and the
// prototype above. We want to catch mismatches with
// a meaningful error messages, as they can cause
// difficult to debug crashes or memory corruption.

#define LATTEQEQ_ABIVERSION 20180622

/* ---------------------------------------------------------------------- */

FixLatteqeq::FixLatteqeq(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Must use units metal with fix latteqeq command");

  if (comm->nprocs != 1)
    error->all(FLERR,"Fix latteqeq currently runs only in serial");

  if (LATTEQEQ_ABIVERSION != latteqeq_abiversion())
    error->all(FLERR,"LAMMPS is linked against incompatible LATTE library");

  if (narg != 4) error->all(FLERR,"Illegal fix latteqeq command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  virial_global_flag = 1;
  thermo_energy = thermo_virial = 1;

  // store ID of compute pe/atom used to generate Coulomb potential for LATTE
  // null pointer means LATTE will compute Coulombic potential

  coulomb = 0;
  id_pe = nullptr;

  if (strcmp(arg[3],"NULL") != 0) {
    coulomb = 1;
    error->all(FLERR,"Fix latteqeq does not yet support a LAMMPS calculation of a Coulomb potential");

    id_pe = utils::strdup(arg[3]);
    c_pe = modify->get_compute_by_id(id_pe);
    if (!c_pe) error->all(FLERR,"Could not find fix latteqeq compute ID {}", id_pe);
    if (c_pe->peatomflag == 0)
      error->all(FLERR,"Fix latteqeq compute ID does not compute pe/atom");
  }
  
  // Sep6, latteqeq shoud only provide coulomb force
  //coulomb = 1; // added on Sep 6

  // initializations

  nmax = 0;
  currentstep = -1;
  qpotential = nullptr;
  //qxlbo = nullptr;
  //kernel = nullptr;
  flatte = nullptr;
  gradx = nullptr;

  latte_energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixLatteqeq::~FixLatteqeq()
{
  delete[] id_pe;
  memory->destroy(qpotential);
  //memory->destroy(qxlbo);
  //memory->destroy(kernel);
  memory->destroy(flatte);
  memory->destroy(gradx);
}

/* ---------------------------------------------------------------------- */

int FixLatteqeq::setmask()
{
  int mask = 0;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLatteqeq::init()
{
  // error checks

  if (domain->dimension == 2)
    error->all(FLERR,"Fix latte requires 3d problem");

  if (coulomb) {
    if (atom->q_flag == 0 || force->pair == nullptr || force->kspace == nullptr)
      error->all(FLERR,"Fix latte cannot compute Coulomb potential");

    c_pe = modify->get_compute_by_id(id_pe);
    if (!c_pe) error->all(FLERR,"Could not find fix latte compute ID {}", id_pe);
  }

  // must be fully periodic or fully non-periodic

  if (domain->nonperiodic == 0) pbcflag = 1;
  else if (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic)
    pbcflag = 0;
  else error->all(FLERR,"Fix latte requires 3d simulation");

  // create qpotential & flatte if needed
  // for now, assume nlocal will never change

  /*if (qxlbo == nullptr) {
    memory->create(qxlbo,atom->nlocal,8, "latte:qxlbo");
  }
  if (kernel == nullptr) {
    memory->create(kernel,atom->nlocal,atom->nlocal,"latte:kernel");
  }*/

  // the first row is charge negativities 
  if (gradx == nullptr) {
    memory->create(gradx,atom->nlocal*3+1,atom->nlocal,"latte:gradx");
  }


  if (coulomb && qpotential == nullptr) {
    memory->create(qpotential,atom->nlocal,"latte:qpotential");
    memory->create(flatte,atom->nlocal,3,"latte:flatte");
  }
}

/* ---------------------------------------------------------------------- */

void FixLatteqeq::init_list(int /*id*/, NeighList * /*ptr*/)
{
  // list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixLatteqeq::setup(int vflag)
{
  newsystem = 1;
  post_force(vflag);
  newsystem = 0;
}

/* ---------------------------------------------------------------------- */

void FixLatteqeq::min_setup(int vflag)
{
  newsystem = 1;
  post_force(vflag);
  newsystem = 0;
}

/* ---------------------------------------------------------------------- */

void FixLatteqeq::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixLatteqeq::initial_integrate(int /*vflag*/) {}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to tally per-atom energies
------------------------------------------------------------------------- */

void FixLatteqeq::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
}

/* ---------------------------------------------------------------------- */

void FixLatteqeq::post_force(int vflag)
{
  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  // compute Coulombic potential = pe[i]/q[i]
  // invoke compute pe/atom
  // wrap with clear/add and trigger pe/atom calculation every step

  std::clock_t c_start = std::clock();

  if (coulomb) {
    modify->clearstep_compute();

    if (!(c_pe->invoked_flag & Compute::INVOKED_PERATOM)) {
      c_pe->compute_peratom();
      c_pe->invoked_flag |= Compute::INVOKED_PERATOM;
    }

    modify->addstep_compute(update->ntimestep+1);

    double *pe = c_pe->vector_atom;
    double *q = atom->q;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (q[i]) qpotential[i] = pe[i]/q[i];
      else qpotential[i] = 0.0;
  }

  // hardwire these unsupported flags for now

  int coulombflag = 0;
  neighflag = 0;

  // set flags used by LATTE
  // NOTE: LATTE does not compute per-atom energies or virials

  int flags[6];

  flags[0] = pbcflag;         // 1 for fully periodic, 0 for fully non-periodic
  flags[1] = coulombflag;     // 1 for LAMMPS computes Coulombics, 0 for LATTE
  flags[2] = eflag_atom;      // 1 to return per-atom energies, 0 for no
  flags[3] = vflag_global && thermo_virial;    // 1 to return global/per-atom
  flags[4] = vflag_atom && thermo_virial;      //   virial, 0 for no
  flags[5] = neighflag;       // 1 to pass neighbor list to LATTE, 0 for no

  // setup LATTE arguments

  int natoms = atom->nlocal;
  int dbdr_len = 3*natoms * natoms;
  double *coords = &atom->x[0][0];
  int *type = atom->type;
  int ntypes = atom->ntypes;
  double *mass = &atom->mass[1];
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *forces;
  int ncoef; 
  bool latteerror = false;
  if (coulomb) forces = &flatte[0][0];
  else forces = &atom->f[0][0];
  int maxiter = -1;

  for (int l = 0; l < atom->nlocal*3+1; l++){
    for (int m = 0; m < atom->nlocal; m++){
        gradx[l][m] =0.0;
    }
  }
 
// dbdr_legnth: 3N*N
// dBdR = lmp_pace[N;(N+dbdr_length), 3:(nd+3)]
// 
// for k in range(nd):
// for l in range(dbdr_length):
// i = force_indices[l,0]
//
  for (int ic = 0; ic < modify->ncompute; ic++){
    if (strcmp(modify->compute[ic]->style,"pace") == 0 && currentstep >=0){
        ncoef = modify->compute[ic]->size_array_cols - 3;
        if (modify->compute[ic]->local_flag){
        	for (int j = 0; j < atom->nlocal; j++){
        	  for (int k = 3; k < modify->compute[ic]->size_local_cols; k++){
        	     gradx[0][j] += modify->compute[ic]->array[j][k];;
        	  }
        	  //printf("%d %f \n", j, gradx[0][j]);
        	}

		// get the gradients
		//printf("get electronagtivities gradients\n");
		for (int ia=natoms; ia < natoms + dbdr_len; ia++){
		   int i = modify->compute[ic]->array[ia][0];
		   int j = modify->compute[ic]->array[ia][1];
		   int a = modify->compute[ic]->array[ia][2];
		   for (int k = 3; k<modify->compute[ic]->size_array_cols; k++){
		      gradx[3*j+a+1][i] += modify->compute[ic]->array[ia][k];
		      //printf("i=%d j=%d a=%d  %f\n", i, j, a, gradx[3*j+a+1][i]);
		   }
		}
        }
        else{
        	for (int j = 0; j < atom->nlocal; j++){
        	  for (int k = 1; k < modify->compute[ic]->size_local_cols+1; k++){
        	     gradx[0][j] += modify->compute[ic]->array[j][k];;
        	  }
        	  printf("%d %f \n", j, gradx[0][j]);
        	}
        }
    }
  }

  //printf("currentstep= %d\n", currentstep);
  if (currentstep >= 0){
	  latteqeq(flags,&natoms,coords,type,&ntypes,mass,boxlo,boxhi,&domain->xy,
			  &domain->xz,&domain->yz,forces,&maxiter,&latte_energy,
			  &atom->v[0][0],&update->dt,virial,
			  &currentstep, //&qxlbo[0][0],&kernel[0][0],
			  &gradx[0][0], &newsystem,&latteerror);
  }
  currentstep += 1;

  if (latteerror) error->all(FLERR,"Internal LATTE problem");

  // sum LATTE forces to LAMMPS forces
  // e.g. LAMMPS may compute Coulombics at some point

  if (coulomb) {
    double **f = atom->f;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) {
      f[i][0] += flatte[i][0];
      f[i][1] += flatte[i][1];
      f[i][2] += flatte[i][2];
    }
  }
  //printf("test-zy: leaving fix_latteqeq:post_force\n");
  std::clock_t c_end = std::clock();
  double time_elapsed_ms = (c_end-c_start) / CLOCKS_PER_SEC;
  //printf("CPU time used= %f s\n", time_elapsed_ms);
}

/* ---------------------------------------------------------------------- */

void FixLatteqeq::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixLatteqeq::final_integrate() {}

/* ---------------------------------------------------------------------- */

void FixLatteqeq::reset_dt()
{
  //dtv = update->dt;
  //dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   DFTB energy from LATTE
------------------------------------------------------------------------- */

double FixLatteqeq::compute_scalar()
{
  return latte_energy;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double FixLatteqeq::memory_usage()
{
  double bytes = 0.0;
  if (coulomb) bytes += (double)nmax * sizeof(double);
  if (coulomb) bytes += (double)nmax*3 * sizeof(double);
  return bytes;
}
