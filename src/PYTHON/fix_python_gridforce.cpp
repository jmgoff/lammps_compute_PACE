// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: James Goff (Sandia National Laboratories)
------------------------------------------------------------------------- */
#ifdef PYTHON_GRIDFORCE
#include "fix_python_gridforce.h"

#include "error.h"
#include "lmppython.h"
#include "python_compat.h"
#include "python_utils.h"
#include "update.h"

#include "modify.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include "tokenizer.h"

//#include "pair_grid.h"
//#include "pair_sna_grid.h"
#include "sna.h"

#include <cstring>
#include <Python.h>   // IWYU pragma: export

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPythonGridForce::FixPythonGridForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), gridlocal(nullptr), alocal(nullptr), beta(nullptr), radelem(nullptr), wjelem(nullptr), test_fix(nullptr), cutsq(nullptr)
{
  ntypes = atom->ntypes;
  int nargmin_sna = 6 + 2 * ntypes;
  printf("narg: %d  | nargmin: %d \n", narg,nargmin_sna+6);
  if (narg < 6 + nargmin_sna) error->all(FLERR,"Illegal fix python/gridforce command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix python/gridforce command");

  settings( narg, arg);

  rcutfac = utils::numeric(FLERR, arg[10], false, lmp);
  rfac0 = utils::numeric(FLERR, arg[11], false, lmp);
  twojmax = utils::inumeric(FLERR, arg[12], false, lmp);
  //memory->create(radelem,ntypes+1,"fix:python/gridforce:radelem"); // offset by 1 to match up with types
  //memory->create(wjelem,ntypes+1,"fix:python/gridforce:wjelem");

  /*
  printf("ntypes before %d \n",ntypes);
  for (int ii = 0; ii < ntypes; ii++){
    printf("args for types %s \n" , arg[13+ii]);
    radelem[ii + 1] = utils::numeric(FLERR, arg[12 + ii], false, lmp);
    //radelem[ii] = utils::numeric(FLERR, arg[13 + ii], false, lmp);
  }
  for (int ii = 0; ii < ntypes; ii++){
    wjelem[ii + 1] = utils::numeric(FLERR, arg[12 + ntypes + ii], false, lmp);
    //wjelem[ii] = utils::numeric(FLERR, arg[13 + ntypes + ii], false, lmp);
  }
  */
  // ensure Python interpreter is initialized
  python->init();

  if (strcmp(arg[4],"post_force") == 0) {
    selected_callback = POST_FORCE;
  } else if (strcmp(arg[4],"end_of_step") == 0) {
    selected_callback = END_OF_STEP;
  } else if (strcmp(arg[4],"pre_force") == 0) {
    selected_callback = PRE_FORCE;
  } else {
    error->all(FLERR,"Unsupported callback name for fix python/gridforce");
  }

  // get Python function
  PyUtils::GIL lock;

  PyObject *pyMain = PyImport_AddModule("__main__");

  if (!pyMain) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not initialize embedded Python");
  }

  char *fname = arg[5];
  pFunc = PyObject_GetAttrString(pyMain, fname);

  if (!pFunc) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not find Python function");
  }

  lmpPtr = PY_VOID_POINTER(lmp);
  lasttime=-1;

  centroidstressflag = CENTROID_NOTAVAIL;

  // construct cutsq

  /*
  double cut;
  double cutmaxa = 0.0;
  //memory->create(cutsq, ntypes + 1, ntypes + 1, "sna/atom:cutsq");
  for (int i = 1; i <= ntypes; i++) {
    cut = 2.0 * radelem[i] * rcutfac;
    if (cut > cutmaxa) cutmaxa = cut;
    cutsq[i][i] = cut * cut;
    for (int j = i + 1; j <= ntypes; j++) {
      cut = (radelem[i] + radelem[j]) * rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut * cut;
    }
  }
  */
  // set local input checks

  int sinnerflag = 0;
  int dinnerflag = 0;

  // process optional args

  int iarg = nargmin_sna+6;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "rmin0") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      rmin0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "switchflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      switchflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bzeroflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      bzeroflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "quadraticflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      quadraticflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "chem") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      chemflag = 1;
      memory->create(map, ntypes + 1, "compute_sna_grid:map");
      nelements = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      for (int i = 0; i < ntypes; i++) {
        int jelem = utils::inumeric(FLERR, arg[iarg + 2 + i], false, lmp);
        if (jelem < 0 || jelem >= nelements) error->all(FLERR, "Illegal compute {} command", style);
        map[i + 1] = jelem;
      }
      iarg += 2 + ntypes;
    } else if (strcmp(arg[iarg], "bnormflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      bnormflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "wselfallflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      wselfallflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bikflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      bikflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dgradflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Illegal compute snap command");
      dgradflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchinnerflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Illegal compute snap command");
      switchinnerflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "sinner") == 0) {
      iarg++;
      if (iarg + ntypes > narg) error->all(FLERR, "Illegal compute {} command", style);
      memory->create(sinnerelem, ntypes + 1, "snap:sinnerelem");
      for (int i = 0; i < ntypes; i++)
        sinnerelem[i + 1] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
      sinnerflag = 1;
      iarg += ntypes;
    } else if (strcmp(arg[iarg], "dinner") == 0) {
      iarg++;
      if (iarg + ntypes > narg) error->all(FLERR, "Illegal compute {} command", style);
      memory->create(dinnerelem, ntypes + 1, "snap:dinnerelem");
      for (int i = 0; i < ntypes; i++)
        dinnerelem[i + 1] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
      dinnerflag = 1;
      iarg += ntypes;
    } else
      iarg +=1;
      //error->all(FLERR, "Illegal compute {} command", style);
  }

  if (switchinnerflag && !(sinnerflag && dinnerflag))
    error->all(
        FLERR,
        "Illegal compute {} command: switchinnerflag = 1, missing sinner/dinner keyword",
        style);

  if (!switchinnerflag && (sinnerflag || dinnerflag))
    error->all(
        FLERR,
        "Illegal compute {} command: switchinnerflag = 0, unexpected sinner/dinner keyword",
        style);

  if (dgradflag && !bikflag)
    error->all(FLERR,"Illegal compute snap command: dgradflag=1 requires bikflag=1");

  if (dgradflag && quadraticflag)
    error->all(FLERR,"Illegal compute snap command: dgradflag=1 not implemented for quadratic SNAP");

  snaptr = new SNA(lmp, rfac0, twojmax, rmin0, switchflag, bzeroflag, chemflag, bnormflag,wselfallflag, nelements, switchinnerflag);
  ndesc = 0;
  ngridlocal = 0;

  ndesc_base = 6;
  gridlocal_allocated = 0;
  beta_max = 0;
  beta = nullptr;

  //snaptr = nullptr;
}

/* ---------------------------------------------------------------------- */

FixPythonGridForce::~FixPythonGridForce()
{
  if (copymode) return;
  PyUtils::GIL lock;
  Py_CLEAR(lmpPtr);
  deallocate_grid();
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;

  if (chemflag) memory->destroy(map);
}

/* ---------------------------------------------------------------------- */

void FixPythonGridForce::init_sna()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style python/gridforce requires newton pair on");

  snaptr = new SNA(lmp, rfac0, twojmax,
                   rmin0, switchflag, bzeroflag,
                   chemflag, bnormflag, wselfallflag,
           nelements, switchinnerflag);
  ncoeff = snaptr->ncoeff;
  ndesc = ndesc_base + ncoeff;
  snaptr->init();

}
/* ---------------------------------------------------------------------- */

void FixPythonGridForce::grid2x(int ix, int iy, int iz, double *x)
{
  x[0] = ix*delx;
  x[1] = iy*dely;
  x[2] = iz*delz;

  if (triclinic) domain->lamda2x(x, x);
}

/* ---------------------------------------------------------------------- */
void FixPythonGridForce::allocate_grid()
{
  if (nxlo <= nxhi && nylo <= nyhi && nzlo <= nzhi) {
    gridlocal_allocated = 1;
    memory->create4d_offset(gridlocal,ndesc,nzlo,nzhi,nylo,nyhi,
                            nxlo,nxhi,"python/gridforce");
    memory->create(alocal, ngridlocal, ndesc, "python/gridforce");
    //memory->create(beta, ngridlocal, ndesc-ndesc_base, "pair/grid:beta");
  }
}
/* ---------------------------------------------------------------------- */
void FixPythonGridForce::set_grid_global()
{
  // calculate grid layout

  triclinic = domain->triclinic;

  if (triclinic == 0) {
    prd = domain->prd;
    boxlo = domain->boxlo;
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    prd = domain->prd_lamda;
    boxlo = domain->boxlo_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];

  delxinv = nx/xprd;
  delyinv = ny/yprd;
  delzinv = nz/zprd;

  delx = 1.0/delxinv;
  dely = 1.0/delyinv;
  delz = 1.0/delzinv;
}
/* ---------------------------------------------------------------------- */
void FixPythonGridForce::set_grid_local()
{
  // nx,ny,nz = extent of global grid
  // indices into the global grid range from 0 to N-1 in each dim
  // if grid point is inside my sub-domain I own it,
  //   this includes sub-domain lo boundary but excludes hi boundary
  // ixyz lo/hi = inclusive lo/hi bounds of global grid sub-brick I own
  // if proc owns no grid cells in a dim, then ilo > ihi
  // if 2 procs share a boundary a grid point is exactly on,
  //   the 2 equality if tests insure a consistent decision
  //   as to which proc owns it

  double xfraclo,xfrachi,yfraclo,yfrachi,zfraclo,zfrachi;

  if (comm->layout != Comm::LAYOUT_TILED) {
    xfraclo = comm->xsplit[comm->myloc[0]];
    xfrachi = comm->xsplit[comm->myloc[0]+1];
    yfraclo = comm->ysplit[comm->myloc[1]];
    yfrachi = comm->ysplit[comm->myloc[1]+1];
    zfraclo = comm->zsplit[comm->myloc[2]];
    zfrachi = comm->zsplit[comm->myloc[2]+1];
  } else {
    xfraclo = comm->mysplit[0][0];
    xfrachi = comm->mysplit[0][1];
    yfraclo = comm->mysplit[1][0];
    yfrachi = comm->mysplit[1][1];
    zfraclo = comm->mysplit[2][0];
    zfrachi = comm->mysplit[2][1];
  }

  nxlo = static_cast<int> (xfraclo * nx);
  if (1.0*nxlo != xfraclo*nx) nxlo++;
  nxhi = static_cast<int> (xfrachi * nx);
  if (1.0*nxhi == xfrachi*nx) nxhi--;

  nylo = static_cast<int> (yfraclo * ny);
  if (1.0*nylo != yfraclo*ny) nylo++;
  nyhi = static_cast<int> (yfrachi * ny);
  if (1.0*nyhi == yfrachi*ny) nyhi--;

  nzlo = static_cast<int> (zfraclo * nz);
  if (1.0*nzlo != zfraclo*nz) nzlo++;
  nzhi = static_cast<int> (zfrachi * nz);
  if (1.0*nzhi == zfrachi*nz) nzhi--;

  ngridlocal = (nxhi - nxlo + 1) * (nyhi - nylo + 1) * (nzhi - nzlo + 1);
  printf("nxlo, nxhi, nylo, nyhi, nzlo, nzhi, %d %d %d %d %d %d \n",nxlo,nxhi,nylo,nyhi,nzlo,nzhi);
}
/* ---------------------------------------------------------------------- */

void FixPythonGridForce::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}


/* ----------------------------------------------------------------------
   copy coords to local array
------------------------------------------------------------------------- */

void FixPythonGridForce::assign_coords()
{
  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++){
    for (int iy = nylo; iy <= nyhi; iy++){
      for (int ix = nxlo; ix <= nxhi; ix++) {
        printf("igrid,iy' %d, %f \n",igrid,iy);
        alocal[igrid][0] = ix;
        alocal[igrid][1] = iy;
        alocal[igrid][2] = iz;
        double xgrid[3];
        grid2x(ix, iy, iz, xgrid);
        alocal[igrid][3] = xgrid[0];
        alocal[igrid][4] = xgrid[1];
        alocal[igrid][5] = xgrid[2];
        igrid++;
      }
    }
  }
}
/* ---------------------------------------------------------------------- */
void FixPythonGridForce::copy_gridlocal_to_local_array()
{
  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
        for (int icol = ndesc_base; icol < ndesc; icol++)
          alocal[igrid][icol] = gridlocal[icol][iz][iy][ix];
        igrid++;
      }
}

/* ---------------------------------------------------------------------- */
void FixPythonGridForce::deallocate_grid()
{
  if (gridlocal_allocated) {
    gridlocal_allocated = 0;
    memory->destroy4d_offset(gridlocal,nzlo,nylo,nxlo);
    memory->destroy(alocal);
    memory->destroy(beta);
  }
}

/* ---------------------------------------------------------------------- */
void FixPythonGridForce::setup()
{
  deallocate_grid();
  allocate_grid();
  set_grid_global();
  set_grid_local();
  //allocate_grid();
  assign_coords();
}

/* ---------------------------------------------------------------------- */

int FixPythonGridForce::setmask()
{
  return selected_callback;
}

/* ---------------------------------------------------------------------- */

void FixPythonGridForce::compute(int eflag, int vflag)
{
  double fij[3];
  //TODO comment init_sna if needed
  //init_sna();
  setup();
  ev_init(eflag,vflag);

  // compute sna for each gridpoint

  double** const x = atom->x;
  double **f = atom->f;
  const int* const mask = atom->mask;
  int * const type = atom->type;
  const int ntotal = atom->nlocal + atom->nghost;

  //auto py_beta = (PyObject *) beta;
  // insure rij, inside, and typej are of size ntotal
  
  snaptr->grow_rij(ntotal);

  // first generate fingerprint,
  // which allows calculation of beta
  
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	double xgrid[3];
	grid2x(ix, iy, iz, xgrid);
	const double xtmp = xgrid[0];
	const double ytmp = xgrid[1];
	const double ztmp = xgrid[2];

	// currently, all grid points are type 1
	
	const int itype = 1;
	int ielem = 0;
	if (chemflag)
	  ielem = map[itype];
	const double radi = radelem[itype];

	// rij[][3] = displacements between atom I and those neighbors
	// inside = indices of neighbors of I within cutoff
	// typej = types of neighbors of I within cutoff

	int ninside = 0;
	for (int j = 0; j < ntotal; j++) {

	  const double delx = xtmp - x[j][0];
	  const double dely = ytmp - x[j][1];
	  const double delz = ztmp - x[j][2];
	  const double rsq = delx*delx + dely*dely + delz*delz;
	  int jtype = type[j];
	  int jelem = 0;
	  if (chemflag)
	    jelem = map[jtype];
	  if (rsq < cutsq[jtype][jtype] && rsq > 1e-20) {
	    snaptr->rij[ninside][0] = delx;
	    snaptr->rij[ninside][1] = dely;
	    snaptr->rij[ninside][2] = delz;
	    snaptr->inside[ninside] = j;
	    snaptr->wj[ninside] = wjelem[jtype];
	    snaptr->rcutij[ninside] = 2.0*radelem[jtype]*rcutfac;
	    if (chemflag) snaptr->element[ninside] = jelem; // element index for multi-element snap
	    ninside++;
	  }
	}

	snaptr->compute_ui(ninside, ielem);
	snaptr->compute_zi();
	snaptr->compute_bi(ielem);

	// linear contributions

	for (int icoeff = 0; icoeff < ncoeff; icoeff++)
	  gridlocal[ndesc_base+icoeff][iz][iy][ix] = 
	    snaptr->blist[icoeff];

	// quadratic contributions
	// untested

	if (quadraticflag) {
	  int ncount = ncoeff;
	  for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
	    double bveci = snaptr->blist[icoeff];
	    gridlocal[ndesc_base+ncount++][iz][iy][ix] = 
	      0.5*bveci*bveci;
	    for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++)
	      gridlocal[ndesc_base+ncount++][iz][iy][ix] = 
		bveci*snaptr->blist[jcoeff];
	  }
	}
      }

  // this is a proxy for a call to the energy model
  // beta is dE/dB^i, the derivative of the total
  // energy w.r.t. to descriptors of grid point i
  
  //TODO get betas here or pass them in in preforce step
  //compute_beta();
  
  // second compute forces using beta
  
  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	double xgrid[3];
	grid2x(ix, iy, iz, xgrid);
	const double xtmp = xgrid[0];
	const double ytmp = xgrid[1];
	const double ztmp = xgrid[2];

	// currently, all grid points are type 1
	
	const int itype = 1;
	int ielem = 0;
	if (chemflag)
	  ielem = map[itype];
	const double radi = radelem[itype];

	// rij[][3] = displacements between atom I and those neighbors
	// inside = indices of neighbors of I within cutoff
	// typej = types of neighbors of I within cutoff

	int ninside = 0;
	for (int j = 0; j < ntotal; j++) {

	  const double delx = xtmp - x[j][0];
	  const double dely = ytmp - x[j][1];
	  const double delz = ztmp - x[j][2];
	  const double rsq = delx*delx + dely*dely + delz*delz;
	  int jtype = type[j];
	  int jelem = 0;
	  jelem = map[jtype];

	  if (rsq < cutsq[jtype][jtype] && rsq > 1e-20) {
	    snaptr->rij[ninside][0] = delx;
	    snaptr->rij[ninside][1] = dely;
	    snaptr->rij[ninside][2] = delz;
	    snaptr->inside[ninside] = j;
	    snaptr->wj[ninside] = wjelem[jtype];
	    snaptr->rcutij[ninside] = 2.0*radelem[jtype]*rcutfac;
	    if (switchinnerflag) {
	      snaptr->sinnerij[ninside] = 0.5*(sinnerelem[ielem]+sinnerelem[jelem]);
	      snaptr->dinnerij[ninside] = 0.5*(dinnerelem[ielem]+dinnerelem[jelem]);
	    }
	    if (chemflag) snaptr->element[ninside] = jelem;
	    ninside++;
	  }
	}

	// compute Ui, Yi for atom I
	
	if (chemflag)
	  snaptr->compute_ui(ninside, ielem);
	else
	  snaptr->compute_ui(ninside, 0);

	// for neighbors of I within cutoff:
	// compute Fij = dEi/dRj = -dEi/dRi
	// add to Fi, subtract from Fj
	// scaling is that for type I

	snaptr->compute_yi(beta[igrid]);

	for (int jj = 0; jj < ninside; jj++) {
	  int j = snaptr->inside[jj];
	  snaptr->compute_duidrj(jj);

	  snaptr->compute_deidrj(fij);

	  f[j][0] += fij[0];
	  f[j][1] += fij[1];
	  f[j][2] += fij[2];

	  // tally per-atom virial contribution
      /*
	  if (vflag)
	    ev_tally_xyz(-1,j,atom->nlocal,force->newton_pair,0.0,0.0,
			 fij[0],fij[1],fij[2],
			 -snaptr->rij[jj][0],-snaptr->rij[jj][1],
			 -snaptr->rij[jj][2]);
    */
	}

	// tally energy contribution

	if (eflag) {

	  // get descriptors again

	  snaptr->compute_zi();
	  snaptr->compute_bi(ielem);
	  
	  // evdwl = energy of atom I, sum over coeffs_k * Bi_k

	  double evdwl = 0.0;

	  // E = beta.B 

	  for (int icoeff = 0; icoeff < ncoeff; icoeff++)
	    evdwl += beta[igrid][icoeff]*snaptr->blist[icoeff];
	  
	  //ev_tally_full(-1,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);

	}
	igrid++;
      }

  //if (vflag_fdotr) virial_fdotr_compute();
  
}


/* ---------------------------------------------------------------------- */

//void FixPythonGridForce::py_dEdB(int ncoef)
//{
//}

/* ---------------------------------------------------------------------- */

void FixPythonGridForce::pre_force(int vflag)
//void FixPythonGridForce::pre_force( vflag)
{
  if (update->ntimestep % nevery) return;
  if (update->ntimestep == lasttime) return;

  //test_fix =  modify->get_fix_by_style("python/gridforce");
  //printf("in lammps: this fix %s \n",beta);
  // set up ACE/SNAP gradient arrays
  // NOTE THAT THIS IS WHERE ALL GRID AND SNAP INFO NEEDS TO BE COLLECTED


  PyUtils::GIL lock;
  //PyObject * dEdB = PyObject_CallFunction((PyObject*)pFunc, (char *)"O", (PyObject*)lmpPtr);
  //PyObject * result = PyObject_CallFunction((PyObject*)pFunc, (char *)"O", (PyObject*)lmpPtr);
  //double lmp_betas = PyFloat_AsDouble(py_betas);
  //PyObject * result = PyObject_CallFunction((PyObject*)pFunc, (char *)"O", (PyObject*)lmpPtr);
  //char fmtpre[] = "Oi";
  char fmtpre[] = "O";
  PyObject * dEdB = PyObject_CallFunction((PyObject*)pFunc, fmtpre, (PyObject*)lmpPtr, vflag);
  double  lmp_betas = PyFloat_AsDouble(dEdB);
  printf("lammps betas %f \n",lmp_betas);
  //printf("in pre force- result: %s \n",result);
  // result should provide dEdB
  //TODO comment out compute here if needed
  compute(1,0);
  if (!dEdB) {
    //PyUtils::Print_Errors();
    //error->all(FLERR,"Fix python/gridforce pre_force() method failed");
    printf("WARNING: dEdB is %s\n",dEdB);
  }

  Py_CLEAR(dEdB);
}

/* ---------------------------------------------------------------------- */

void FixPythonGridForce::end_of_step()
{
  PyUtils::GIL lock;

  PyObject * result = PyObject_CallFunction((PyObject*)pFunc, (char *)"O", (PyObject*)lmpPtr);

  if (!result) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/gridforce end_of_step() method failed");
  }

  Py_CLEAR(result);
}

/* ---------------------------------------------------------------------- */

void FixPythonGridForce::post_force(int vflag)
{
  if (update->ntimestep % nevery != 0) return;
  if (update->ntimestep == lasttime) return;
  lasttime = update->ntimestep;

  PyUtils::GIL lock;
  char fmt[] = "Oi";

  PyObject * result = PyObject_CallFunction((PyObject*)pFunc, fmt, (PyObject*)lmpPtr, vflag);

  if (!result) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/gridforce post_force() method failed");
  }

  Py_CLEAR(result);
}



/* ----------------------------------------------------------------------
   grid settings
------------------------------------------------------------------------- */
void FixPythonGridForce::grid_settings(int narg, char ** arg)
{
  if (narg < 4) error->all(FLERR,"Illegal pair style command");
  int iarg0 = 0;
  int iarg = iarg0 + 6;
  printf("gridsettings iarg, arg %d %s \n", iarg, arg[iarg]);
  if (strcmp(arg[iarg],"grid") == 0) {
    if (iarg+4 > narg) error->all(FLERR,"Illegal pair grid command");
    nx = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
    ny = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
    nz = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
    if (nx <= 0 || ny <= 0 || nz <= 0)
      error->all(FLERR,"All grid/local dimensions must be positive");
    iarg += 4;
  } else error->all(FLERR,"Illegal pair grid command");
  nargbase = iarg - iarg0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void FixPythonGridForce::settings(int narg, char ** arg)
{


  FixPythonGridForce::grid_settings(narg, arg);
    
  // skip over arguments used by base class
  // so that argument positions are identical to
  // regular per-atom compute
  
  arg += nargbase;
  narg -= nargbase;

  int ntypes = atom->ntypes;
  int nargmin = 3+2*ntypes;

  if (narg < nargmin) error->all(FLERR,"Illegal pair python/gridforce command");

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  bikflag = 0;
  dgradflag = 0;
  quadraticflag = 0;
  chemflag = 0;
  bnormflag = 0;
  quadraticflag = 0;
  wselfallflag = 0;
  switchinnerflag = 0;
  nelements = 1;
  
  // process required arguments

  memory->create(radelem,ntypes+1,"fix:python/gridforce:radelem"); // offset by 1 to match up with types
  memory->create(wjelem,ntypes+1,"fix:python/gridforce:wjelem");

  rcutfac = atof(arg[0]);
  rfac0 = atof(arg[1]);
  twojmax = atoi(arg[2]);


  printf("ntypes before %d narg %d \n",ntypes,narg);
  for (int itst = 0; itst < narg; itst ++){
    printf("itst %d | arg[itst] %s \n", itst,arg[itst]);
  }
  for (int ii = 0; ii < ntypes; ii++){
    printf("args for types %s \n" , arg[13+ii]);
    radelem[ii + 1] = utils::numeric(FLERR, arg[12 + ii], false, lmp);
    //radelem[ii] = utils::numeric(FLERR, arg[13 + ii], false, lmp);
  }
  for (int ii = 0; ii < ntypes; ii++){
    wjelem[ii + 1] = utils::numeric(FLERR, arg[12 + ntypes + ii], false, lmp);
    //wjelem[ii] = utils::numeric(FLERR, arg[13 + ntypes + ii], false, lmp);
  }
  /*
  for(int i = 0; i < ntypes; i++)
    radelem[i+1] = atof(arg[3+i]);
  for(int i = 0; i < ntypes; i++)
    wjelem[i+1] = atof(arg[3+ntypes+i]);
  */
  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq,ntypes+1,ntypes+1,"fix:python/gridforce:cutsq");
  for(int i = 1; i <= ntypes; i++) {
    cut = 2.0*radelem[i]*rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut*cut;
    for(int j = i+1; j <= ntypes; j++) {
      cut = (radelem[i]+radelem[j])*rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }

  // set local input checks

  int sinnerflag = 0;
  int dinnerflag = 0;

  // process optional args

  int iarg = nargmin;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rmin0") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair python/gridforce command");
      rmin0 = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair python/gridforce command");
      switchflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bzeroflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair python/gridforce command");
      bzeroflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"quadraticflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair python/gridforce command");
      quadraticflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"chem") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair python/gridforce command");
      chemflag = 1;
      memory->create(map,ntypes+1,"fix:python/gridforce:map");
      nelements = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      for (int i = 0; i < ntypes; i++) {
        int jelem = utils::inumeric(FLERR,arg[iarg+2+i],false,lmp);
        if (jelem < 0 || jelem >= nelements)
          error->all(FLERR,"Illegal pair python/gridforce command");
        map[i+1] = jelem;
      }
      iarg += 2+ntypes;
    } else if (strcmp(arg[iarg],"bnormflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair python/gridforce command");
      bnormflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"wselfallflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair python/gridforce command");
      wselfallflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchinnerflag") == 0) {
      if (iarg+2 > narg)
	error->all(FLERR,"Illegal pair python/gridforce command");
      switchinnerflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"sinner") == 0) {
      iarg++;
      if (iarg+ntypes > narg)
	error->all(FLERR,"Illegal pair python/gridforce command");
      memory->create(sinnerelem,ntypes+1,"snap:sinnerelem");
      for (int i = 0; i < ntypes; i++)
        sinnerelem[i+1] = utils::numeric(FLERR,arg[iarg+i],false,lmp);
      sinnerflag = 1;
      iarg += ntypes;
    } else if (strcmp(arg[iarg],"dinner") == 0) {
      iarg++;
      if (iarg+ntypes > narg)
	error->all(FLERR,"Illegal pair python/gridforce command");
      memory->create(dinnerelem,ntypes+1,"snap:dinnerelem");
      for (int i = 0; i < ntypes; i++)
        dinnerelem[i+1] = utils::numeric(FLERR,arg[iarg+i],false,lmp);
      dinnerflag = 1;
      iarg += ntypes;
    } else error->all(FLERR,"Illegal pair python/gridforce command");

  }

}

#endif
