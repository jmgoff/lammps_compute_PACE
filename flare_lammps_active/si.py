from lammps import lammps
import numpy as np

import ase
import ase.io
from ase.calculators.lammpsrun import LAMMPS

from flare.bffs.sgp._C_flare import SparseGP, NormalizedDotProduct, B2
from flare.bffs.sgp.sparse_gp import optimize_hyperparameters
from flare.learners.lmpotf import LMPOTF

import os

os.system("rm -rf sidft && mkdir sidft")

n = 8
l = 4
rcut = 3.77118
n_species = 1

# sparse GP with noise hyperparameters for
# energy, force and stress
sigma_e = 0.002 * 64
sigma_f = 0.05
sigma_s = 0.003
kernel = NormalizedDotProduct(2.0, 2)  # sigma, power
sparse_gp = SparseGP([kernel], sigma_e, sigma_f, sigma_s)
sparse_gp.Kuu_jitter = 1e-8

lmp_path = '/home/jmgoff/Software/new_flare/lammps_compute_PACE'

os.environ["ASE_LAMMPSRUN_COMMAND"] = "%s/build/lmp" %lmp_path
dftcalc = LAMMPS(
    specorder=["Si"],
    keep_tmp_files=True,
    tmp_dir="tmp",
    keep_alive=False,
    **{
        "pair_style": "sw",
        "pair_coeff": [f"* * %s/potentials/Si.sw Si" % lmp_path]
    }
)

import sys
sys.argv = ["in.si"]
#import wandb
#wandb.init(anonymous="allow")

def opt_hyps(lmpotf, lmp, step):
    return 10 <= lmpotf.dft_calls <= 40 and lmpotf.dft_calls % 5 == 0

lmpotf = LMPOTF(
    sparse_gp=sparse_gp,
    descriptors=B2("chebyshev", "quadratic", [0.0, rcut], [], [n_species, n, l]),
    rcut=rcut,
    type2number=14,
    dftcalc=dftcalc,
    dft_call_threshold=0.0005,
    dft_add_threshold=0.0001,
    energy_correction=0, #QE: -125.1, # JDFTx: -102.6,
    #wandb=wandb,
    wandb=None,
    std_xyz_fname="sidft/std.*.xyz",
    dft_xyz_fname="sidft/dft.*.xyz",
    model_fname="si.otf.flare",
    log_fname="otf.si.log",
    hyperparameter_optimization=opt_hyps,
)

otf = lmpotf.step
