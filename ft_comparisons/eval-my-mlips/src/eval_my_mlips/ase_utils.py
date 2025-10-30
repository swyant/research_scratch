"""
mp_mace_minimize.py

Workflow:
  1) download conventional unit cell for an MP id (default mp-352)
  2) convert to ASE Atoms
  3) optionally relax (local minimization) with a chosen ASE calculator (MACE examples provided)
  4) create a supercell
  5) write XYZ and LAMMPS DATA files

Requirements:
  - ase
  - pymatgen
  - mace-torch (if using MACE calculators)
  - torch (if using MACE with CUDA)
  - (optionally) mp-api / pymatgen MPRester (we use pymatgen.ext.matproj.MPRester)

Usage example:
  export MAPI_KEY="your_materialsproject_api_key"
  python mp_mace_minimize.py --mp-id mp-352 --supercell 2 2 2 --out-prefix mp352_mace
"""

import os
import argparse
from pathlib import Path

# ASE, pymatgen, MACE imports
from ase.io import write as ase_write
from ase.optimize import BFGS, FIRE
from ase import Atoms

# pymatgen API
from pymatgen.ext.matproj import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor

# Provide optional imports for mace - import lazily to allow running the script without MACE installed
MACE_AVAILABLE = True
try:
    # two common patterns used in MACE examples:
    #  - MACECalculator (full class)
    #  - convenience factories like mace_mp or mace_off
    from mace.calculators import MACECalculator, mace_mp, mace_off  # noqa: F401
except Exception:
    MACE_AVAILABLE = False


def get_api_key_from_env_or_arg(env_var="MAPI_KEY", arg_key=None):
    """
    Prefer environment variable (common practice), else fallback to passed argument.
    Raises RuntimeError if not found.
    """
    if arg_key:
        return arg_key
    api = os.environ.get(env_var, None)
    if api:
        return api
    raise RuntimeError(
        f"Materials Project API key not found. Set environment variable {env_var} or pass --api-key."
    )


def fetch_conventional_structure(mp_id: str, api_key: str):
    """
    Fetch structure for mp_id using pymatgen's MPRester,
    then convert to *conventional* standard cell using SpacegroupAnalyzer.
    Returns a pymatgen Structure object.
    """
    with MPRester(api_key) as m:
        # get the structure that MP stores for this id
        struct = m.get_structure_by_material_id(mp_id)
    # get conventional standard cell (this relies on pymatgen's symmetry analyzer)
    sga = SpacegroupAnalyzer(struct, symprec=1e-3)
    conv = sga.get_conventional_standard_structure()
    return conv


def pymatgen_to_ase(structure):
    """Convert pymatgen Structure to ase.Atoms"""
    return AseAtomsAdaptor.get_atoms(structure)


def setup_mace_calculator(choice: str = "mace_mp", *,
                          model_path: str = None,
                          device: str = "cpu",
                          extra_kwargs: dict = None):
    """
    Return an ASE calculator object configured for MACE.

    Parameters:
      - choice: 'mace_mp', 'mace_off', or 'mace_modelpath'
          - 'mace_mp' : convenience medium/large model interface (if available)
          - 'mace_off': another convenience factory (example only)
          - 'mace_modelpath': use explicit serialized model via MACECalculator(model_path=...)
      - model_path: path to a saved MACE model file (required for 'mace_modelpath')
      - device: 'cpu' or 'cuda'
      - extra_kwargs: dict of other kwargs to pass to the calculator

    NOTE: This function tries not to import MACE until called.
    """
    if not MACE_AVAILABLE:
        raise RuntimeError("MACE is not installed in this Python environment. Install mace-torch first.")

    extra_kwargs = extra_kwargs or {}
    choice = choice.lower()
    if choice == "mace_mp":
        # convenience: choose built-in medium model (depends on MACE installation)
        return mace_mp(model="medium", device=device, **extra_kwargs)
    if choice == "mace_off":
        return mace_off(model="medium", device=device, **extra_kwargs)
    if choice == "mace_modelpath":
        if model_path is None:
            raise ValueError("model_path is required for choice='mace_modelpath'")
        # The exact class name can differ by mace versions; many versions expose MACECalculator.
        # MACECalculator takes model_path and device in current docs.
        return MACECalculator(model_path=model_path, device=device, **extra_kwargs)
    raise ValueError("Unknown MACE calculator choice: " + choice)


def relax_atoms(atoms: Atoms, fmax: float = 0.05, optimizer: str = "BFGS", maxsteps: int = 200):
    """
    Run a geometry optimization on `atoms`. Uses ASE optimizers (BFGS or FIRE).
    Returns the optimizer object after running.
    """
    if optimizer.lower() == "bfgs":
        opt = BFGS(atoms, logfile=None)
    else:
        opt = FIRE(atoms, logfile=None)
    opt.run(fmax=fmax, steps=maxsteps)
    return opt


def make_supercell(atoms: Atoms, repeat_tuple=(2, 2, 2)):
    """Return new Atoms object repeated by repeat_tuple (nx,ny,nz)."""
    nx, ny, nz = repeat_tuple
    return atoms.repeat((nx, ny, nz))


def write_outputs(atoms: Atoms, prefix: str = "out"):
    p_xyz = Path(f"{prefix}.xyz")
    p_data = Path(f"{prefix}.data")

    # write xyz
    ase_write(str(p_xyz), atoms, format="xyz")
    # write lammps data -- specify 'lammps-data' format and a typical style
    # Note: for some systems you may want to pass additional keywords required by ASE's LAMMPS-data writer.
    ase_write(str(p_data), atoms, format="lammps-data", lammps_keyword="atomic", atom_style="atomic")
    return p_xyz, p_data


def main(args):
    api_key = get_api_key_from_env_or_arg(arg_key=args.api_key)
    print(f"Using Materials Project API key from environment/argument (masked).")

    print(f"Fetching conventional unit cell for {args.mp_id} ...")
    pmg_struct = fetch_conventional_structure(args.mp_id, api_key)
    print("Fetched structure. Converted to conventional standard cell.")

    ase_atoms = pymatgen_to_ase(pmg_struct)
    print(f"Converted to ASE Atoms. N_atoms = {len(ase_atoms)}")

    # center / periodic check (optional)
    # ase_atoms.center()  # not necessary for periodic solids

    # Attach calculator if requested (and MACE present)
    if args.use_mace:
        try:
            calc = setup_mace_calculator(choice=args.mace_choice,
                                         model_path=args.model_path,
                                         device=args.device,
                                         extra_kwargs={})
        except Exception as e:
            raise RuntimeError(f"Failed to create MACE calculator: {e}")
        ase_atoms.set_calculator(calc)
        print(f"Attached MACE calculator ({args.mace_choice}). Device={args.device}")

        # run minimization
        print("Running geometry optimization (local relaxation) ...")
        relax_atoms(ase_atoms, fmax=args.fmax, optimizer=args.optimizer, maxsteps=args.maxsteps)
        print("Relaxation done.")
    else:
        print("No calculator attached; skipping relaxation step (use --use-mace to enable).")

    # Make supercell
    sc = make_supercell(ase_atoms, repeat_tuple=(args.nx, args.ny, args.nz))
    print(f"Created supercell with repeat {args.nx}x{args.ny}x{args.nz}. N_atoms_total = {len(sc)}")

    # Write outputs
    xyz_path, data_path = write_outputs(sc, prefix=f"{args.out_prefix}")
    print(f"Wrote XYZ -> {xyz_path}")
    print(f"Wrote LAMMPS DATA -> {data_path}")
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch MP structure, relax with MACE, create supercell, write outputs")
    parser.add_argument("--api-key", default=None, help="Materials Project API key (optional; default reads $MAPI_KEY)")
    parser.add_argument("--mp-id", default="mp-352", help="Materials Project id (default mp-352)")
    parser.add_argument("--use-mace", action="store_true", help="Attach a MACE ASE calculator and relax")
    parser.add_argument("--mace-choice", default="mace_mp", choices=["mace_mp", "mace_off", "mace_modelpath"],
                        help="Which MACE flavor to use (mace_mp/mace_off convenience factories, or mace_modelpath)")
    parser.add_argument("--model-path", default=None, help="Path to MACE model (required if mace_choice=mace_modelpath)")
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"], help="Device for MACE ('cpu' or 'cuda')")
    parser.add_argument("--fmax", default=0.05, type=float, help="Force tolerance (eV/Å) for relax")
    parser.add_argument("--optimizer", default="BFGS", choices=["BFGS", "FIRE"], help="ASE optimizer to use")
    parser.add_argument("--maxsteps", default=200, type=int, help="Max optimizer steps")
    parser.add_argument("--nx", default=2, type=int, help="Supercell repeat in x")
    parser.add_argument("--ny", default=2, type=int, help="Supercell repeat in y")
    parser.add_argument("--nz", default=2, type=int, help="Supercell repeat in z")
    parser.add_argument("--out-prefix", default="mp_structure", help="Prefix for output files (prefix.xyz, prefix.data)")

    args = parser.parse_args()
    main(args)

