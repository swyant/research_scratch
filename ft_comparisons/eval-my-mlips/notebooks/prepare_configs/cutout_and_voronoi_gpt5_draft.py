"""
cutout + Voronoi interstitials + mapping back to original coords

Requirements:
  - pymatgen
  - pymatgen-analysis-defects (the VoronoiInterstitialGenerator lives here)
  - numpy

Usage example:
  python find_gb_interstitials.py \
    --structure big_supercell.vasp \
    --center-index 12345 \
    --box-size 15.0 \
    --min-sep 2.0 \
    --out-interstitials interstitials.vasp
"""

import numpy as np
from pymatgen.core import Structure, Lattice, Element
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
import argparse
import sys

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--structure", "-s", required=True,
                   help="Path to structure file (POSCAR/vasp/xyz etc).")
    p.add_argument("--center-index", "-i", required=True, type=int,
                   help="Index (0-based) of the grain-boundary atom to center on.")
    p.add_argument("--box-size", "-b", default=15.0, type=float,
                   help="Cubic box size in Å (default 15.0).")
    p.add_argument("--min-sep", default=2.0, type=float,
                   help="Minimum distance from any periodic boundary in cutout; atoms closer are removed (Å).")
    p.add_argument("--insert-elements", "-e", nargs="+", default=None,
                   help="List of element symbols to consider as inserted elements (e.g. Li H). If omitted, generator will be called without explicit elements.")
    p.add_argument("--out-interstitials", "-o", default=None,
                   help="If set, write POSCAR (VASP) with interstitial sites placed as single-atom structure for inspection.")
    return p.parse_args()

def nearest_image_cart_coords(structure: Structure, center_cart: np.ndarray):
    """
    For each site in `structure`, return the cartesian coordinates of the periodic-image
    that is closest to center_cart (minimum image). Also return fractional coords of that
    chosen image (useful if mapping back).
    """
    lattice = structure.lattice
    frac_center = lattice.get_fractional_coords(center_cart)
    coords_list = []
    frac_list = []
    for site in structure:
        f = np.array(site.frac_coords)
        # delta in fractional coordinates wrapped into [-0.5, 0.5)
        delta = f - frac_center
        delta = delta - np.rint(delta)   # subtract integer rounding
        chosen_frac = frac_center + delta
        chosen_cart = lattice.get_cartesian_coords(chosen_frac)
        coords_list.append(np.array(chosen_cart))
        frac_list.append(chosen_frac)
    return np.array(coords_list), np.array(frac_list)

def create_cutout_structure(structure: Structure, center_index: int,
                            box_size: float=15.0, min_sep: float=2.0):
    """
    Create a cubic cutout structure centered on atom `center_index` (0-based).
    The new Structure has lattice vectors forming an orthogonal box of `box_size` and
    positions translated so center atom is at box center (box_size/2, box_size/2, box_size/2).

    Atoms that end up within `min_sep` of any box face are removed (so when we re-enable
    periodic boundary conditions for the generator, there won't be overlapping images close
    to opposite faces).
    Returns:
      cut_struct: pymatgen Structure (with periodic=True)
      transform_info: dict with keys to map coordinates back:
         { 'center_cart': center_cart_in_original,
           'box_size': box_size }
    """
    # original lattice and positions
    center_cart = np.array(structure[center_index].coords)

    # compute nearest-image cart coords for all atoms
    coords_nearest, frac_chosen = nearest_image_cart_coords(structure, center_cart)

    half = box_size / 2.0
    # choose atoms whose cart coords fall within the cubic box centered at center_cart
    rel = coords_nearest - center_cart  # relative vector to center (may be negative)
    mask_in_box = np.all(np.abs(rel) <= half + 1e-8, axis=1)
    if not np.any(mask_in_box):
        raise RuntimeError("No atoms found within requested box. Check center index / box size.")

    indices_in_box = np.where(mask_in_box)[0]
    rel_in_box = rel[indices_in_box]

    # compute new cart coords inside box coordinates: shift so box origin at (0,0,0)
    # box center will be at (half, half, half)
    new_cart_coords = rel_in_box + half

    # remove atoms that are too close to box faces
    # i.e., if any coordinate < min_sep or > box_size - min_sep, drop them
    near_face = np.any((new_cart_coords < min_sep) | (new_cart_coords > (box_size - min_sep)), axis=1)
    kept_mask = ~near_face
    kept_indices = indices_in_box[kept_mask]
    new_cart_kept = new_cart_coords[kept_mask]

    if kept_indices.size == 0:
        raise RuntimeError("After trimming atoms near faces, no atoms remain in cutout. Try smaller min_sep or bigger box.")

    # Build new Structure: cubic/orthogonal lattice
    lattice = Lattice.cubic(box_size)
    species = [structure[i].specie for i in kept_indices]
    cut_struct = Structure(lattice, [s.symbol for s in species], new_cart_kept, coords_are_cartesian=True)
    # keep PBC True to use Voronoi generator (but we removed atoms near faces)
    cut_struct.set_pbc((True, True, True))

    transform_info = {
        "center_cart": center_cart.copy(),
        "box_size": box_size,
        "box_half": half
    }

    return cut_struct, transform_info

def run_voronoi_and_extract(cut_struct: Structure, insert_elements=None):
    """
    Run the VoronoiInterstitialGenerator on cut_struct.
    We try to be tolerant of the generator's return type:
      - It may return a list of Defect objects (with .site or .structure)
      - It may return structures directly
    We'll extract cartesian coordinates of proposed interstitial sites.

    Returns list of (element_symbol, cart_coord_in_cut_box_coords)
    """
    gen = VoronoiInterstitialGenerator()
    # Many examples show calling .generate(structure, insert_elems)
    # but the API might accept either; try both.
    if insert_elements:
        raw = gen.generate(cut_struct, insert_elements)
    else:
        raw = gen.generate(cut_struct)

    interstitial_coords = []

    # raw may be list of objects. Inspect and try to extract coordinates robustly.
    for item in raw:
        # item might be a pymatgen Structure (one-site) or a Defect object
        # We'll try a few heuristics
        extracted = False
        # 1) if item has attribute 'structure' (Defect-like), try to find the interstitial site
        if hasattr(item, "structure"):
            s = item.structure
            # try to find a site in s that is not in cut_struct by comparing coordinates
            # relaxed tolerance
            tol = 1e-3
            if len(s) == len(cut_struct) + 1:
                # find the site whose coords are not present in the original site list
                orig_coords = np.round(np.array(cut_struct.cart_coords) / tol).astype(int)
                new_coords = np.round(np.array([site.coords for site in s]) / tol).astype(int)
                # identify index in s that is not in orig_coords
                for idx, nc in enumerate(new_coords):
                    # boolean match against all orig
                    if not any(np.all(nc == oc) for oc in orig_coords):
                        # found the added interstitial
                        interstitial_coords.append((s[idx].species_string, np.array(s[idx].coords)))
                        extracted = True
                        break
            else:
                # fallback: if item has attribute 'site', take it
                if hasattr(item, "site"):
                    site = item.site
                    interstitial_coords.append((site.specie.symbol if hasattr(site, "specie") else str(site.species), np.array(site.coords)))
                    extracted = True

        if extracted:
            continue

        # 2) if item is a Structure with one extra site
        from pymatgen.core import Structure as PMGStructure
        if isinstance(item, PMGStructure):
            s = item
            # try to find site not present in cut_struct
            tol = 1e-3
            orig_coords = np.round(np.array(cut_struct.cart_coords) / tol).astype(int)
            new_coords = np.round(np.array([site.coords for site in s]) / tol).astype(int)
            found = False
            if len(s) == len(cut_struct) + 1:
                for idx, nc in enumerate(new_coords):
                    if not any(np.all(nc == oc) for oc in orig_coords):
                        interstitial_coords.append((s[idx].species_string, np.array(s[idx].coords)))
                        found = True
                        break
            if found:
                continue

        # 3) item may be a dict or tuple with 'site' or coords
        if isinstance(item, (list, tuple)) and len(item) >= 1:
            # maybe (structure, something) or (coords,)
            # find any numpy array-like inside
            for part in item:
                if hasattr(part, "__iter__") and len(part) == 3 and all(isinstance(x, (float, int)) for x in part):
                    interstitial_coords.append(("X", np.array(part)))
                    extracted = True
                    break
            if extracted:
                continue

        if isinstance(item, dict):
            for k in ("site", "coords", "position", "cart_coords"):
                if k in item:
                    c = np.array(item[k])
                    interstitial_coords.append((item.get("element", "X"), c))
                    extracted = True
                    break
            if extracted:
                continue

        # 4) last resort: string/repr parsing to find coords (unlikely)
        # Skip if nothing matched
    # Deduplicate near-duplicate coords (Voronoi often returns many near-identical sites)
    unique = []
    tol = 1e-3
    for el, c in interstitial_coords:
        found = False
        for (el2, c2) in unique:
            if np.linalg.norm(c - c2) < 1e-2:
                found = True
                break
        if not found:
            unique.append((el, c))
    return unique

def map_back_to_original(cart_in_cut: np.ndarray, transform_info: dict, wrap_to_original_cell=True, original_structure: Structure=None):
    """
    Map a cartesian coordinate given in cut-structure coordinates (0..box_size)
    back to the original full-system cart coords.

    In the cut creation we did:
       rel = coords_nearest - center_cart
       new_cart = rel + box_half  <- this is cart_in_cut

    So inverse:
       rel = cart_in_cut - box_half
       original_cart = rel + center_cart

    If wrap_to_original_cell=True and original_structure is provided, we map the resulting
    coordinate back into [0,1) fractional coordinates of original_structure.lattice.
    """
    half = transform_info["box_half"]
    center_cart = transform_info["center_cart"]
    rel = np.array(cart_in_cut) - half
    original_cart = rel + center_cart
    if wrap_to_original_cell and (original_structure is not None):
        lat = original_structure.lattice
        frac = lat.get_fractional_coords(original_cart)
        frac_wrapped = frac - np.floor(frac)  # bring into [0,1)
        original_cart = lat.get_cartesian_coords(frac_wrapped)
    return np.array(original_cart)

def main():
    args = parse_args()
    struct = Structure.from_file(args.structure)
    if args.center_index < 0 or args.center_index >= len(struct):
        print("center-index out of range (0..{})".format(len(struct)-1), file=sys.stderr)
        sys.exit(1)

    cut_struct, transform_info = create_cutout_structure(struct, args.center_index, box_size=args.box_size, min_sep=args.min_sep)
    print(f"Cutout created: {len(cut_struct)} atoms in a {args.box_size:.2f} A box (center atom index {args.center_index}).")

    # run generator
    inters = run_voronoi_and_extract(cut_struct, insert_elements=args.insert_elements)
    print(f"Voronoi generator returned {len(inters)} unique candidate interstitial(s) in cutout coords.")

    mapped = []
    for el, c in inters:
        orig_cart = map_back_to_original(c, transform_info, wrap_to_original_cell=True, original_structure=struct)
        mapped.append((el, orig_cart))
    # print results
    for i, (el, cart) in enumerate(mapped):
        print(f"Interstitial {i:3d}: element={el:3s}  cart (Å)= {cart[0]:8.4f} {cart[1]:8.4f} {cart[2]:8.4f}")

    # optionally write an output structure with all interstitials as single atoms
    if args.out_interstitials:
        # create structure with original lattice (so coordinates fall into original cell)
        lat = struct.lattice
        frac_coords = [lat.get_fractional_coords(cart) for _, cart in mapped]
        # use a dummy element (if not provided in inters) choose first provided or 'X'
        species = [el if el is not None and el != "X" else "X" for el, _ in mapped]
        out_struct = Structure(lat, species, frac_coords, coords_are_cartesian=False)
        out_struct.to(filename=args.out_interstitials)
        print(f"Wrote {len(mapped)} interstitials to {args.out_interstitials}")

if __name__ == "__main__":
    main()
