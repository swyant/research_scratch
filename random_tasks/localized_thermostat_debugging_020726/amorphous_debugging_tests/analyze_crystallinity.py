"""Analyze crystallinity of LAMMPS dump trajectories.

Distinguishes between:
  - A crystalline structure with thermal vibrations (atoms oscillate around lattice sites)
  - A truly amorphous/disordered structure (atoms have drifted far from lattice sites)

Metrics computed:
  1. Per-atom displacement from initial (reference) positions — the most direct measure.
     In a crystal at high T, displacements are small (~0.1-0.3 Å). If atoms have moved
     >1 Å from their initial sites, the structure is becoming disordered.
  2. Radial distribution function g(r) — crystalline structures show sharp peaks;
     amorphous structures show broad, smeared features.
  3. Lindemann ratio — ratio of RMS displacement to nearest-neighbor distance.
     Lindemann criterion: melting occurs when this exceeds ~0.1-0.15.

Usage:
    uv run python amorphous_debugging_tests/analyze_crystallinity.py <dump_file> [--ref-frame 0] [--data-file path]

If --data-file is provided, displacements are computed relative to the original DATA file
positions (true lattice sites). Otherwise, frame 0 of the dump is used as reference.
"""

import argparse
import sys
import os
import numpy as np

# Add parent dir to path so we can import the dump parser
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from localized_thermostat.diagnostics.dump_parser import parse_lammpstrj, get_frame_count


def parse_data_file_positions(path):
    """Extract atom positions from a LAMMPS DATA file (atom_style atomic).

    Returns:
        ids: array of atom IDs (1-based)
        positions: (N, 3) array of positions
        box: (3, 2) array of box bounds [[xlo,xhi],[ylo,yhi],[zlo,zhi]]
    """
    ids = []
    positions = []
    box = np.zeros((3, 2))

    in_atoms = False
    n_atoms = 0
    box_idx = 0

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                if in_atoms and len(positions) >= n_atoms:
                    break
                continue

            if 'atoms' in line and 'atom' not in line.split()[1:]:
                # "4861 atoms"
                parts = line.split()
                if parts[1] == 'atoms':
                    n_atoms = int(parts[0])
                continue

            if 'xlo' in line:
                parts = line.split()
                box[0] = [float(parts[0]), float(parts[1])]
                continue
            if 'ylo' in line:
                parts = line.split()
                box[1] = [float(parts[0]), float(parts[1])]
                continue
            if 'zlo' in line:
                parts = line.split()
                box[2] = [float(parts[0]), float(parts[1])]
                continue

            if line.startswith('Atoms'):
                in_atoms = True
                continue

            if in_atoms:
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        atom_id = int(parts[0])
                        # atom_style atomic: id type x y z
                        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                        ids.append(atom_id)
                        positions.append([x, y, z])
                    except (ValueError, IndexError):
                        continue

    return np.array(ids), np.array(positions), box


def compute_displacements_mic(pos, ref_pos, box_bounds):
    """Compute displacements with minimum image convention."""
    delta = pos - ref_pos
    box_lengths = box_bounds[:, 1] - box_bounds[:, 0]
    for i in range(3):
        delta[:, i] -= box_lengths[i] * np.round(delta[:, i] / box_lengths[i])
    return delta


def compute_rdf(positions, box_bounds, r_max=8.0, n_bins=200, sample_fraction=0.2):
    """Compute radial distribution function g(r).

    Uses a random subsample of atoms as centers for efficiency.
    """
    box_lengths = box_bounds[:, 1] - box_bounds[:, 0]
    n_atoms = len(positions)

    # Subsample centers for speed
    n_centers = max(50, int(n_atoms * sample_fraction))
    center_indices = np.random.choice(n_atoms, size=min(n_centers, n_atoms), replace=False)

    dr = r_max / n_bins
    r_edges = np.linspace(0, r_max, n_bins + 1)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    hist = np.zeros(n_bins)

    volume = np.prod(box_lengths)
    rho = n_atoms / volume

    for ci in center_indices:
        delta = positions - positions[ci]
        for i in range(3):
            delta[:, i] -= box_lengths[i] * np.round(delta[:, i] / box_lengths[i])
        dists = np.linalg.norm(delta, axis=1)
        # Exclude self
        dists = dists[dists > 1e-10]
        counts, _ = np.histogram(dists, bins=r_edges)
        hist += counts

    # Normalize
    n_centers_used = len(center_indices)
    for i in range(n_bins):
        r_lo, r_hi = r_edges[i], r_edges[i + 1]
        shell_vol = (4.0 / 3.0) * np.pi * (r_hi**3 - r_lo**3)
        ideal_count = rho * shell_vol
        if ideal_count > 0:
            hist[i] /= (n_centers_used * ideal_count)

    return r_centers, hist


def analyze_frame(frame_data, ref_positions, ref_ids, box_bounds, nn_distance=3.2):
    """Analyze a single frame for crystallinity metrics.

    Args:
        frame_data: dict from parse_lammpstrj
        ref_positions: reference (lattice) positions, indexed by atom ID
        ref_ids: atom IDs for reference positions
        box_bounds: (3,2) box bounds
        nn_distance: nearest-neighbor distance for Lindemann ratio (Hf-Hf ~ 3.2 Å)

    Returns dict of metrics.
    """
    columns = frame_data['columns']
    data = frame_data['data']

    id_col = columns.index('id')
    x_col = columns.index('x')
    y_col = columns.index('y')
    z_col = columns.index('z')

    dump_ids = data[:, id_col].astype(int)
    dump_pos = data[:, [x_col, y_col, z_col]]

    # Build ID→position map for reference
    ref_map = {int(aid): ref_positions[i] for i, aid in enumerate(ref_ids)}

    # Match dump atoms to reference by ID
    matched_ref = []
    matched_dump = []
    for i, aid in enumerate(dump_ids):
        aid = int(aid)
        if aid in ref_map:
            matched_ref.append(ref_map[aid])
            matched_dump.append(dump_pos[i])

    matched_ref = np.array(matched_ref)
    matched_dump = np.array(matched_dump)

    # Displacements from reference
    delta = compute_displacements_mic(matched_dump, matched_ref, box_bounds)
    disp_magnitudes = np.linalg.norm(delta, axis=1)

    rms_disp = np.sqrt(np.mean(disp_magnitudes**2))
    max_disp = np.max(disp_magnitudes)
    mean_disp = np.mean(disp_magnitudes)
    median_disp = np.median(disp_magnitudes)

    # Lindemann ratio
    lindemann = rms_disp / nn_distance

    # Fraction of atoms displaced > threshold
    frac_gt_05 = np.mean(disp_magnitudes > 0.5)
    frac_gt_1 = np.mean(disp_magnitudes > 1.0)
    frac_gt_2 = np.mean(disp_magnitudes > 2.0)

    # RDF
    r_centers, g_r = compute_rdf(matched_dump, box_bounds)

    # RDF crystallinity score: ratio of first peak height to first minimum
    # Crystalline: sharp peaks → high ratio; amorphous: broad → low ratio
    # Look for first peak in range 2.5-4.0 Å (Hf-Hf NN)
    mask_peak = (r_centers > 2.5) & (r_centers < 4.0)
    mask_min = (r_centers > 4.0) & (r_centers < 5.5)

    first_peak = np.max(g_r[mask_peak]) if np.any(mask_peak) else 1.0
    first_min = np.min(g_r[mask_min]) if np.any(mask_min) else 1.0
    rdf_sharpness = first_peak / max(first_min, 0.01)

    return {
        'n_atoms_analyzed': len(matched_dump),
        'rms_displacement': rms_disp,
        'mean_displacement': mean_disp,
        'median_displacement': median_disp,
        'max_displacement': max_disp,
        'lindemann_ratio': lindemann,
        'frac_displaced_gt_0.5A': frac_gt_05,
        'frac_displaced_gt_1.0A': frac_gt_1,
        'frac_displaced_gt_2.0A': frac_gt_2,
        'rdf_first_peak_height': first_peak,
        'rdf_first_min_depth': first_min,
        'rdf_sharpness_ratio': rdf_sharpness,
        'r_centers': r_centers,
        'g_r': g_r,
        'displacements': disp_magnitudes,
    }


def print_diagnosis(metrics, timestep):
    """Print human-readable diagnosis."""
    print(f"\n{'='*60}")
    print(f"  Crystallinity Analysis — Timestep {timestep}")
    print(f"{'='*60}")
    print(f"  Atoms analyzed:           {metrics['n_atoms_analyzed']}")
    print(f"")
    print(f"  --- Displacement from reference ---")
    print(f"  RMS displacement:         {metrics['rms_displacement']:.4f} Å")
    print(f"  Mean displacement:        {metrics['mean_displacement']:.4f} Å")
    print(f"  Median displacement:      {metrics['median_displacement']:.4f} Å")
    print(f"  Max displacement:         {metrics['max_displacement']:.4f} Å")
    print(f"  Lindemann ratio:          {metrics['lindemann_ratio']:.4f}  (melting ~ 0.10-0.15)")
    print(f"")
    print(f"  Fraction displaced > 0.5 Å: {metrics['frac_displaced_gt_0.5A']:.1%}")
    print(f"  Fraction displaced > 1.0 Å: {metrics['frac_displaced_gt_1.0A']:.1%}")
    print(f"  Fraction displaced > 2.0 Å: {metrics['frac_displaced_gt_2.0A']:.1%}")
    print(f"")
    print(f"  --- RDF sharpness ---")
    print(f"  First peak height:        {metrics['rdf_first_peak_height']:.2f}")
    print(f"  First minimum depth:      {metrics['rdf_first_min_depth']:.2f}")
    print(f"  Sharpness ratio:          {metrics['rdf_sharpness_ratio']:.2f}  (crystal > 5, amorphous < 3)")
    print(f"")

    # Diagnosis
    lindemann = metrics['lindemann_ratio']
    rdf_sharp = metrics['rdf_sharpness_ratio']
    frac_1A = metrics['frac_displaced_gt_1.0A']

    if lindemann < 0.10 and frac_1A < 0.05:
        verdict = "CRYSTALLINE — normal thermal vibrations"
    elif lindemann < 0.15 and frac_1A < 0.20:
        verdict = "BORDERLINE — significant disorder, near melting"
    else:
        verdict = "AMORPHOUS/MELTED — structure has lost crystallinity"

    print(f"  VERDICT: {verdict}")
    print(f"{'='*60}")


def save_rdf_plot(all_metrics, output_path):
    """Save RDF comparison plot across frames."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not available, skipping RDF plot)")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # RDF plot
    ax = axes[0]
    for label, metrics in all_metrics.items():
        ax.plot(metrics['r_centers'], metrics['g_r'], label=label, alpha=0.8)
    ax.set_xlabel('r (Å)')
    ax.set_ylabel('g(r)')
    ax.set_title('Radial Distribution Function')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 8)
    ax.axhline(1.0, color='gray', ls='--', alpha=0.5)

    # Displacement histogram
    ax = axes[1]
    for label, metrics in all_metrics.items():
        ax.hist(metrics['displacements'], bins=50, alpha=0.5, label=label, density=True)
    ax.set_xlabel('Displacement from reference (Å)')
    ax.set_ylabel('Probability density')
    ax.set_title('Displacement Distribution')
    ax.axvline(0.5, color='orange', ls='--', label='0.5 Å threshold')
    ax.axvline(1.0, color='red', ls='--', label='1.0 Å threshold')
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"  Plot saved: {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Analyze crystallinity of LAMMPS dump trajectory')
    parser.add_argument('dump_file', help='Path to .lammpstrj dump file')
    parser.add_argument('--data-file', default=None,
                        help='LAMMPS DATA file for reference positions (default: use frame 0 of dump)')
    parser.add_argument('--frames', default='0,-1',
                        help='Comma-separated frame indices to analyze (default: "0,-1")')
    parser.add_argument('--all-frames', action='store_true',
                        help='Analyze all frames')
    parser.add_argument('--nn-distance', type=float, default=3.2,
                        help='Nearest-neighbor distance in Å for Lindemann ratio (default: 3.2 for Hf)')
    parser.add_argument('--output-plot', default=None,
                        help='Path for output plot (default: <dump_dir>/crystallinity_analysis.png)')

    args = parser.parse_args()

    dump_path = args.dump_file
    n_frames = get_frame_count(dump_path)
    print(f"Dump file: {dump_path}")
    print(f"Total frames: {n_frames}")

    # Get reference positions
    if args.data_file:
        print(f"Reference: DATA file {args.data_file}")
        ref_ids, ref_positions, box_bounds = parse_data_file_positions(args.data_file)
    else:
        print(f"Reference: frame 0 of dump")
        ref_frame = parse_lammpstrj(dump_path, frame=0)
        columns = ref_frame['columns']
        id_col = columns.index('id')
        x_col = columns.index('x')
        y_col = columns.index('y')
        z_col = columns.index('z')
        ref_ids = ref_frame['data'][:, id_col].astype(int)
        ref_positions = ref_frame['data'][:, [x_col, y_col, z_col]]
        box_bounds = ref_frame['box_bounds']

    # Determine which frames to analyze
    if args.all_frames:
        frame_indices = list(range(n_frames))
    else:
        frame_indices = []
        for s in args.frames.split(','):
            s = s.strip()
            idx = int(s)
            if idx < 0:
                idx = n_frames + idx
            frame_indices.append(idx)

    # Analyze each frame
    all_metrics = {}
    for fi in frame_indices:
        frame = parse_lammpstrj(dump_path, frame=fi)
        metrics = analyze_frame(frame, ref_positions, ref_ids, box_bounds, args.nn_distance)
        label = f"frame {fi} (step {frame['timestep']})"
        all_metrics[label] = metrics
        print_diagnosis(metrics, frame['timestep'])

    # Save plot
    if args.output_plot is None:
        dump_dir = os.path.dirname(dump_path) or '.'
        args.output_plot = os.path.join(dump_dir, 'crystallinity_analysis.png')

    save_rdf_plot(all_metrics, args.output_plot)

    # Summary table
    print(f"\n{'='*80}")
    print(f"  Summary Table")
    print(f"{'='*80}")
    print(f"  {'Frame':<25} {'RMS(Å)':>8} {'Max(Å)':>8} {'Lindemann':>10} {'>1Å':>8} {'RDF sharp':>10}")
    print(f"  {'-'*25} {'-'*8} {'-'*8} {'-'*10} {'-'*8} {'-'*10}")
    for label, m in all_metrics.items():
        print(f"  {label:<25} {m['rms_displacement']:>8.3f} {m['max_displacement']:>8.3f} "
              f"{m['lindemann_ratio']:>10.4f} {m['frac_displaced_gt_1.0A']:>7.1%} "
              f"{m['rdf_sharpness_ratio']:>10.2f}")


if __name__ == '__main__':
    main()
