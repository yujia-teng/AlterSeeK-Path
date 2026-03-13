#!/usr/bin/env python3
"""
Batch IBZ Centroid Calculator (hybrid version)
Usage:
    python batch_centroid_hybrid.py <directory> [--output summary.csv]
    python batch_centroid_hybrid.py file1.vasp file2.cif
"""
import os, sys, glob, csv, traceback
from compute_centroid_hybrid import run

def process_directory(directory, pattern="*.vasp", output_csv=None):
    files = sorted(glob.glob(os.path.join(directory, '**', pattern), recursive=True))
    if not files:
        for ext in ['*.cif', '*.POSCAR', '*.poscar', '*POSCAR*']:
            files = sorted(glob.glob(os.path.join(directory, '**', ext), recursive=True))
            if files: break
    if not files:
        print(f"No structure files found in {directory}"); return []

    print(f"Found {len(files)} structure files")
    results_list = []
    for i, fp in enumerate(files):
        print(f"\n{'#'*80}\n# [{i+1}/{len(files)}] {fp}\n{'#'*80}")
        try:
            r = run(fp, output_dir=directory)
            r['filename'] = os.path.basename(fp); r['status'] = 'OK'
            results_list.append(r)
        except Exception as e:
            print(f"ERROR: {e}"); traceback.print_exc()
            results_list.append({'filename': os.path.basename(fp),
                                 'status': f'ERROR: {e}', 'sc_type': 'N/A'})
    if output_csv: _write_csv(results_list, output_csv)
    _print_summary(results_list)
    return results_list

def _write_csv(results, path):
    fields = ['filename','status','sc_type','seekpath_bravais','spacegroup',
              'sg_symbol','point_group','n_symmetry_ops','ibz_volume',
              'centroid_b1','centroid_b2','centroid_b3']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for r in results:
            row = {k: r.get(k, '') for k in fields}
            if 'centroid_frac' in r:
                row['centroid_b1'] = f"{r['centroid_frac'][0]:.6f}"
                row['centroid_b2'] = f"{r['centroid_frac'][1]:.6f}"
                row['centroid_b3'] = f"{r['centroid_frac'][2]:.6f}"
            if 'ibz_volume' in r: row['ibz_volume'] = f"{r['ibz_volume']:.6e}"
            w.writerow(row)
    print(f"\nSaved: {path}")

def _print_summary(results):
    print(f"\n{'='*110}\nSUMMARY\n{'='*110}")
    print(f"{'File':<30s} {'SC Type':>8s} {'SG':>5s} {'Ops':>5s} "
          f"{'Centroid (frac)':<40s} {'Status':>8s}")
    print("-" * 110)
    for r in results:
        cf = r.get('centroid_frac')
        c = f"[{cf[0]:.4f}, {cf[1]:.4f}, {cf[2]:.4f}]" if cf is not None else "N/A"
        print(f"{r.get('filename','?'):<30s} {r.get('sc_type','?'):>8s} "
              f"{r.get('spacegroup','?'):>5} {r.get('n_symmetry_ops','?'):>5} "
              f"{c:<40s} {r.get('status','?'):>8s}")
    print("=" * 110)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python batch_centroid_hybrid.py <directory> [--output summary.csv]"); sys.exit(1)
    output_csv = None; args = sys.argv[1:]
    if '--output' in args:
        idx = args.index('--output'); output_csv = args[idx+1]; args = args[:idx] + args[idx+2:]
    if len(args) == 1 and os.path.isdir(args[0]):
        process_directory(args[0], output_csv=output_csv)
    else:
        results_list = []
        for f in args:
            try:
                r = run(f); r['filename'] = os.path.basename(f); r['status'] = 'OK'; results_list.append(r)
            except Exception as e:
                print(f"ERROR: {f}: {e}")
                results_list.append({'filename':os.path.basename(f), 'status':'ERROR', 'sc_type':'N/A'})
        if output_csv: _write_csv(results_list, output_csv)
        _print_summary(results_list)