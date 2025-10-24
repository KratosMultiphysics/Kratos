#!/usr/bin/env python3
"""
Scan a folder for HDF5 files matching pattern Adj_*_T_*.h5, read
ResultsData/NodalSolutionStepData/ADJOINT_DISPLACEMENT and
ResultsData/NodalSolutionStepData/ADJOINT_ROTATION from each file,
compute the L2 norm (Frobenius) of either the combined arrays (if rotation
exists) or only displacement (if rotation missing), and plot the norm
versus the T index (the numeric suffix).

This version does NOT use CLI. Set FOLDER and OUT below.
"""
import re
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

H5PATH_DISP = "ResultsData/NodalSolutionStepData/ADJOINT_DISPLACEMENT"
H5PATH_ROT  = "ResultsData/NodalSolutionStepData/ADJOINT_ROTATION"

# --- CONFIG: set these before running the script ---
FOLDER = "/home/ansari/phd/1_Research/1_System_Identification/2_Transient_SI/fresh/cube/system_identification/hdf5_output"
OUT = None                      # <-- set to "out.png" to save, or None to show interactively
# ---------------------------------------------------

def find_files(folder):
    p = Path(folder)
    # only files that start with 'Adj_' and contain the _T_<n> suffix
    files = list(p.glob("Adj_*_T_*.h5")) + list(p.glob("Adj_*_T_*.hdf5"))
    return files

def extract_index(fname):
    m = re.search(r"_T_(\d+)(?:\.h5|\.hdf5)?$", str(fname))
    if m:
        return int(m.group(1))
    m2 = re.search(r"_T_(\d+)", str(fname))
    return int(m2.group(1)) if m2 else None

def read_dataset(h5file, dataset_path):
    with h5py.File(h5file, "r") as f:
        parts = dataset_path.split("/")
        node = f
        try:
            for p in parts:
                node = node[p]
            return np.array(node[()])
        except Exception:
            return None

def l2_norm(arr):
    return float(np.linalg.norm(arr.ravel(), ord=2))

def main(folder, out):
    files = find_files(folder)
    if not files:
        print("No matching h5 files found in", folder, file=sys.stderr)
        return

    entries = []
    used_combined = False
    for f in files:
        idx = extract_index(f)
        if idx is None:
            print("Skipping file with no _T_<n> index:", f, file=sys.stderr)
            continue
        arr_disp = read_dataset(f, H5PATH_DISP)
        arr_rot  = read_dataset(f, H5PATH_ROT)

        # If rotation exists use combined (even if displacement missing)
        if arr_rot is not None:
            used_combined = True
            if arr_disp is None:
                norm = l2_norm(arr_rot)
            else:
                combined = np.concatenate((arr_disp.ravel(), arr_rot.ravel()))
                norm = float(np.linalg.norm(combined, ord=2))
        else:
            # rotation missing -> use displacement only (skip if displacement missing)
            if arr_disp is None:
                print(f"Skipping {f}: neither ADJOINT_ROTATION nor ADJOINT_DISPLACEMENT found", file=sys.stderr)
                continue
            norm = l2_norm(arr_disp)

        entries.append((idx, norm, str(f)))

    if not entries:
        print("No valid entries read.", file=sys.stderr)
        return

    entries.sort(key=lambda x: x[0])
    indices = [e[0] for e in entries]
    norms = np.array([e[1] for e in entries], dtype=float)

    # avoid zeros for log scale
    eps = max(np.finfo(float).tiny, 1e-30)
    norms = np.maximum(norms, eps)

    plt.figure(figsize=(8,4.5))
    #plt.yscale("log")
    plt.plot(indices, norms, marker="o", linestyle="-")
    plt.xlabel("T index")
    if used_combined:
        plt.ylabel("L2 norm of (ADJOINT_DISPLACEMENT + ADJOINT_ROTATION) (log scale)")
        plt.title("Combined ADJOINT variables L2 norm vs T index")
    else:
        plt.ylabel("L2 norm of ADJOINT_DISPLACEMENT (log scale)")
        plt.title("ADJOINT_DISPLACEMENT L2 norm vs T index")
    plt.grid(True, which="both", ls=":")

    if out:
        plt.tight_layout()
        plt.savefig(out, dpi=200)
        print("Saved plot to", out)
    else:
        plt.show()

if __name__ == "__main__":
    if not FOLDER:
        print("Please set FOLDER at the top of the script to the folder containing the H5 files.", file=sys.stderr)
        sys.exit(1)
    main(FOLDER, OUT)