"""
Plot 1D profiles from MHD shock tube problems.

Supports overlaying multiple solver runs and a high-res reference solution.
Accepts both AMRVAC .dat files and previously saved .csv files.

Workflow
--------
  # 1. Run simulation (edit par file to select solver/resolution)
  mpirun -np 1 ./amrvac -i briowu.par

  # 2. Extract 1D profile to CSV
  python plot_riemann.py --save-csv results/briowu_tvdlf.csv tvdlf:output/bw0001.dat

  # 3. Repeat for other solvers, then plot comparison
  python plot_riemann.py tvdlf:results/briowu_tvdlf.csv hll:results/briowu_hll.csv

  # With a high-res reference underneath
  python plot_riemann.py --ref results/briowu_ref.csv tvdlf:results/briowu_tvdlf.csv

  # Use solver:path to label each dataset (e.g. hll:path/to/file)
"""

import argparse
import os

import numpy as np
import yt
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
#  Data loading
# ---------------------------------------------------------------------

def load_dat(filename):
    """Load an AMRVAC .dat file and extract a 1D ray along x."""
    ds = yt.load(filename)

    ymid = 0.5 * (ds.domain_left_edge[1] + ds.domain_right_edge[1])
    zmid = 0.5 * (ds.domain_left_edge[2] + ds.domain_right_edge[2])

    ray = ds.ortho_ray("x", (ymid, zmid))
    isort = np.argsort(ray["x"])

    x   = np.array(ray["x"][isort])
    rho = np.array(ray[("amrvac", "rho")][isort])
    m1  = np.array(ray[("amrvac", "m1")][isort])
    m2  = np.array(ray[("amrvac", "m2")][isort])
    m3  = np.array(ray[("amrvac", "m3")][isort])
    e   = np.array(ray[("amrvac", "e")][isort])
    b1  = np.array(ray[("amrvac", "b1")][isort])
    b2  = np.array(ray[("amrvac", "b2")][isort])
    b3  = np.array(ray[("amrvac", "b3")][isort])

    vx = m1 / rho
    vy = m2 / rho
    vz = m3 / rho
    gamma = ds.parameters.get("gamma", 5.0 / 3.0)
    vmag2 = vx**2 + vy**2 + vz**2
    Bmag2 = b1**2 + b2**2 + b3**2
    p = (gamma - 1.0) * (e - 0.5 * rho * vmag2 - 0.5 * Bmag2)

    return {
        "x": x, "rho": rho, "p": p,
        "vx": vx, "vy": vy, "By": b2,
        "time": float(ds.current_time),
        "nx": len(x),
    }


def load_csv(csvfile):
    """Load a profile from CSV."""
    data = np.genfromtxt(csvfile, delimiter=",", names=True)
    prof = {name: data[name] for name in data.dtype.names}
    prof.setdefault("time", 0.0)
    prof.setdefault("nx", len(prof["x"]))
    return prof


def save_csv(profile, csvfile):
    """Save a 1D profile as CSV."""
    os.makedirs(os.path.dirname(csvfile) or ".", exist_ok=True)
    header = "x,rho,p,vx,vy,By"
    arr = np.column_stack([
        profile["x"], profile["rho"], profile["p"],
        profile["vx"], profile["vy"], profile["By"],
    ])
    np.savetxt(csvfile, arr, delimiter=",", header=header, comments="")
    print(f"Saved CSV: {csvfile}")


# ---------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------

PANELS = [
    ("rho", r"$\rho$"),
    ("p",   r"$p$"),
    ("vx",  r"$v_x$"),
    ("vy",  r"$v_y$"),
    ("By",  r"$B_y$"),
]

SOLVER_STYLES = {
    "hll":   {"colour": "C0", "marker": ".", "ms": 1.5, "label": "HLL"},
    "tvdlf": {"colour": "C1", "marker": ".", "ms": 1.5, "label": "TVDLF"},
    "hlld":  {"colour": "C2", "marker": ".", "ms": 1.5, "label": "HLLD"},
}

DEFAULT_STYLE = {"colour": "C3", "marker": ".", "ms": 1.5}


def plot_profiles(profiles, reference=None, savepath="riemann.png"):
    """Multi-panel plot with optional reference and solver comparison."""
    fig, axes = plt.subplots(len(PANELS), 1, figsize=(8, 12), sharex=True)

    p0 = profiles[0]
    solvers = ", ".join(p.get("label", "?") for p in profiles)
    fig.suptitle(
        f"MHD shock tube  [{solvers}]  N = {p0['nx']}",
        fontsize=13,
    )

    # Reference as solid line
    if reference is not None:
        for ax, (key, _) in zip(axes, PANELS):
            if key in reference:
                ax.plot(reference["x"], reference[key], "k-", lw=0.8,
                        alpha=0.5, label="reference", zorder=0)

    # Each solver run
    for prof in profiles:
        label = prof.get("label", "unknown")
        style = SOLVER_STYLES.get(label.lower(), DEFAULT_STYLE)

        for ax, (key, _) in zip(axes, PANELS):
            ax.plot(prof["x"], prof[key],
                    marker=style["marker"], ms=style["ms"],
                    color=style["colour"], ls="none",
                    label=style.get("label", label), zorder=1)

    for ax, (_, ylabel) in zip(axes, PANELS):
        ax.set_ylabel(ylabel, fontsize=12)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("x", fontsize=12)

    # Legend (deduplicated)
    handles, labels = axes[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    axes[0].legend(by_label.values(), by_label.keys(),
                   loc="upper right", fontsize=10)

    fig.tight_layout()
    fig.savefig(savepath, dpi=150)
    print(f"Saved: {savepath}")
    plt.close(fig)


# ---------------------------------------------------------------------
#  CLI
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot 1D MHD shock tube profiles.")
    parser.add_argument("inputs", nargs="+",
                        help="solver:path pairs (.dat or .csv), "
                             "e.g. hll:output/bw0001.dat tvdlf:results/tvdlf.csv")
    parser.add_argument("--ref", help="Reference CSV to plot as solid line")
    parser.add_argument("--save-csv", help="Save first profile as CSV")
    parser.add_argument("-o", "--output", default="riemann.png",
                        help="Output plot filename (default: riemann.png)")
    args = parser.parse_args()

    profiles = []
    for entry in args.inputs:
        if ":" in entry and not entry.startswith("/"):
            label, f = entry.split(":", 1)
        else:
            label, f = "unknown", entry

        print(f"Loading: {f}")
        prof = load_csv(f) if f.endswith(".csv") else load_dat(f)
        prof["label"] = label
        profiles.append(prof)

    if args.save_csv:
        save_csv(profiles[0], args.save_csv)

    reference = None
    if args.ref and os.path.exists(args.ref):
        print(f"Loading reference: {args.ref}")
        reference = load_csv(args.ref)

    if not args.save_csv or args.output != "riemann.png":
        plot_profiles(profiles, reference=reference, savepath=args.output)


if __name__ == "__main__":
    main()
