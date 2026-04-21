"""
Читает silica_reflection.csv и строит R(λ): FDTD vs Теория Френеля.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = "."

def load_csv(path):
    import csv
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        cols = {k: [] for k in reader.fieldnames}
        for row in reader:
            for k, v in row.items():
                cols[k].append(float(v))
    return {k: np.array(v) for k, v in cols.items()}

def plot_silica():
    fname = os.path.join(DATA_DIR, "C:\\Users\\6anna\\CLionProjects\\Sem6Laba3_FDTD\\cmake-build-debug\\silica_reflection.csv")
    if not os.path.exists(fname):
        print(f"[warn] {fname} not found")
        return

    d = load_csv(fname)
    lam  = d["lambda_nm"]
    Rfdtd  = d["R_fdtd"]
    Rtheory = d["R_theory"]

    # Теоретическое значение Френеля (константа)
    n1, n2 = 1.0, 1.45
    R_fresnel_val = ((n1 - n2) / (n1 + n2))**2
    print(f"Theoretical Fresnel R = {R_fresnel_val:.6f}  (~{R_fresnel_val*100:.3f}%)")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # ─── График 1: R vs λ ─────────────────────────────────────
    ax1.plot(lam, Rfdtd,   color="#3498db", lw=2, label="FDTD")
    ax1.axhline(R_fresnel_val, color="#e74c3c", lw=1.5,
                ls="--", label=f"Fresnel theory  R={R_fresnel_val:.4f}")
    ax1.set_xlabel(r"Wavelength $\lambda$ (nm)")
    ax1.set_ylabel("Reflection coefficient $R$")
    ax1.set_title("Reflection from SiO\u2082 half-space")
    ax1.set_xlim(lam.min(), lam.max())
    ax1.set_ylim(0, max(Rfdtd.max(), R_fresnel_val) * 1.3)
    ax1.legend()
    ax1.grid(True, ls="--", alpha=0.5)

    # ─── График 2: отклонение от теории ──────────────────────
    diff = Rfdtd - R_fresnel_val
    ax2.plot(lam, diff, color="#8e44ad", lw=2)
    ax2.axhline(0, color="gray", lw=1, ls="--")
    ax2.set_xlabel(r"Wavelength $\lambda$ (nm)")
    ax2.set_ylabel(r"$R_{FDTD} - R_{theory}$")
    ax2.set_title("Deviation from Fresnel theory")
    ax2.set_xlim(lam.min(), lam.max())
    ax2.grid(True, ls="--", alpha=0.5)

    plt.tight_layout()
    out = os.path.join(DATA_DIR, "plot_silica_reflection.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

if __name__ == "__main__":
    print("=== Plotting SiO2 reflection ===")
    plot_silica()
    print("Done.")
