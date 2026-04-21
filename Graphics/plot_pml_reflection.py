"""
plot_pml_reflection.py
Читает CSV-файлы, записанные C++-программой, и строит графики:
  1) R(f) для трёх профилей PML при фиксированной ширине
  2) R_max(width) для трёх профилей PML
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ─── Параметры ──────────────────────────────────────────────────
PROFILE_NAMES = {0: "constant (m=0)", 2: "quadratic (m=2)", 3: "cubic (m=3)"}
PROFILE_COLORS = {0: "#e74c3c", 2: "#2ecc71", 3: "#3498db"}
FIXED_WIDTH = 20        # ширина PML для графика спектра
DATA_DIR = "."          # директория с CSV файлами (где лежит этот скрипт)

# ─── Вспомогательные функции ─────────────────────────────────────
def load_csv(path):
    """Загружает CSV с заголовком, возвращает dict{col: array}."""
    import csv
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        cols = {k: [] for k in reader.fieldnames}
        for row in reader:
            for k, v in row.items():
                cols[k].append(float(v))
    return {k: np.array(v) for k, v in cols.items()}

def file_path(*parts):
    return os.path.join(DATA_DIR, *parts)

# ─── График 1: R(f) для фиксированной ширины, три профиля ────────
def plot_spectra_vs_profile(fixed_width=FIXED_WIDTH):
    fig, ax = plt.subplots(figsize=(8, 5))

    for m, name in PROFILE_NAMES.items():
        fname = file_path(f"C:\\Users\\6anna\\CLionProjects\\Sem6Laba3_FDTD\\cmake-build-debug\\reflection_profile_{m}_width_{fixed_width}.csv")
        if not os.path.exists(fname):
            print(f"  [warn] {fname} not found, skipping")
            continue
        d = load_csv(fname)
        freq = d["freq"]
        R    = d["R"]
        ax.semilogy(freq, R, color=PROFILE_COLORS[m], lw=2, label=name)

    ax.set_xlabel("Frequency $f$ (normalised units)")
    ax.set_ylabel("Reflection coefficient $R(f)$")
    ax.set_title(f"PML Reflection Spectrum  (width = {fixed_width} cells)")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext())
    plt.tight_layout()
    out = file_path("plot_R_vs_freq.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

# ─── График 2: R_max(width) для трёх профилей ────────────────────
def plot_maxR_vs_width():
    fname = file_path("reflection_maxR_vs_width.csv")
    if not os.path.exists(fname):
        print(f"  [warn] {fname} not found, skipping")
        return

    d = load_csv(fname)
    widths = d["width"].astype(int)

    fig, ax = plt.subplots(figsize=(8, 5))

    for m, name in PROFILE_NAMES.items():
        col = f"R_profile{m}"
        if col not in d:
            print(f"  [warn] column {col} not found")
            continue
        ax.semilogy(widths, d[col], "o-", color=PROFILE_COLORS[m],
                    lw=2, ms=6, label=name)

    ax.set_xlabel("PML width (cells)")
    ax.set_ylabel("Max reflection coefficient $R_{\\max}$")
    ax.set_title("Maximum PML Reflection vs Width")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext())
    plt.tight_layout()
    out = file_path("plot_maxR_vs_width.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

# ─── Бонус: R(f) для кубического профиля при разных ширинах ──────
def plot_spectra_vs_width(profile=3):
    widths_to_check = [5, 10, 20, 30]
    fig, ax = plt.subplots(figsize=(8, 5))

    cmap = plt.get_cmap("plasma")
    colors = [cmap(i / (len(widths_to_check) - 1)) for i in range(len(widths_to_check))]

    for w, c in zip(widths_to_check, colors):
        fname = file_path(f"C:\\Users\\6anna\\CLionProjects\\Sem6Laba3_FDTD\\cmake-build-debug\\reflection_cubic_width_{w}.csv")
        if not os.path.exists(fname):
            print(f"  [warn] {fname} not found")
            continue
        d = load_csv(fname)
        ax.semilogy(d["freq"], d["R"], color=c, lw=2, label=f"width={w}")

    ax.set_xlabel("Frequency $f$ (normalised units)")
    ax.set_ylabel("Reflection coefficient $R(f)$")
    ax.set_title(f"PML Reflection: cubic profile, varying width")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    out = file_path("plot_R_vs_freq_cubic.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

# ─── Точка входа ─────────────────────────────────────────────────
if __name__ == "__main__":
    print("=== Plotting PML reflection results ===")
    plot_spectra_vs_profile(FIXED_WIDTH)
    plot_maxR_vs_width()
    plot_spectra_vs_width(3)
    print("Done.")