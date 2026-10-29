"""Generate an informative 2-D sketch of the periodically replicated MD system."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import numpy as np


OUTPUT_DIR = Path(__file__).resolve().parent / "final_results"
OUTPUT_DIR.mkdir(exist_ok=True)

L = 10.0
SIGMA = 1.0
R_C = 3.0
K = 1600.0
MASS = 1.0
PARTICLE_RADIUS = SIGMA / 2
RC_COLOR = "#555555"


def draw_spring(ax, start, end, coils=8, amplitude=0.11, color="black"):
    """Draw a spring from start to end with a short straight lead at each end."""
    start = np.asarray(start, dtype=float)
    end = np.asarray(end, dtype=float)
    vector = end - start
    length = np.linalg.norm(vector)
    direction = vector / length
    normal = np.array([-direction[1], direction[0]])
    lead = min(0.28, length * 0.18)
    first = start + lead * direction
    last = end - lead * direction
    phase = np.linspace(0.0, 2.0 * np.pi * coils, 180)
    axis = np.linspace(0.0, 1.0, phase.size)
    points = first[:, None] + (last - first)[:, None] * axis
    points += amplitude * np.sin(phase)[None, :] * normal[:, None]
    ax.plot(
        [start[0], first[0]],
        [start[1], first[1]],
        color=color,
        lw=1.8,
        solid_capstyle="round",
        zorder=5,
    )
    ax.plot(points[0], points[1], color=color, lw=1.8, zorder=5)
    ax.plot(
        [last[0], end[0]],
        [last[1], end[1]],
        color=color,
        lw=1.8,
        solid_capstyle="round",
        zorder=5,
    )


def add_particles(ax, positions, color, edgecolor, alpha=1.0, zorder=3):
    for x, y, radius in positions:
        ax.add_patch(
            Circle(
                (x, y),
                radius,
                facecolor=color,
                edgecolor=edgecolor,
                linewidth=1.0,
                alpha=alpha,
                zorder=zorder,
            )
        )


def main():
    fig, ax = plt.subplots(figsize=(10.4, 7.8), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_aspect("equal")

    # The anchor is the center of the image. Only half of each neighboring
    # periodic box is visible.
    ax.set_xlim(-7.7, 7.7)
    ax.set_ylim(-7.7, 7.7)
    for dx in (-L, 0.0, L):
        for dy in (-L, 0.0, L):
            is_central = dx == 0.0 and dy == 0.0
            ax.add_patch(
                Rectangle(
                    (-L / 2 + dx, -L / 2 + dy),
                    L,
                    L,
                    fill=False,
                    edgecolor="black" if is_central else "#888888",
                    linewidth=2.2 if is_central else 1.0,
                    linestyle="-" if is_central else "--",
                    alpha=1.0 if is_central else 0.8,
                    zorder=1,
                )
            )

    # Sparse environment: enough particles to communicate the fluid background.
    environment = np.array(
        [
            [-3.8, -3.6], [-2.0, -3.2], [-0.1, -3.8], [2.0, -3.4],
            [3.8, -3.6], [-3.4, -1.7], [-1.6, -1.0], [2.7, -1.5],
            [3.8, 0.0], [-3.5, 1.4], [-1.2, 2.1], [1.8, 2.3],
            [3.6, 3.1], [-3.7, 3.8], [-0.4, 3.5], [2.2, 3.8],
        ]
    )
    base_positions = [(x, y, PARTICLE_RADIUS) for x, y in environment]
    for dx in (-L, 0.0, L):
        for dy in (-L, 0.0, L):
            if dx == 0.0 and dy == 0.0:
                continue
            shifted = [(x + dx, y + dy, r) for x, y, r in base_positions]
            add_particles(ax, shifted, "#d0d0d0", "#a0a0a0", alpha=0.6, zorder=2)
    add_particles(ax, base_positions, "#b8b8b8", "#4a4a4a", zorder=3)

    # The anchored position and a slightly displaced elastic particle.
    anchor = np.array([0.0, 0.0])
    particle = np.array([1.25, 0.9])
    draw_spring(ax, anchor, particle)
    ax.add_patch(Circle(tuple(anchor), 0.10, color="black", zorder=6))
    ax.plot(
        [anchor[0], particle[0]],
        [anchor[1], particle[1]],
        color="black",
        lw=1.0,
        linestyle=":",
        alpha=0.8,
        zorder=4,
    )
    ax.add_patch(
        Circle(
            tuple(particle),
            SIGMA / 2,
            facecolor="#969696",
            edgecolor="#202020",
            linewidth=1.6,
            zorder=7,
        )
    )

    # Cutoff is deliberately drawn only around the elastic particle.
    ax.add_patch(
        Circle(
            tuple(particle),
            R_C,
            fill=False,
            edgecolor=RC_COLOR,
            linewidth=2.2,
            linestyle=(0, (6, 4)),
            zorder=4,
        )
    )
    # Radius marker for the cutoff: the label sits directly on the segment.
    rc_angle = np.deg2rad(35.0)
    rc_direction = np.array([np.cos(rc_angle), np.sin(rc_angle)])
    rc_end = particle + R_C * rc_direction
    ax.plot(
        [particle[0], rc_end[0]],
        [particle[1], rc_end[1]],
        color="black",
        lw=1.0,
        zorder=8,
    )
    rc_midpoint = particle + 0.58 * R_C * rc_direction
    ax.text(
        rc_midpoint[0] + 0.15,
        rc_midpoint[1] + 0.35,
        r"$r_c$",
        ha="center",
        va="center",
        color="black",
        fontsize=14,
        zorder=9,
    )

    dimension_y = -5.85
    dimension_x0 = -L / 2
    dimension_x1 = L / 2
    ax.plot(
        [dimension_x0, dimension_x1],
        [dimension_y, dimension_y],
        color="black",
        lw=1.4,
        zorder=8,
    )
    ax.plot(
        [dimension_x0, dimension_x0],
        [dimension_y - 0.16, dimension_y + 0.16],
        color="black",
        lw=1.4,
        zorder=8,
    )
    ax.plot(
        [dimension_x1, dimension_x1],
        [dimension_y - 0.16, dimension_y + 0.16],
        color="black",
        lw=1.4,
        zorder=8,
    )
    ax.text(
        x=0.0,
        y=dimension_y - 0.38,
        s=r"$L=\sqrt{V}$",
        fontsize=15,
        color="black",
        ha="center",
        va="top",
        zorder=9,
    )

    ax.axis("off")

    png_path = OUTPUT_DIR / "md_system_scheme.png"
    pdf_path = OUTPUT_DIR / "md_system_scheme.pdf"
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"Saved {png_path}")
    print(f"Saved {pdf_path}")
    plt.show()


if __name__ == "__main__":
    main()
