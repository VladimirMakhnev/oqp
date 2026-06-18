"""Render an attachment/detachment density figure for formaldehyde S1 (n->pi*).

Projects the densities onto the molecular yz-plane (summing over x) and draws
filled contours side by side, with nuclei marked.  Output: PNG for the report.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from oqp.pyoqp import Runner
from oqp.analysis import MRSFExcitedStates, AOBasis, make_box_grid, attachment_detachment
from oqp.export.cubegen import _density_values

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "_runs", "figs")
os.makedirs(OUT, exist_ok=True)

r = Runner(project="ch2o_mrsf", input_file=os.path.join(HERE, "inputs", "ch2o_mrsf.inp"),
           log=os.path.join(HERE, "_runs", "ch2o_fig.log"), silent=1, usempi=False)
r.run()
st = MRSFExcitedStates(r.mol)
ao = AOBasis(r.mol)
ad = attachment_detachment(st, 1, ref=0)

lo, n, dvec, pts = make_box_grid(ao.coords, padding=4.0, spacing=0.12)
nx, ny, nz = n
det = _density_values(ao, ad["D_ao"], pts).reshape(nx, ny, nz).sum(axis=0) * dvec[0]
att = _density_values(ao, ad["A_ao"], pts).reshape(nx, ny, nz).sum(axis=0) * dvec[0]
yax = lo[1] + dvec[1] * np.arange(ny)
zax = lo[2] + dvec[2] * np.arange(nz)
Z, Y = np.meshgrid(zax, yax)
coords = ao.coords

fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), constrained_layout=True)
for axc, dat, title in ((axes[0], det, "Detachment (hole)"),
                        (axes[1], att, "Attachment (particle)")):
    vmax = float(np.abs(dat).max())
    cf = axc.contourf(Z, Y, dat, levels=20, cmap="viridis", vmin=0, vmax=vmax)
    axc.scatter(coords[:, 2], coords[:, 1], c="white", edgecolors="k", s=80, zorder=5)
    for c in coords:
        axc.annotate("", (c[2], c[1]))
    axc.set_title(title)
    axc.set_xlabel("z (bohr)"); axc.set_ylabel("y (bohr)")
    axc.set_aspect("equal")
    fig.colorbar(cf, ax=axc, shrink=0.8)
fig.suptitle(r"Formaldehyde S$_1$ (n$\to\pi^*$) unrelaxed attachment/detachment densities "
             r"(MRSF-BHHLYP/6-31G$^*$)")
path = os.path.join(OUT, "ch2o_S1_attach_detach.png")
fig.savefig(path, dpi=150)
print("wrote", path)
