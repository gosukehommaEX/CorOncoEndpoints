# Compares the C++ transformation core (src/transform_core.h) with reference.py.
# Run in this folder:
#   g++ -O2 -std=c++11 -o transform_core_harness transform_core_harness.cpp
#   python3 compare_transform_core.py
# The comparison writes temporary text files in the current folder.
import sys, subprocess, numpy as np
sys.path.insert(0, ".")
import reference as R
rng = np.random.default_rng(5)
scen = {
 "A": dict(pfs_median=6, orr=0.3, death_prop=0.15, resp_cor=0.4, os_model="idm", os_median=15, kappa=0.6),
 "B": dict(pfs_median=6, orr=0.3, death_prop=0.15, resp_pfs_median=11, os_model="idm", pps_median=8, kappa=1.0, tau=1.5),
 "C": dict(pfs_median=6, orr=0.3, death_prop=0.15, resp_cor=0.4, os_model="expexp", os_median=15, tau=1.5, ttr_median=2.5),
 "D": dict(pfs_median=6, orr=0.3, death_prop=0.4, resp_cor=0.4, os_model="expexp", os_median=15),
}
n = 200000
for k, kw in scen.items():
    A = R.make_arm(**kw)
    z1 = rng.standard_normal(n); z1[:6] = [-9, -5, 0, 5, 9, 30]
    w = rng.standard_normal(n); ut = rng.random(n); ep = rng.exponential(1.0, n); uttr = rng.random(n)
    np.savetxt("inputs.txt", np.column_stack([z1, w, ut, ep, uttr]), fmt="%.17g")
    if A.os_model == "idm":
        grid = np.array([0.0, 0.0]); gh = 1.0; lamO = 1.0; cdec = 0.0; g0, g1 = A.gam0, A.gam1
    else:
        grid = A.gG; gh = A.gt[1] - A.gt[0]; lamO = A.lamO; cdec = A.cdec; g0 = g1 = 1.0
    np.savetxt("grid.txt", grid, fmt="%.17g")
    ttr_mode = 0 if A.ttr_rate is None else 1
    with open("params.txt", "w") as f:
        f.write(" ".join("%.17g" % v for v in [A.lamP, A.theta, A.c, A.ztau if np.isfinite(A.ztau) else -1e300, A.tau,
                 0 if A.os_model == "idm" else 1, A.pi, g0, g1, lamO, cdec, gh, ttr_mode, A.ttr_rate or 1.0]))
    subprocess.run(["./transform_core_harness", "params.txt", "grid.txt", "inputs.txt", "out.txt", "acc.txt"], check=True)
    out = np.loadtxt("out.txt")
    pfs, os_, resp, prog, ttr = R.transform(A, z1, w, ut, ep, uttr)
    # numpy uses the grid with uniform spacing gt; C++ uses (i)*h -- compare
    d_pfs = np.max(np.abs(out[:, 0] - pfs) / np.maximum(1, pfs))
    d_os = np.max(np.abs(out[:, 1] - os_) / np.maximum(1, os_))
    print(k, "max rel diff pfs %.2e os %.2e" % (d_pfs, d_os),
          "prog mismatches", int(np.sum(out[:, 2] != prog)), "resp mismatches", int(np.sum(out[:, 3] != resp)),
          "ttr max diff %.2e" % (np.nanmax(np.abs(out[:, 4] - ttr)) if ttr_mode else 0),
          "ttr NaN pattern equal", bool(np.all(np.isnan(out[:, 4]) == np.isnan(ttr))))
    print("   extreme z1 rows pfs:", out[:6, 0], "python:", pfs[:6])
acc = np.loadtxt("acc.txt")
print("accrual inverse:", acc[:, 1].tolist())
