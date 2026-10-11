# Expected values of ExpectedEvents(), AverageHR() and RequiredEvents() used in
# tests/testthat, by adaptive quadrature (run: python3 design_ref.py). The last
# section covers the accrual and censoring options of version 0.3.0.
import numpy as np
from scipy import integrate, optimize, special
import reference as R
ctl = R.make_arm(pfs_median=6, orr=0.3, death_prop=0.15, resp_cor=0.4, os_model="idm", os_median=15, kappa=0.6)
trt = R.make_arm(pfs_median=6/0.7, orr=0.45, death_prop=0.15, resp_cor=0.4, os_model="idm", pps_hazard=ctl.gam0, kappa=0.6)
print("trt gam0", trt.gam0, "trt mOS", R.quantile(0.5, trt, "os"), "theta_trt", trt.theta)
A = 24.0
PA = lambda x: np.clip(x / A, 0, 1)
def dens(t, arm, endpoint):
    return R.dens_os(t, arm) if endpoint == "os" else arm.lamP*np.exp(-arm.lamP*t)
def surv(t, arm, endpoint):
    return R.surv_os(t, arm) if endpoint == "os" else np.exp(-arm.lamP*t)
def E(T, n, d, endpoint):
    tot = 0
    for arm, nj, dj in [(ctl, n[0], d[0]), (trt, n[1], d[1])]:
        f = lambda t: dens(t, arm, endpoint)*np.exp(-dj*t)*PA(T-t)
        tot += nj*integrate.quad(f, 0, T, points=[max(T-A,0)] if T > A else None, limit=300, epsabs=1e-10)[0]
    return tot
def ahr(T, n, d, endpoint):
    def hz(t, arm): return dens(t, arm, endpoint)/surv(t, arm, endpoint)
    w = lambda t: sum(nj*dens(t, arm, endpoint)*np.exp(-dj*t)*PA(T-t) for arm, nj, dj in [(ctl,n[0],d[0]),(trt,n[1],d[1])])
    num = integrate.quad(lambda t: np.log(hz(t,trt)/hz(t,ctl))*w(t), 1e-9, T, limit=300, epsabs=1e-10)[0]
    den = integrate.quad(w, 1e-9, T, limit=300, epsabs=1e-10)[0]
    return np.exp(num/den)
def required(n, d, endpoint, alpha=0.025, power=0.8):
    rc = n[0]/sum(n); rt = 1-rc
    z = special.ndtri(1-alpha)+special.ndtri(power)
    D = None; a = ahr(36.0, n, d, endpoint)
    for it in range(50):
        Dstar = z**2/(rc*rt*np.log(a)**2)
        Dint = int(np.ceil(Dstar))
        T = optimize.brentq(lambda T: E(T, n, d, endpoint)-Dint, 1e-3, 500, xtol=1e-10)
        a_new = ahr(T, n, d, endpoint)
        if D == Dint and abs(a_new-a) < 1e-10: break
        D, a = Dint, a_new
    return Dint, T, a, Dstar
n = (350, 350)
for d in [(0.0, 0.0), (0.01, 0.01)]:
    print("dropout", d)
    for T in [12, 24, 36]:
        print("  E_os(%g)=%.6f  E_pfs(%g)=%.6f  AHR_os(%g)=%.6f" % (T, E(T,n,d,"os"), T, E(T,n,d,"pfs"), T, ahr(T,n,d,"os")))
    print("  required OS:", required(n, d, "os"))
    print("  required PFS:", required(n, d, "pfs"))
# ---------------------------------------------------------------------------
# Expected values of the design functions with the accrual and censoring
# options of version 0.3.0 (rates per interval with an open final interval,
# delayed start, piecewise exponential dropout and censoring of PFS).
def accrual_cdf(knots, probs):
    # piecewise-linear distribution function of the entry time
    cum = np.concatenate([[0.0], np.cumsum(probs)])
    return lambda x: np.interp(x, knots, cum, left=0.0, right=1.0)
def cumhaz(times, hazards):
    # piecewise exponential cumulative hazard with pieces starting at times
    times = np.asarray(times, float); hazards = np.asarray(hazards, float)
    cum = np.concatenate([[0.0], np.cumsum(hazards[:-1] * np.diff(times))])
    def H(x):
        i = np.searchsorted(times, x, side="right") - 1
        return cum[i] + hazards[i] * (x - times[i])
    return H
def weight(T, x, arm, nj, Hd, cj, G, endpoint):
    return nj * dens(x, arm, endpoint) * np.exp(-Hd(x) - cj * x) * G(T - x)
def E_gen(T, n, Hds, cs, G, kinks, endpoint):
    tot = 0.0
    pts = sorted(set([p for p in kinks(T) if 0 < p < T]))
    for arm, nj, Hd, cj in [(ctl, n[0], Hds[0], cs[0]), (trt, n[1], Hds[1], cs[1])]:
        tot += integrate.quad(lambda x: weight(T, x, arm, nj, Hd, cj, G, endpoint), 0, T,
                              points=pts or None, limit=500, epsabs=1e-11, epsrel=1e-12)[0]
    return tot
def ahr_gen(T, n, Hds, cs, G, kinks, endpoint):
    pts = sorted(set([p for p in kinks(T) if 0 < p < T]))
    hz = lambda x, arm: dens(x, arm, endpoint) / surv(x, arm, endpoint)
    w = lambda x: sum(weight(T, x, arm, nj, Hd, cj, G, endpoint) for arm, nj, Hd, cj in
                      [(ctl, n[0], Hds[0], cs[0]), (trt, n[1], Hds[1], cs[1])])
    num = integrate.quad(lambda x: np.log(hz(x, trt) / hz(x, ctl)) * w(x), 1e-9, T,
                         points=pts or None, limit=500, epsabs=1e-11)[0]
    den = integrate.quad(w, 1e-9, T, points=pts or None, limit=500, epsabs=1e-11)[0]
    return np.exp(num / den)
def required_gen(n, Hds, cs, G, kinks, endpoint, alpha=0.025, power=0.8, T0=36.0):
    rc = n[0] / sum(n); rt = 1 - rc
    z = special.ndtri(1 - alpha) + special.ndtri(power)
    D = None; a = ahr_gen(T0, n, Hds, cs, G, kinks, endpoint)
    for it in range(50):
        Dstar = z ** 2 / (rc * rt * np.log(a) ** 2)
        Dint = int(np.ceil(Dstar))
        T = optimize.brentq(lambda T: E_gen(T, n, Hds, cs, G, kinks, endpoint) - Dint, 1e-3, 500, xtol=1e-10)
        a_new = ahr_gen(T, n, Hds, cs, G, kinks, endpoint)
        if D == Dint and abs(a_new - a) < 1e-10: break
        D, a = Dint, a_new
    return Dint, T, a, Dstar
print("version 0.3.0 options")
n = (350, 350)
# accrual: 10 per month on [0, 6], then 40 per month until 700 patients (end 22)
a_end = 6 + (700 - 60) / 40
G1 = accrual_cdf([0, 6, a_end], [60 / 700, 640 / 700])
# piecewise exponential dropout: hazard 0.02 on [0, 12), 0.005 after; censoring
# of PFS 0.01 (control) and 0.03 (treatment)
Hd = cumhaz([0, 12], [0.02, 0.005])
k1 = lambda T: [T - 6, T - a_end, 12]
print("  end of accrual", a_end)
for T in [12, 24, 36]:
    print("  E_os(%g)=%.6f  E_pfs_cens(%g)=%.6f  AHR_os(%g)=%.6f" % (
        T, E_gen(T, n, (Hd, Hd), (0, 0), G1, k1, "os"),
        T, E_gen(T, n, (Hd, Hd), (0.01, 0.03), G1, k1, "pfs"),
        T, ahr_gen(T, n, (Hd, Hd), (0, 0), G1, k1, "os")))
print("  AHR_pfs_cens(24)=%.6f" % ahr_gen(24.0, n, (Hd, Hd), (0.01, 0.03), G1, k1, "pfs"))
print("  required OS (power 0.85):", required_gen(n, (Hd, Hd), (0, 0), G1, k1, "os", power=0.85))
print("  required PFS with censoring:", required_gen(n, (Hd, Hd), (0.01, 0.03), G1, k1, "pfs"))
# delayed start: uniform accrual over [3, 27], no dropout
G2 = accrual_cdf([3, 27], [1.0])
k2 = lambda T: [T - 3, T - 27]
zero = cumhaz([0], [0.0])
for T in [12, 24, 36]:
    print("  delayed E_os(%g)=%.6f  check shift %.6f" % (T, E_gen(T, n, (zero, zero), (0, 0), G2, k2, "os"), E(T - 3, n, (0, 0), "os")))
