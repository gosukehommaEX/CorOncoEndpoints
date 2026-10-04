# Expected values of ExpectedEvents(), AverageHR() and RequiredEvents() used in
# tests/testthat, by adaptive quadrature (run: python3 design_ref.py).
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
