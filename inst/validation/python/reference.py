# Independent Python reference implementation of the CorOncoEndpoints model.
# Used only to compute expected values for the R tests (tests/testthat) and to
# check the C++ transform core. It is written separately from the R code.
#
# Model for one arm
#   (Z1, Z2) standard bivariate normal with correlation theta
#   PFS = -log(Phi(-Z1)) / lamP                                   (exactly Exp(lamP))
#   R   = 1{PFS > tau} 1{Z2 > c},  c chosen so that P(R = 1) = p   (landmark when tau > 0)
#   Death before progression:
#     idm    : with probability pi, independent of PFS
#     expexp : with probability (lamO / lamP) exp(-cdec * PFS),  cdec = lamO / pi - lamP
#   Post-progression survival D:
#     idm    : Exp(gam_R), gam_1 = kappa * gam_0 (clock reset)
#     expexp : Markov hazard h12(t) = lamO (1 - exp(-k t)) / (1 - exp(-b t)),
#              k = lamP + cdec - lamO, b = lamP - lamO
#   OS = PFS + Delta * D
#   TTR (responders): tau + Exp(rho) truncated to [0, PFS - tau), rho = log 2 / (ttr_median - tau)

import numpy as np
from scipy import integrate, optimize, special, stats

LOG2 = np.log(2.0)


def bvn_upper(a, c, theta):
    # P(Z1 > a, Z2 > c) via Owen's T (exact), handling infinite a
    if a == -np.inf:
        return special.ndtr(-c)
    if a == np.inf:
        return 0.0
    # P(Z1 > a, Z2 > c) = Phi2(-a, -c; theta)
    h, k = -a, -c
    if abs(h) < 1e-14:
        h = 1e-14
    if abs(k) < 1e-14:
        k = 1e-14
    r = np.sqrt(1 - theta ** 2)
    ah, ak = (k - theta * h) / (h * r), (h - theta * k) / (k * r)
    beta = 0.0 if h * k > 0 else 0.5
    return 0.5 * special.ndtr(h) + 0.5 * special.ndtr(k) - special.owens_t(h, ah) - special.owens_t(k, ak) - beta


class Arm:
    pass


def z_of_t(t, lamP):
    return special.ndtri(-np.expm1(-lamP * np.asarray(t, dtype=float)))


def resp_threshold(theta, p, ztau):
    if ztau == -np.inf:
        return special.ndtri(1 - p)
    return optimize.brentq(lambda c: bvn_upper(ztau, c, theta) - p, -12, 12, xtol=1e-14)


def q_of_z(z, theta, c, ztau):
    s = np.sqrt(1 - theta ** 2)
    return np.where(z > ztau, special.ndtr((theta * z - c) / s), 0.0)


def corr_pfs_r_std(theta, p, ztau):
    c = resp_threshold(theta, p, ztau)
    lo = max(ztau, -12.0)
    f = lambda z: (-special.log_ndtr(-z) - 1.0) * float(q_of_z(z, theta, c, ztau)) * stats.norm.pdf(z)
    pts = [c / theta] if abs(theta) > 1e-8 and lo < c / theta < 12 else None
    v = integrate.quad(f, lo, 12.0, points=pts, limit=500, epsabs=1e-13, epsrel=1e-12)[0]
    return v / np.sqrt(p * (1 - p)), c


def surv_pfs_resp(t, A):
    # responders' PFS survival
    z = max(float(z_of_t(t, A.lamP)), A.ztau)
    return bvn_upper(z, A.c, A.theta) / A.p


def make_arm(pfs_median, orr, death_prop, resp_cor=None, resp_pfs_median=None, os_model="idm",
             os_median=None, pps_median=None, pps_hazard=None, kappa=1.0, tau=0.0, ttr_median=None):
    A = Arm()
    A.lamP = LOG2 / pfs_median
    A.p = orr
    A.pi = death_prop
    A.tau = tau
    A.ztau = -np.inf if tau == 0 else float(z_of_t(tau, A.lamP))
    A.os_model = os_model
    A.kappa = kappa
    A.ttr_rate = None if ttr_median is None else LOG2 / (ttr_median - tau)
    if resp_cor is not None:
        A.theta = optimize.brentq(lambda th: corr_pfs_r_std(th, orr, A.ztau)[0] - resp_cor, -0.9999, 0.9999, xtol=1e-13)
    else:
        def g(th):
            A.theta = th
            A.c = resp_threshold(th, orr, A.ztau)
            return surv_pfs_resp(resp_pfs_median, A) - 0.5
        A.theta = optimize.brentq(g, -0.9999, 0.9999, xtol=1e-13)
    A.c = resp_threshold(A.theta, orr, A.ztau)
    if os_model == "idm":
        if os_median is not None:
            A.gam0 = 1.0
            A.gam0 = np.exp(optimize.brentq(lambda lg: surv_os(os_median, A, g0=np.exp(lg)) - 0.5, np.log(1e-8), np.log(1e3), xtol=1e-14))
        elif pps_median is not None:
            A.gam0 = LOG2 / pps_median
        else:
            A.gam0 = pps_hazard
        A.gam1 = kappa * A.gam0
    else:
        A.lamO = LOG2 / os_median
        A.cdec = A.lamO / death_prop - A.lamP
        A.k = A.lamP + A.cdec - A.lamO
        A.b = A.lamP - A.lamO
        tmax = 40.0 / A.lamO
        A.gt = np.linspace(0, tmax, 200001)
        h = h12(A.gt, A)
        A.gG = np.concatenate([[0.0], np.cumsum(0.5 * (h[1:] + h[:-1]) * np.diff(A.gt))])
    return A


def h12(t, A):
    t = np.asarray(t, dtype=float)
    lim0 = A.lamO * A.k / A.b
    with np.errstate(invalid="ignore", divide="ignore"):
        v = A.lamO * (-np.expm1(-A.k * t)) / (-np.expm1(-A.b * t))
    return np.where(t < 1e-12, lim0, v)


def Gfun(t, A):
    t = np.asarray(t, dtype=float)
    out = np.interp(t, A.gt, A.gG)
    over = t > A.gt[-1]
    return np.where(over, A.gG[-1] + A.lamO * (t - A.gt[-1]), out)


def q_of_s(s, A):
    return q_of_z(z_of_t(s, A.lamP), A.theta, A.c, A.ztau)


def surv_os(t, A, g0=None, response="all"):
    lamP, pi = A.lamP, A.pi
    if response == "all":
        S_p, wfun, norm = np.exp(-lamP * t), None, 1.0
    elif response == "resp":
        S_p, norm = surv_pfs_resp(t, A), A.p
    else:
        S_p, norm = (np.exp(-lamP * t) - A.p * surv_pfs_resp(t, A)) / (1 - A.p), 1 - A.p
    if t <= 0:
        return 1.0
    pts = [A.tau] if 0 < A.tau < t else None
    if A.os_model == "idm":
        g0 = A.gam0 if g0 is None else g0
        g1 = A.kappa * g0
        def f(s):
            q = float(q_of_s(s, A))
            fP = lamP * np.exp(-lamP * s)
            if response == "all":
                return fP * ((1 - q) * np.exp(-g0 * (t - s)) + q * np.exp(-g1 * (t - s)))
            if response == "resp":
                return fP * q * np.exp(-g1 * (t - s)) / norm
            return fP * (1 - q) * np.exp(-g0 * (t - s)) / norm
        return S_p + (1 - pi) * integrate.quad(f, 0, t, points=pts, limit=500, epsabs=1e-13, epsrel=1e-12)[0]
    # expexp
    if response == "all":
        return np.exp(-A.lamO * t)
    Gt = float(Gfun(t, A))
    def f(s):
        q = float(q_of_s(s, A))
        w = q if response == "resp" else 1 - q
        fP = lamP * np.exp(-lamP * s)
        pis = (A.lamO / lamP) * np.exp(-A.cdec * s)
        return fP * w / norm * (1 - pis) * np.exp(-(Gt - float(Gfun(s, A))))
    return S_p + integrate.quad(f, 0, t, points=pts, limit=500, epsabs=1e-12, epsrel=1e-11)[0]


def dens_os(t, A):
    lamP, pi = A.lamP, A.pi
    fP_t = lamP * np.exp(-lamP * t)
    if A.os_model == "idm":
        g0, g1 = A.gam0, A.gam1
        def f(s):
            q = float(q_of_s(s, A))
            return lamP * np.exp(-lamP * s) * ((1 - q) * g0 * np.exp(-g0 * (t - s)) + q * g1 * np.exp(-g1 * (t - s)))
        pts = [A.tau] if 0 < A.tau < t else None
        v = integrate.quad(f, 0, t, points=pts, limit=500, epsabs=1e-13, epsrel=1e-12)[0] if t > 0 else 0.0
        return pi * fP_t + (1 - pi) * v
    return A.lamO * np.exp(-A.lamO * t)


def quantile(prob, A, endpoint="os", response="all"):
    if endpoint == "pfs":
        if response == "all":
            S = lambda t: np.exp(-A.lamP * t)
        elif response == "resp":
            S = lambda t: surv_pfs_resp(t, A)
        else:
            S = lambda t: (np.exp(-A.lamP * t) - A.p * surv_pfs_resp(t, A)) / (1 - A.p)
    else:
        S = lambda t: surv_os(t, A, response=response)
    return optimize.brentq(lambda t: S(t) - (1 - prob), 1e-9, 1e4, xtol=1e-12)


def correlations(A):
    p, lamP = A.p, A.lamP
    rho_pr, _ = corr_pfs_r_std(A.theta, p, A.ztau)
    cov_pr = rho_pr * np.sqrt(p * (1 - p)) / lamP
    var_p = 1 / lamP ** 2
    if A.os_model == "idm":
        m0, m1 = (1 - A.pi) / A.gam0, (1 - A.pi) / A.gam1
        s0, s1 = 2 * (1 - A.pi) / A.gam0 ** 2, 2 * (1 - A.pi) / A.gam1 ** 2
        d = m1 - m0
        var_dd = p * s1 + (1 - p) * s0 - (p * m1 + (1 - p) * m0) ** 2
        cov_or = cov_pr + p * (1 - p) * d
        cov_po = var_p + d * cov_pr
        var_o = var_p + var_dd + 2 * d * cov_pr
    else:
        # m(s) = E[D | progression at s] = exp(G(s)) int_s^inf exp(-G(t)) dt, tabulated then interpolated
        t, G = A.gt, A.gG
        eG = np.exp(-G)
        I = np.concatenate([np.cumsum((0.5 * (eG[1:] + eG[:-1]) * np.diff(t))[::-1])[::-1], [0.0]])
        I = I + eG[-1] / A.lamO
        m = np.exp(G) * I
        pis = (A.lamO / lamP) * np.exp(-A.cdec * t)
        gtab = (1 - pis) * m
        g = lambda x: np.interp(x, t, gtab)
        fP = lambda x: lamP * np.exp(-lamP * x)
        top = t[-1]
        pts = [A.tau] if A.tau > 0 else None
        qd = lambda f: integrate.quad(f, 0, top, points=pts, limit=1000, epsabs=1e-12, epsrel=1e-11)[0]
        e_pdd = qd(lambda x: x * fP(x) * g(x))
        e_dd = qd(lambda x: fP(x) * g(x))
        cov_po = var_p + e_pdd - (1 / lamP) * e_dd
        cov_or = cov_pr + qd(lambda x: g(x) * (float(q_of_s(x, A)) - p) * fP(x))
        var_o = 1 / A.lamO ** 2
    var_r = p * (1 - p)
    r_pr = cov_pr / np.sqrt(var_p * var_r)
    r_or = cov_or / np.sqrt(var_o * var_r)
    r_po = cov_po / np.sqrt(var_p * var_o)
    partial = (r_or - r_pr * r_po) / np.sqrt((1 - r_pr ** 2) * (1 - r_po ** 2))
    return dict(cor_pfs_resp=r_pr, cor_os_resp=r_or, cor_pfs_os=r_po, partial_os_resp=partial)


# ---------------------------------------------------------------------------
def simulate(A, n, rng):
    z1 = rng.standard_normal(n)
    w = rng.standard_normal(n)
    u_type = rng.random(n)
    e_pps = rng.exponential(1.0, n)
    u_ttr = rng.random(n)
    return transform(A, z1, w, u_type, e_pps, u_ttr)


def transform(A, z1, w, u_type, e_pps, u_ttr):
    s = np.sqrt(1 - A.theta ** 2)
    pfs = -special.log_ndtr(-z1) / A.lamP
    z2 = A.theta * z1 + s * w
    resp = ((z1 > A.ztau) & (z2 > A.c)).astype(int)
    if A.os_model == "idm":
        death = u_type < A.pi
        rate = np.where(resp == 1, A.gam1, A.gam0)
        os_ = np.where(death, pfs, pfs + e_pps / rate)
    else:
        death = u_type < (A.lamO / A.lamP) * np.exp(-A.cdec * pfs)
        target = Gfun(pfs, A) + e_pps
        tt = np.interp(target, A.gG, A.gt)
        over = target > A.gG[-1]
        tt = np.where(over, A.gt[-1] + (target - A.gG[-1]) / A.lamO, tt)
        os_ = np.where(death, pfs, tt)
    ttr = np.full(len(z1), np.nan)
    if A.ttr_rate is not None:
        rr = resp == 1
        width = pfs[rr] - A.tau
        ttr[rr] = A.tau - np.log1p(-u_ttr[rr] * (-np.expm1(-A.ttr_rate * width))) / A.ttr_rate
    return pfs, os_, resp, (~death).astype(int), ttr
