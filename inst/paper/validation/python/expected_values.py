# Independent computation (Python, loop-based) of the expected values used in
# inst/paper/validation/validate_example_functions.R for a small hand-made
# data set with tied and censored times. Run: python3 expected_values.py
import math
import numpy as np
from scipy import stats

trt = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
os_t = [5.2, 8.1, 3.3, 12.0, 7.5, 2.0, 4.4, 8.1, 6.0, 1.5, 9.9, 3.3, 10.4]
os_e = [1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1]
pfs_t = [3.1, 8.1, 2.0, 6.5, 7.5, 2.0, 2.2, 5.0, 6.0, 1.5, 4.0, 3.3, 6.6]
pfs_e = [1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1]
resp = [1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0]


def logrank(time, event, group):
    """Two-group log-rank: O - E and variance for group == 1."""
    o_e, var = 0.0, 0.0
    for t in sorted(set(ti for ti, ei in zip(time, event) if ei == 1)):
        n = sum(1 for ti in time if ti >= t)
        n1 = sum(1 for ti, gi in zip(time, group) if ti >= t and gi == 1)
        d = sum(1 for ti, ei in zip(time, event) if ti == t and ei == 1)
        d1 = sum(1 for ti, ei, gi in zip(time, event, group) if ti == t and ei == 1 and gi == 1)
        o_e += d1 - d * n1 / n
        if n > 1:
            var += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    return o_e, var, -o_e / math.sqrt(var)


def gehan(ti, ei, tj, ej):
    """Score of patient i against patient j: 1 win, -1 loss, 0 undecided."""
    if ej == 1 and ti > tj:
        return 1
    if ei == 1 and tj > ti:
        return -1
    return 0


def win_stats():
    I = [k for k in range(len(trt)) if trt[k] == 1]
    J = [k for k in range(len(trt)) if trt[k] == 0]
    W = np.zeros((len(I), len(J)))
    L = np.zeros((len(I), len(J)))
    layer = {"os": [0, 0], "pfs": [0, 0], "resp": [0, 0]}
    marg = {"os": [0, 0], "pfs": [0, 0], "resp": [0, 0]}
    for a, i in enumerate(I):
        for b, j in enumerate(J):
            s_os = gehan(os_t[i], os_e[i], os_t[j], os_e[j])
            s_pfs = gehan(pfs_t[i], pfs_e[i], pfs_t[j], pfs_e[j])
            s_r = (resp[i] > resp[j]) - (resp[i] < resp[j])
            for nm, s in (("os", s_os), ("pfs", s_pfs), ("resp", s_r)):
                if s == 1:
                    marg[nm][0] += 1
                elif s == -1:
                    marg[nm][1] += 1
            for nm, s in (("os", s_os), ("pfs", s_pfs), ("resp", s_r)):
                if s != 0:
                    W[a, b] = s == 1
                    L[a, b] = s == -1
                    layer[nm][0 if s == 1 else 1] += 1
                    break
    n1, n0 = len(I), len(J)
    npair = n1 * n0
    pw, pl = W.mean(), L.mean()
    aw, al, bw, bl = W.mean(1), L.mean(1), W.mean(0), L.mean(0)
    cv = lambda x, y: np.mean((x - x.mean()) * (y - y.mean()))
    vw = cv(aw, aw) / n1 + cv(bw, bw) / n0
    vl = cv(al, al) / n1 + cv(bl, bl) / n0
    cwl = cv(aw, al) / n1 + cv(bw, bl) / n0
    nb = pw - pl
    vnb = vw + vl - 2 * cwl
    vlwr = vw / pw ** 2 + vl / pl ** 2 - 2 * cwl / (pw * pl)
    vlwo = vnb * (2 / (1 - nb ** 2)) ** 2
    out = dict(p_win=pw, p_loss=pl, p_tie=1 - pw - pl, wr=pw / pl, wo=(1 + nb) / (1 - nb),
               nb=nb, se_log_wr=math.sqrt(vlwr), se_log_wo=math.sqrt(vlwo), se_nb=math.sqrt(vnb))
    for nm in ("os", "pfs", "resp"):
        out["win_" + nm] = layer[nm][0] / npair
        out["loss_" + nm] = layer[nm][1] / npair
        out["mwin_" + nm] = marg[nm][0] / npair
        out["mloss_" + nm] = marg[nm][1] / npair
    return out


def orr_z():
    x1 = sum(r for r, g in zip(resp, trt) if g == 1)
    x0 = sum(r for r, g in zip(resp, trt) if g == 0)
    n1, n0 = sum(trt), len(trt) - sum(trt)
    p = (x1 + x0) / (n1 + n0)
    return (x1 / n1 - x0 / n0) / math.sqrt(p * (1 - p) * (1 / n1 + 1 / n0))


def independence(w, l):
    t = [1 - a - b for a, b in zip(w, l)]
    pw = w[0] + t[0] * w[1] + t[0] * t[1] * w[2]
    pl = l[0] + t[0] * l[1] + t[0] * t[1] * l[2]
    return pw, pl, t[0] * t[1] * t[2]


def power(pw, pl, pt, n, measure):
    if measure == "wr":
        est, s2 = math.log(pw / pl), 4 * (1 + pt) / (3 * 0.25 * (1 - pt))
    else:
        est = math.log((pw + pt / 2) / (pl + pt / 2))
        s2 = 4 * (1 + pt) * (1 - pt) / (3 * 0.25)
    return stats.norm.cdf(math.sqrt(n / s2) * abs(est) - stats.norm.ppf(0.975))


if __name__ == "__main__":
    print("logrank OS  (o_e, var, z):", ["%.15g" % v for v in logrank(os_t, os_e, trt)])
    print("logrank PFS (o_e, var, z):", ["%.15g" % v for v in logrank(pfs_t, pfs_e, trt)])
    print("orr z: %.15g" % orr_z())
    for k, v in win_stats().items():
        print("win %s: %.15g" % (k, v))
    ind = independence([0.30, 0.20, 0.10], [0.20, 0.15, 0.05])
    print("independence:", ["%.15g" % v for v in ind])
    print("power wr: %.15g" % power(0.45, 0.35, 0.20, 700, "wr"))
    print("power wo: %.15g" % power(0.45, 0.35, 0.20, 700, "wo"))
    s2 = 4 * 1.1 / (3 * 0.25 * 0.9)
    print("Yu and Ganju N: %.15g" % (s2 * (stats.norm.ppf(0.975) + stats.norm.ppf(0.9)) ** 2 / math.log(1.5) ** 2))
