"""Survival analysis primitives, implemented here because lifelines / scikit-survival are not installed.

Everything below is standard, but a silent bug in any of it would poison every number the survival layer
reports, so `python coxph.py` runs known-answer checks: a simulated cohort with a known coefficient, a
C-index whose value is analytically obvious, and a Kaplan-Meier / log-rank pair checked against a
hand-worked example.

The one thing that matters more than the maths: **`overallSurvival` for a living patient is follow-up
time, not survival time.** Every function here takes (time, event) and treats event=0 as right-censored.
Regressing on time directly, or dropping the censored patients, would both be wrong — the first pretends
censored patients died at their last contact, the second throws away the good outcomes.
"""
from __future__ import annotations
import numpy as np


# ------------------------------------------------------------------ Cox ----
def _order(time, event):
    """Sort by time ascending; deaths before censorings at a tie (Breslow convention)."""
    idx = np.lexsort((event == 0, time))
    return idx


def cox_neg_loglik(beta, X, time, event, alpha):
    """Ridge-penalised negative Breslow partial log-likelihood and its gradient."""
    eta = X @ beta
    eta = eta - eta.max()                      # shift for numerical stability; cancels in the ratio
    exp_eta = np.exp(eta)
    # risk set sums computed by reverse cumulative sum over time-sorted data
    csum = np.cumsum(exp_eta[::-1])[::-1]
    ev = event.astype(bool)
    ll = np.sum(eta[ev] - np.log(csum[ev]))
    # gradient: sum_i d_i * (x_i - weighted mean of risk set)
    wx = np.cumsum((exp_eta[:, None] * X)[::-1], axis=0)[::-1]
    mean_x = wx / csum[:, None]
    grad = np.sum(X[ev] - mean_x[ev], axis=0)
    n_ev = max(1, int(ev.sum()))
    return (-ll / n_ev + alpha * beta @ beta,
            -grad / n_ev + 2.0 * alpha * beta)


class CoxPH:
    """Ridge-penalised Cox proportional hazards with a Breslow baseline hazard.

    `alpha` is per-event-normalised, so the same value behaves comparably across blocks of very
    different width — which matters here because the RNA block is 100 PCs and the clinical block is 6.
    """

    def __init__(self, alpha=0.1, max_iter=400):
        self.alpha = float(alpha)
        self.max_iter = int(max_iter)
        self.beta = None
        self.mu = None
        self.sd = None
        self.baseline_t = None
        self.baseline_H = None

    def fit(self, X, time, event):
        from scipy.optimize import minimize
        X = np.asarray(X, float)
        self.mu, self.sd = X.mean(0), X.std(0)
        self.sd[self.sd == 0] = 1.0
        Z = (X - self.mu) / self.sd
        t, e = np.asarray(time, float), np.asarray(event, int)
        o = _order(t, e)
        Z, t, e = Z[o], t[o], e[o]
        b0 = np.zeros(Z.shape[1])
        r = minimize(cox_neg_loglik, b0, args=(Z, t, e, self.alpha), jac=True,
                     method="L-BFGS-B", options={"maxiter": self.max_iter})
        self.beta = r.x
        self._breslow(Z, t, e)
        return self

    def _breslow(self, Z, t, e, ):
        """Cumulative baseline hazard H0(t) = sum over event times of d_j / sum_{risk set} exp(eta)."""
        eta = Z @ self.beta
        exp_eta = np.exp(eta - eta.max())
        csum = np.cumsum(exp_eta[::-1])[::-1]
        ts, Hs, H = [], [], 0.0
        for i in range(len(t)):
            if e[i]:
                H += 1.0 / csum[i]
                ts.append(t[i]); Hs.append(H)
        # exp(-eta.max()) was folded into exp_eta, so fold it back out of H0
        scale = np.exp(-eta.max())
        self.baseline_t = np.asarray(ts, float)
        self.baseline_H = np.asarray(Hs, float) * scale

    def risk(self, X):
        """Linear predictor (log hazard ratio). Higher = worse."""
        Z = (np.asarray(X, float) - self.mu) / self.sd
        return Z @ self.beta

    def survival(self, X, times):
        """S(t | x) = exp(-H0(t) * exp(eta)) for each requested t -> (n_samples, n_times)."""
        eta = self.risk(X)
        H0 = np.interp(np.asarray(times, float), self.baseline_t, self.baseline_H,
                       left=0.0, right=self.baseline_H[-1] if len(self.baseline_H) else 0.0)
        return np.exp(-np.outer(np.exp(eta), H0))

    def median_survival(self, X, cap=None):
        """First time S(t|x) drops to 0.5. NaN when the curve never reaches it inside follow-up —
        which is the honest answer, not a number extrapolated past the data."""
        if self.baseline_t is None or not len(self.baseline_t):
            return np.full(len(X), np.nan)
        grid = self.baseline_t
        S = self.survival(X, grid)
        out = np.full(S.shape[0], np.nan)
        for i in range(S.shape[0]):
            below = np.where(S[i] <= 0.5)[0]
            if len(below):
                out[i] = grid[below[0]]
        if cap is not None:
            out[out > cap] = np.nan
        return out


# -------------------------------------------------------------- metrics ----
def c_index(time, event, risk):
    """Harrell's concordance. Comparable pairs are those where the earlier time is an event."""
    t, e, r = (np.asarray(x, float) for x in (time, event, risk))
    n = len(t)
    conc = disc = tied = 0.0
    for i in range(n):
        if not e[i]:
            continue
        m = t > t[i]                                    # i died first -> i should have the higher risk
        if not m.any():
            continue
        conc += np.sum(r[i] > r[m])
        disc += np.sum(r[i] < r[m])
        tied += np.sum(r[i] == r[m])
    denom = conc + disc + tied
    return float((conc + 0.5 * tied) / denom) if denom else float("nan")


def km(time, event, grid=None):
    """Kaplan-Meier estimate. Returns (times, survival)."""
    t, e = np.asarray(time, float), np.asarray(event, int)
    o = np.argsort(t)
    t, e = t[o], e[o]
    uniq = np.unique(t[e == 1])
    S, out = 1.0, []
    for u in uniq:
        at_risk = np.sum(t >= u)
        d = np.sum((t == u) & (e == 1))
        if at_risk > 0:
            S *= (1.0 - d / at_risk)
        out.append(S)
    ts, Ss = np.asarray(uniq, float), np.asarray(out, float)
    if grid is None:
        return ts, Ss
    return np.asarray(grid, float), np.interp(grid, ts, Ss, left=1.0,
                                              right=Ss[-1] if len(Ss) else 1.0)


def logrank(time, event, group):
    """Two-group log-rank chi-square (1 df) and its p-value."""
    from scipy.stats import chi2
    t, e, g = np.asarray(time, float), np.asarray(event, int), np.asarray(group, int)
    O_E, V = 0.0, 0.0
    for u in np.unique(t[e == 1]):
        n1 = np.sum((t >= u) & (g == 1)); n0 = np.sum((t >= u) & (g == 0))
        n = n1 + n0
        d1 = np.sum((t == u) & (e == 1) & (g == 1)); d = np.sum((t == u) & (e == 1))
        if n <= 1 or d == 0:
            continue
        O_E += d1 - d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    stat = (O_E ** 2 / V) if V > 0 else 0.0
    return float(stat), float(chi2.sf(stat, 1))


def _censoring_km(time, event):
    """KM of the CENSORING distribution G(t) — the inverse-probability weights for Brier/AUC."""
    return km(time, 1 - np.asarray(event, int))


def ipcw_brier(time, event, surv_prob, horizon):
    """IPCW Brier score at one horizon. `surv_prob` = predicted P(alive at horizon)."""
    t, e = np.asarray(time, float), np.asarray(event, int)
    gt, gs = _censoring_km(t, e)
    def G(x):
        return np.interp(x, gt, gs, left=1.0, right=gs[-1] if len(gs) else 1.0)
    Gh = max(G(horizon), 1e-8)
    num = 0.0
    for i in range(len(t)):
        if t[i] <= horizon and e[i] == 1:                       # died before the horizon
            gi = max(G(t[i]), 1e-8)
            num += (surv_prob[i] - 0.0) ** 2 / gi
        elif t[i] > horizon:                                    # known alive at the horizon
            num += (surv_prob[i] - 1.0) ** 2 / Gh
        # censored before the horizon contributes nothing (its weight is carried by the others)
    return float(num / len(t))


def td_auc(time, event, risk, horizon):
    """Time-dependent AUC at a horizon: cases = died by then, controls = known alive past it.

    Patients censored before the horizon are excluded rather than guessed at, so the effective n is
    reported alongside — a time-dependent AUC on 40 evaluable patients means much less than on 300.
    """
    from sklearn.metrics import roc_auc_score
    t, e, r = np.asarray(time, float), np.asarray(event, int), np.asarray(risk, float)
    case = (t <= horizon) & (e == 1)
    ctrl = t > horizon
    m = case | ctrl
    y = case[m].astype(int)
    if y.sum() < 5 or (1 - y).sum() < 5:
        return None, int(m.sum())
    return float(roc_auc_score(y, r[m])), int(m.sum())


# --------------------------------------------------------- self-checks ----
def _selftest():
    rng = np.random.RandomState(0)
    n, true_beta = 4000, 1.0
    x = rng.normal(size=n)
    # exponential survival with hazard exp(beta*x); independent uniform censoring
    T = rng.exponential(1.0 / np.exp(true_beta * x))
    C = rng.exponential(2.0, size=n)
    t = np.minimum(T, C); e = (T <= C).astype(int)
    m = CoxPH(alpha=1e-6).fit(x.reshape(-1, 1), t, e)
    print("1. Cox recovers a known coefficient: beta=%.3f (true %.1f), %d events"
          % (m.beta[0], true_beta, e.sum()))
    assert abs(m.beta[0] - true_beta) < 0.12, m.beta

    # C-index: risk perfectly ordered with (negative) survival time -> 1.0; reversed -> 0.0
    t2 = np.arange(1, 51, dtype=float); e2 = np.ones(50, int)
    print("2. C-index perfect / reversed / random: %.2f / %.2f / %.2f"
          % (c_index(t2, e2, -t2), c_index(t2, e2, t2), c_index(t2, e2, np.zeros(50))))
    assert c_index(t2, e2, -t2) == 1.0 and c_index(t2, e2, t2) == 0.0
    assert c_index(t2, e2, np.zeros(50)) == 0.5

    # Kaplan-Meier, hand-worked: times 1,2,3,4,5 with a censoring at 3
    t3 = np.array([1., 2., 3., 4., 5.]); e3 = np.array([1, 1, 0, 1, 1])
    ts, Ss = km(t3, e3)
    expect = [4 / 5, 4 / 5 * 3 / 4, 4 / 5 * 3 / 4 * 1 / 2, 4 / 5 * 3 / 4 * 1 / 2 * 0.0]
    print("3. KM at event times %s -> %s (expected %s)"
          % (ts.tolist(), np.round(Ss, 4).tolist(), np.round(expect, 4).tolist()))
    assert np.allclose(Ss, expect), (Ss, expect)

    # log-rank: identical groups -> large p; strongly separated -> tiny p
    g = np.repeat([0, 1], 60)
    same = np.concatenate([rng.exponential(1, 60), rng.exponential(1, 60)])
    diff = np.concatenate([rng.exponential(1, 60), rng.exponential(6, 60)])
    ev = np.ones(120, int)
    print("4. log-rank p, same groups %.2f ; separated %.2e"
          % (logrank(same, ev, g)[1], logrank(diff, ev, g)[1]))
    assert logrank(same, ev, g)[1] > 0.05 and logrank(diff, ev, g)[1] < 1e-4

    # survival curves: a higher-risk patient must sit below a lower-risk one everywhere
    S = m.survival(np.array([[-1.0], [1.0]]), np.linspace(0.05, 1.5, 20))
    print("5. S(t) low-risk vs high-risk at t=0.5: %.3f vs %.3f" % (S[0][6], S[1][6]))
    assert np.all(S[0] >= S[1])

    # IPCW Brier: a perfect predictor should beat a constant-0.5 predictor
    sp_perfect = (t > 0.5).astype(float)
    b_perfect = ipcw_brier(t, e, sp_perfect, 0.5)
    b_flat = ipcw_brier(t, e, np.full(n, 0.5), 0.5)
    print("6. IPCW Brier at t=0.5: perfect %.3f < uninformative %.3f" % (b_perfect, b_flat))
    assert b_perfect < b_flat
    print("\nALL SURVIVAL SELF-CHECKS PASSED")


if __name__ == "__main__":
    _selftest()
