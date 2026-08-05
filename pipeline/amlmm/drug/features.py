"""Patient feature blocks for the drug-response model, in a space computable on BOTH assay types.

Everything here is a function of a single expression vector in the shared 14,237-ENSG space plus
optional metadata, which is what lets one trained model score a BeatAML bulk specimen, an uploaded
patient's single-cell bulk-equivalent, and one cell state's pseudobulk without retraining.

Blocks (each independently ablatable, because "which modality is carrying this?" is a question the
validation has to answer, not assert):

  rna     PCA of the cohort-z-scored transcriptome -- the bulk view
  state   10 atlas-derived lineage-signature scores + a primitive->mature axis.
          BeatAML2's own analysis found differentiation state to be a broad determinant of ex-vivo
          response, so this block exists to give the model that axis explicitly rather than hoping the
          RNA PCs happen to encode it.
  mut     variant-level driver indicators. OBSERVED for BeatAML; at inference for a new single-cell
          patient these are the mutation layer's calibrated probabilities, which is why the block is
          kept separate and why every headline number is also reported without it.
  clin    age, sex, blast counts, WBC, specimen type, disease stage, prior-disease flags, ELN2017.

The PCA and all imputation constants are FIT ON TRAINING SPECIMENS ONLY and carried in the object; a
new sample is transformed, never re-fit.
"""
from __future__ import annotations
import os, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SIG_PATH = os.path.join(ROOT, "labels", "cellstate_signatures.json")

BLOCKS = ["rna", "state", "mut", "clin"]

_CLIN_NUM = ["ageAtSpecimenAcquisition", "%.Blasts.in.BM", "%.Blasts.in.PB", "wbcCount"]
_CLIN_CAT = {
    "consensus_sex": ["Male"],
    "specimenType": ["Bone Marrow Aspirate", "Peripheral Blood", "Leukapheresis"],
    "diseaseStageAtSpecimenCollection": ["Initial Diagnosis", "Relapse", "Residual"],
    "ELN2017": ["Favorable", "Intermediate", "Adverse"],
    "isRelapse": ["TRUE"], "isDenovo": ["TRUE"], "isTransformed": ["TRUE"], "priorMDS": ["y"],
}


def load_signatures(path=SIG_PATH):
    with open(path) as f:
        return json.load(f)


class FeatureSpace:
    """Fit on training specimens; transform anything in the same gene space."""

    def __init__(self, genes, n_pc=100, var_genes=4000, sig_path=SIG_PATH):
        self.genes = list(genes)
        self.n_pc = n_pc
        self.var_genes = var_genes
        self.sig = load_signatures(sig_path)
        self.sig_idx = {}
        gidx = {g: i for i, g in enumerate(self.genes)}
        for name, ens in self.sig["signatures"].items():
            ii = [gidx[e] for e in ens if e in gidx]
            if ii:
                self.sig_idx[name] = np.array(ii)
        self.state_names = sorted(self.sig_idx)
        self.ref = {}            # cohort -> (mu, sd) on the common log2 scale
        self.pca = None
        self.sel = None
        self.clin_med = {}
        self.mut_cols = []

    # ---------------------------------------------------------- scaling ----
    @staticmethod
    def clog(x_lin):
        return np.log2(np.clip(np.asarray(x_lin, float), 0, None) + 1.0)

    def add_reference(self, name, X_lin):
        """Store a cohort's own mean/sd. Cross-assay batch is handled the same way the mutation caller
        handles it: put every cohort on a common log2 scale, then z-score WITHIN cohort."""
        cl = self.clog(X_lin)
        mu, sd = cl.mean(0), cl.std(0)
        sd[sd == 0] = 1.0
        self.ref[name] = (mu, sd)
        return (cl - mu) / sd

    def z(self, X_lin, ref):
        mu, sd = self.ref[ref]
        return (self.clog(X_lin) - mu) / sd

    # ------------------------------------------------------------ blocks ----
    def _state_block(self, Z):
        cols, names = [], []
        for n in self.state_names:
            cols.append(Z[:, self.sig_idx[n]].mean(1)); names.append("state_" + n)
        S = np.vstack(cols).T
        by = {n: S[:, i] for i, n in enumerate(names)}
        prim = np.mean([by["state_" + k] for k in self.sig["axes"]["primitive"] if "state_" + k in by], 0)
        mat = np.mean([by["state_" + k] for k in self.sig["axes"]["mature"] if "state_" + k in by], 0)
        extra = np.vstack([prim, mat, prim - mat]).T          # the differentiation axis, made explicit
        return np.hstack([S, extra]), names + ["axis_primitive", "axis_mature", "axis_prim_minus_mature"]

    def fit(self, Z_train, meta_train=None, mut_train=None):
        from sklearn.decomposition import PCA
        v = Z_train.var(0)
        self.sel = np.argsort(v)[::-1][:min(self.var_genes, Z_train.shape[1])]
        npc = int(min(self.n_pc, Z_train.shape[0] - 1, len(self.sel)))
        self.pca = PCA(n_components=npc, random_state=0).fit(Z_train[:, self.sel])
        if meta_train is not None:
            for c in _CLIN_NUM:
                v = pd.to_numeric(meta_train.get(c), errors="coerce") if c in meta_train else None
                self.clin_med[c] = float(np.nanmedian(v)) if v is not None and v.notna().any() else 0.0
        if mut_train is not None:
            self.mut_cols = list(mut_train.columns)
        return self

    def _clin_block(self, meta, n):
        cols, names = [], []
        for c in _CLIN_NUM:
            v = pd.to_numeric(meta[c], errors="coerce") if (meta is not None and c in meta) else pd.Series([np.nan] * n)
            med = self.clin_med.get(c, 0.0)
            miss = v.isna().astype(float).values
            cols += [v.fillna(med).values.astype(float), miss]
            names += [c, c + "__missing"]
        for c, levels in _CLIN_CAT.items():
            s = meta[c].astype(str) if (meta is not None and c in meta) else pd.Series(["nan"] * n)
            for lv in levels:
                cols.append(s.str.lower().eq(str(lv).lower()).astype(float).values)
                names.append("%s=%s" % (c, lv))
        C = np.vstack(cols).T
        # standardise the numeric columns so ridge treats age and %blasts comparably
        for i, nm in enumerate(names):
            if nm in _CLIN_NUM:
                sd = C[:, i].std() or 1.0
                C[:, i] = (C[:, i] - C[:, i].mean()) / sd
        return C, names

    def transform(self, Z, meta=None, mut=None, blocks=BLOCKS):
        """Z: already cohort-z-scored expression (n x genes). Returns (X, names, block_slices)."""
        parts, names = [], []
        n = Z.shape[0]
        if "rna" in blocks:
            P = self.pca.transform(Z[:, self.sel])
            parts.append(P); names += ["rna_pc%d" % (i + 1) for i in range(P.shape[1])]
        if "state" in blocks:
            S, sn = self._state_block(Z)
            parts.append(S); names += sn
        if "mut" in blocks:
            if mut is None:
                M = np.zeros((n, len(self.mut_cols)))
            else:
                M = mut.reindex(columns=self.mut_cols).astype(float).fillna(0.0).values
            parts.append(M); names += ["mut_" + c for c in self.mut_cols]
        if "clin" in blocks:
            C, cn = self._clin_block(meta, n)
            parts.append(C); names += ["clin_" + c for c in cn]
        X = np.hstack(parts) if parts else np.zeros((n, 0))
        sl, i = {}, 0
        for b, p in zip([b for b in BLOCKS if b in blocks], parts):
            sl[b] = slice(i, i + p.shape[1]); i += p.shape[1]
        return X, names, sl

    # -------------------------------------------------- drug interaction ----
    def target_features(self, Z, targets, sym2ens=None):
        """z-expression of one drug's annotated target genes: the (patient x drug) interaction term
        that lets a single shared model generalise across inhibitors instead of memorising each one."""
        gidx = getattr(self, "_gidx", None)
        if gidx is None:
            gidx = self._gidx = {g: i for i, g in enumerate(self.genes)}
        ii = []
        for t in targets:
            e = t if t in gidx else (sym2ens or {}).get(t)
            if e in gidx:
                ii.append(gidx[e])
        if not ii:
            return np.zeros((Z.shape[0], 3)), ["tgt_mean", "tgt_max", "tgt_n"]
        sub = Z[:, np.array(ii)]
        return (np.vstack([sub.mean(1), sub.max(1), np.full(Z.shape[0], len(ii) / 10.0)]).T,
                ["tgt_mean", "tgt_max", "tgt_n"])
