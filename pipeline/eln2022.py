#!/usr/bin/env python3
"""Assign ELN 2022 risk to BeatAML2 specimens — ONLY where it is determinable with high confidence.

Scope decision (deliberate): a cohort-wide ELN 2022 re-derivation is NOT possible here. ELN risk is a
diagnosis-time label, and the adverse tier turns on ISCN karyotype parsing (complex / monosomal /
-5,del(5q) / -7 / -17,abn(17p)) that BeatAML only provides as free-text strings — 500/842 are non-trivial,
multi-clone, sometimes sex-discordant. BeatAML's own expert curators left 27 eligible cases ambiguous.
So this ASSIGNS where the answer is forced by data we hold cleanly, and ABSTAINS everywhere else.

The key insight that makes a useful subset possible: adverse MOLECULAR markers are decisive regardless of
karyotype. If TP53 (VAF>=10%) or an MDS-related gene is mutated and no favorable-defining lesion and no
t(9;11) is present, unparsed cytogenetics can only reinforce Adverse — never rescue it. Likewise the
curated `consensusAMLFusions` column gives the class-defining rearrangements without ISCN parsing.

Rules + footnotes taken verbatim from Dohner/Lowenberg, Blood 2022 (in the supplied zip):
  * Favorable   : t(8;21)/RUNX1::RUNX1T1; inv(16)|t(16;16)/CBFB::MYH11; NPM1mut w/o FLT3-ITD;
                  bZIP in-frame CEBPA
  * Intermediate: NPM1mut with FLT3-ITD; NPM1wt with FLT3-ITD (no adverse lesions);
                  t(9;11)/MLLT3::KMT2A; anything not favorable or adverse
  * Adverse     : t(6;9)/DEK::NUP214; t(v;11q23.3)/KMT2A-r (excl. PTD); t(9;22)/BCR::ABL1;
                  t(8;16)/KAT6A::CREBBP; inv(3)|t(3;3)/GATA2,MECOM; t(3q26.2;v)/MECOM-r;
                  -5|del(5q); -7; -17|abn(17p); complex; monosomal;
                  mutated ASXL1/BCOR/EZH2/RUNX1/SF3B1/SRSF2/STAG2/U2AF1/ZRSR2; TP53 VAF>=10%
  footnotes honoured:
    (S) NPM1 + adverse-risk cytogenetics -> Adverse  => a favourable NPM1 call needs a NORMAL karyotype
    (||) only IN-FRAME bZIP CEBPA counts, mono- or biallelic
    (P) t(9;11) takes precedence over rare concurrent adverse-risk GENE mutations
    (++) MDS-related gene mutations are NOT adverse when they co-occur with a favorable-risk subtype
    (‡) concurrent KIT/FLT3 does not alter CBF risk
    (a) TP53 at VAF >= 10%, irrespective of allelic status

  python eln2022.py            # writes labels/eln2022_beataml.tsv + prints the audit
"""
import os, re, sys, warnings
warnings.filterwarnings("ignore")
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BA = os.path.join(ROOT, "data", "external", "beataml")
OUT = os.path.join(ROOT, "labels", "eln2022_beataml.tsv")

MDS_GENES = ["ASXL1", "BCOR", "EZH2", "RUNX1", "SF3B1", "SRSF2", "STAG2", "U2AF1", "ZRSR2"]
FAV_FUS = {"CBFB-MYH11", "RUNX1-RUNX1T1"}
ADV_FUS = {"GATA2-MECOM", "DEK-NUP214", "BCR-ABL1", "KMT2A_re"}      # KMT2A_re: t(v;11q23.3); PTD is not a fusion
T911 = {"MLLT3-KMT2A"}
APL = {"PML-RARA"}
NORMAL_KARYO = re.compile(r"^46,X[XY]\[\d+\]$")                       # a single normal clone, unambiguous
CEBPA_BZIP = (278, 358)                                              # basic leucine zipper region


def _pos(h):
    m = re.match(r"^p\.[A-Za-z*](\d+)", str(h))
    return int(m.group(1)) if m else None


def load():
    cl = pd.read_excel(os.path.join(BA, "clinical.xlsx"), "summary")
    mut = pd.read_csv(os.path.join(BA, "mutations.txt"), sep="\t",
                      usecols=["dbgap_sample_id", "symbol", "hgvsp_short", "variant_classification", "t_vaf"])
    mut["dbgap_sample_id"] = mut["dbgap_sample_id"].astype(str)
    return cl, mut


def build(cl, mut):
    wes = set(mut["dbgap_sample_id"].unique())
    by_s = {s: g for s, g in mut.groupby("dbgap_sample_id")}

    rows = []
    for _, r in cl.iterrows():
        dna = str(r.get("dbgap_dnaseq_sample"))
        sid = str(r.get("dbgap_subject_id"))
        rna = str(r.get("dbgap_rnaseq_sample"))
        stage = str(r.get("diseaseStageAtSpecimenCollection"))
        fus = str(r.get("consensusAMLFusions"))
        ky = str(r.get("karyotype")).strip()
        rec = {"dbgap_dnaseq_sample": dna, "dbgap_rnaseq_sample": rna, "dbgap_subject_id": sid,
               "ELN2017": str(r.get("ELN2017")), "fusion": fus, "karyotype_normal": bool(NORMAL_KARYO.match(ky)),
               "ELN2022": None, "basis": None, "abstain_reason": None}

        # ---- eligibility ----------------------------------------------------------------
        if str(r.get("ELN2017")) == "NonAML":
            rec["abstain_reason"] = "not AML"; rows.append(rec); continue
        if "Initial Diagnosis" not in stage:
            rec["abstain_reason"] = "not an initial-diagnosis specimen (ELN is a diagnosis-time label)"
            rows.append(rec); continue
        if dna not in wes:
            rec["abstain_reason"] = "no WES"; rows.append(rec); continue
        if fus in APL:
            rec["abstain_reason"] = "APL (PML-RARA) - separate entity, not ELN-risk classified"
            rows.append(rec); continue

        g = by_s.get(dna)
        sym = set(g["symbol"].astype(str)) if g is not None else set()

        # ---- molecular facts ------------------------------------------------------------
        npm1 = str(r.get("NPM1")).lower() == "positive"
        itd = str(r.get("FLT3-ITD")).lower() == "positive"
        tp53 = False
        if g is not None:
            t = g[g["symbol"] == "TP53"]
            tp53 = bool((t["t_vaf"] >= 0.10).any())                   # footnote a: VAF >= 10%
        cebpa_bzip = False
        if g is not None:
            c = g[g["symbol"] == "CEBPA"].copy()
            if len(c):
                c["p"] = c["hgvsp_short"].map(_pos)
                inframe = c["variant_classification"].isin(
                    ["missense_variant", "inframe_insertion", "inframe_deletion"])
                cebpa_bzip = bool((inframe & c["p"].between(*CEBPA_BZIP)).any())   # footnote ||
        mds = sorted(sym & set(MDS_GENES))
        cbf = fus in FAV_FUS
        t911 = fus in T911
        adv_fus = fus in ADV_FUS
        fav_lesion = cbf or (npm1 and not itd) or cebpa_bzip

        def done(risk, basis):
            rec["ELN2022"] = risk; rec["basis"] = basis; rows.append(rec)

        def skip(why):
            rec["abstain_reason"] = why; rows.append(rec)

        # ---- classification (precedence follows the guideline) ---------------------------
        # 1. adverse class-defining fusion: decisive, karyotype cannot rescue it
        if adv_fus:
            if fav_lesion:
                skip("conflict: adverse fusion %s + favorable lesion - guideline does not resolve" % fus)
            else:
                done("Adverse", "adverse class-defining fusion %s" % fus)
            continue
        # 2. t(9;11): Intermediate, takes precedence over adverse GENE mutations (footnote P)
        if t911:
            if tp53:
                skip("t(9;11) + TP53 - footnote P covers gene mutations; TP53 precedence unclear")
            else:
                done("Intermediate", "t(9;11)/MLLT3::KMT2A (footnote P: precedence over adverse gene muts)")
            continue
        # 3. TP53 VAF>=10% -> Adverse (unparsed cytogenetics can only reinforce)
        if tp53:
            if fav_lesion:
                skip("TP53 VAF>=10%% + favorable lesion - guideline does not resolve precedence")
            else:
                done("Adverse", "TP53 mutated at VAF>=10%")
            continue
        # 4. CBF -> Favorable (footnote ‡ KIT/FLT3 irrelevant; ++ MDS genes do not override)
        if cbf:
            done("Favorable", "CBF fusion %s%s" % (fus, " (MDS genes present but footnote ++ applies)" if mds else ""))
            continue
        # 5. NPM1 w/o ITD, or bZIP CEBPA -> Favorable, but footnote S needs NO adverse cytogenetics,
        #    so require a definitively NORMAL karyotype
        if fav_lesion:
            if rec["karyotype_normal"]:
                why = "NPM1mut without FLT3-ITD" if (npm1 and not itd) else "bZIP in-frame CEBPA"
                done("Favorable", "%s + normal karyotype%s" % (why, " (footnote ++ over %s)" % ",".join(mds) if mds else ""))
            else:
                skip("favorable molecular lesion but karyotype not definitively normal - "
                     "footnote S (NPM1 + adverse cytogenetics = Adverse) cannot be excluded")
            continue
        # 6. adverse MDS-related gene mutation -> Adverse (no favorable lesion, no t(9;11) reached here)
        if mds:
            done("Adverse", "MDS-related gene mutation: %s" % ",".join(mds))
            continue
        # 7. nothing favorable or adverse molecularly -> Intermediate ONLY if karyotype is normal
        if rec["karyotype_normal"]:
            b = "normal karyotype, no favorable/adverse marker"
            if npm1 and itd:
                b = "NPM1mut with FLT3-ITD + normal karyotype"
            elif itd:
                b = "FLT3-ITD, NPM1 wild-type + normal karyotype, no adverse lesion"
            done("Intermediate", b)
        else:
            skip("no decisive molecular marker and karyotype not definitively normal - "
                 "complex/monosomal/-5/-7/-17 cannot be excluded without ISCN parsing")
    return pd.DataFrame(rows)


def main():
    cl, mut = load()
    df = build(cl, mut)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)

    n = len(df)
    ok = df[df["ELN2022"].notna()]
    print("BeatAML2 specimens: %d" % n)
    print("ELN 2022 ASSIGNED (high confidence): %d  (%d unique patients)"
          % (len(ok), ok["dbgap_subject_id"].nunique()))
    print("ABSTAINED: %d" % (n - len(ok)))
    print()
    print("=== assigned risk distribution ===")
    print(ok["ELN2022"].value_counts().to_string())
    print()
    print("=== what forced each call ===")
    kind = ok["basis"].str.replace(r"[:(].*", "", regex=True).str.strip()
    print(kind.value_counts().to_string())
    print()
    print("=== why the rest abstained ===")
    ab = df[df["ELN2022"].isna()]["abstain_reason"].str.replace(r" -.*", "", regex=True)
    print(ab.value_counts().to_string())
    print()
    print("=== concordance with BeatAML's own ELN2017 (expect real 2017->2022 shifts) ===")
    m = ok[ok["ELN2017"].isin(["Favorable", "Intermediate", "Adverse"])]
    print(pd.crosstab(m["ELN2017"], m["ELN2022"]).to_string())
    agree = (m["ELN2017"] == m["ELN2022"]).mean() if len(m) else float("nan")
    print("\nagreement with ELN2017: %.1f%% of %d comparable" % (100 * agree, len(m)))
    print("\nwrote %s" % OUT)


if __name__ == "__main__":
    sys.exit(main())
