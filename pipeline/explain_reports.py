#!/usr/bin/env python3
"""Add the tumor-board AGENT's explanation to every report's therapy panel.

WHY THIS IS A SEPARATE PASS (and not just part of ingest):
the LiteLLM gateway listens on the LOGIN node (0.0.0.0:4000, answers 200 there), but COMPUTE nodes are
firewalled off port 4000. The pipeline runs via bsub on compute nodes, so its LLM calls time out and it
silently degrades to the deterministic explanation. (The arbiter's narration has the same blind spot.)
So the narration runs HERE, on the login node, where the gateway is actually reachable.

It needs no numpy: therapy.py and llm.py are pure json+urllib, and they're loaded by FILE PATH to skip
`amlmm/__init__` (which imports pandas/numpy and would die on the login node's old glibc). That means it
runs with the login node's system python3:

    /usr/local/anaconda3-2020/bin/python explain_reports.py [runs_dir]

(the login node system python3 is 3.6 — too old for `from __future__ import annotations`; the conda
python is 3.12 and its stdlib works fine on the login node, only numpy breaks there.)

Idempotent: rows already explained by the agent are left alone unless FORCE=1.
"""
import os, sys, json, glob, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(HERE), "runs")
FORCE = os.environ.get("FORCE") == "1"
ONLY = os.environ.get("ONLY")


def _load(name, path):
    """Import a module straight from its file, bypassing the package __init__ (which needs numpy)."""
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


therapy = _load("_therapy", os.path.join(HERE, "amlmm", "therapy.py"))
llmmod = _load("_llm", os.path.join(HERE, "amlmm", "llm.py"))

client = None
try:
    client = llmmod.LLMClient()
    client.chat("Reply with the single word: ok", max_tokens=64)
    print("LLM gateway reachable: %s (%s)" % (client.base_url, client.model))
except Exception as e:
    print("LLM gateway NOT reachable (%s) — nothing to upgrade; reports keep their deterministic text."
          % str(e)[:90])
    sys.exit(1)

done, skipped, failed = [], 0, []
for f in sorted(glob.glob(os.path.join(RUNS, "*", "patient_report.json"))):
    run = os.path.basename(os.path.dirname(f))
    if ONLY and not run.startswith(ONLY):
        continue
    if run.startswith("leucegene_") and not (ONLY and ONLY.startswith("leucegene")):
        continue                                        # 367 validation samples: keep deterministic text
    try:
        with open(f, encoding="utf-8") as fh:
            rep = json.load(fh)
    except Exception as e:
        failed.append((run, "unreadable: %s" % e)); continue
    tp = rep.get("treatment_panel")
    if not isinstance(tp, list) or not tp:
        continue
    if not FORCE and all(t.get("explained_by") == "agent" for t in tp):
        skipped += 1; continue
    ctx = {"sample": rep.get("sample_key"), "dataset": rep.get("dataset"),
           "specimen_class": rep.get("specimen_class"),
           "input": (rep.get("ingest") or {}).get("input_kind") or "single-cell"}
    try:
        therapy.explain_treatments(client, tp, ctx)      # rules already decided; agent only explains
    except Exception as e:
        failed.append((run, "explain: %s" % e)); continue
    rep["treatment_panel"] = tp
    with open(f, "w", encoding="utf-8") as fh:
        json.dump(rep, fh, default=str, indent=1)
    n_agent = sum(1 for t in tp if t.get("explained_by") == "agent")
    done.append((run, n_agent, len(tp)))

print("\nexplained %d report(s) (%d already done, skipped):" % (len(done), skipped))
for run, na, n in done:
    print("   %-44s agent=%d/%d" % (run, na, n))
if failed:
    print("\n%d FAILED:" % len(failed))
    for run, why in failed:
        print("   %-44s %s" % (run, why))
print("\nEXPLAIN OK")
