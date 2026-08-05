#!/usr/bin/env python3
"""Wait for one ingest run to finish, then add the AGENT's therapy explanations to its report.

Why this exists: the pipeline is dispatched via bsub onto a COMPUTE node, and compute nodes are
firewalled from the LLM gateway (port 4000 answers only from the login node). So a freshly ingested
patient would land with deterministic [RULES] text instead of the [AGENT] narration. gui_server runs ON
the login node, so it spawns this watcher; the watcher polls the run's status.json and, once the report
exists, runs the explanation pass for just that run — where the gateway is reachable.

  python await_and_explain.py <run_dir> <run_id> [timeout_s]
"""
import os, sys, json, time, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
run_dir = sys.argv[1]
run_id = sys.argv[2]
limit = int(sys.argv[3]) if len(sys.argv) > 3 else 5400          # give LSF time to queue + run

report = os.path.join(run_dir, "patient_report.json")
status = os.path.join(run_dir, "status.json")
t0 = time.time()
while time.time() - t0 < limit:
    st = None
    try:
        with open(status, encoding="utf-8") as fh:
            st = (json.load(fh) or {}).get("state")
    except Exception:
        pass
    if st == "error":
        sys.exit(0)                                              # nothing to explain
    if os.path.exists(report) and st == "done":
        break
    time.sleep(10)
else:
    sys.exit(0)                                                  # timed out; report keeps its rules text

runs_root = os.path.dirname(os.path.abspath(run_dir))
# no FORCE: explain_reports skips rows already agent-explained, so this stays idempotent
env = dict(os.environ, ONLY=run_id)
try:
    subprocess.run([sys.executable, os.path.join(HERE, "explain_reports.py"), runs_root],
                   env=env, timeout=1800,
                   stdout=open(os.path.join(run_dir, "_explain.log"), "w"), stderr=subprocess.STDOUT)
except Exception:
    pass                                                         # report already has deterministic text
