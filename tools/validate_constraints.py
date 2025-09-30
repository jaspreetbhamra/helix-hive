import json
import sys

exp = json.load(open(sys.argv[1]))
act = json.load(open(sys.argv[2]))

ok = True

# cysteines: must match exactly
if act.get("cysteines", []) != exp.get("cysteines", []):
    print("❌ cysteines mismatch")
    ok = False

# minimal conserved: expected ⊆ actual
exp_cons = set(exp.get("conserved_sites", []))
act_cons = set(act.get("conserved_sites", []))
missing = exp_cons - act_cons
if missing:
    print(f"❌ conserved_sites missing required: {sorted(missing)}")
    ok = False

# signal peptide & TM spans exact
for k in ("signal_peptide",):
    if act.get(k) != exp.get(k):
        print(f"❌ {k} mismatch")
        ok = False
if act.get("transmembrane_spans", []) != exp.get("transmembrane_spans", []):
    print("❌ transmembrane_spans mismatch")
    ok = False

print("✅ PASS" if ok else "❌ FAIL")
sys.exit(0 if ok else 1)
