"""Per-chain observation analysis + copy-number sanity check.

Two fixes, both prompted by 30ZS:

1. UNOBSERVED STATS PER CHAIN, NOT UNION.
   Taking the union of unobserved ranges across chain copies is the wrong metric.
   30ZS chain A is 223/246 observed, while copies C and D have 3 residues each;
   the union of their gaps covers all 246 positions, so the union reported "100%
   unobserved" for a structure that is in fact nearly complete. Reported now:
   the BEST-resolved instance of each entity, which is what a student compares to.

2. GHOST COPIES IN THE BIOLOGICAL ASSEMBLY.
   30ZS declares 3 copies in assembly 1, but two are near-empty. Trypsin is a
   monomer; those copies are lattice artifacts that got swept into the assembly
   definition. Modelling 3 full-length chains against that is meaningless, so a
   copy whose best instance is mostly unobserved is not counted as a real copy.

Threshold: an instance counts as a real copy if >=50% of its construct residues
have coordinates. Chains failing this are reported explicitly, not silently
dropped.
"""
import json
import os

OUT = "/Users/lilindu/pdb-structures-for-homework/alphafold_homework_41"
raw = json.load(open("seqgap_raw.json"))
recs = json.load(open("targets41_full.json"))

MIN_OBSERVED_FRACTION = 0.50


def instance_stats(eid):
    """entity_id -> list of {chain, n_unobs, construct_len, frac_observed}"""
    out = {}
    e = raw.get(eid)
    if not e:
        return out
    for pe in e.get("polymer_entities") or []:
        ent = pe["rcsb_id"].split("_")[-1]
        ep = pe.get("entity_poly") or {}
        seq = (ep.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").strip()
        clen = len(seq)
        insts = []
        for pi in pe.get("polymer_entity_instances") or []:
            ch = ((pi.get("rcsb_polymer_entity_instance_container_identifiers") or {})
                  .get("auth_asym_id"))
            miss = 0
            ranges = []
            for f in pi.get("rcsb_polymer_instance_feature") or []:
                if f.get("type") == "UNOBSERVED_RESIDUE_XYZ":
                    for fp in f.get("feature_positions") or []:
                        b = fp.get("beg_seq_id")
                        en = fp.get("end_seq_id") or b
                        if b:
                            miss += en - b + 1
                            ranges.append((b, en))
            obs = clen - miss
            insts.append({"chain": ch, "n_unobs": miss, "construct_len": clen,
                          "observed": obs,
                          "frac": (obs / clen) if clen else 0.0,
                          "ranges": ranges})
        insts.sort(key=lambda x: -x["frac"])
        out[ent] = {"construct_len": clen, "instances": insts,
                    "polymer_type": ep.get("rcsb_entity_polymer_type")}
    return out


print("=== chains whose copies are not all real (ghost copies in assembly) ===")
fixes = []
report = {}
for r in recs:
    st = instance_stats(r["id"])
    report[r["id"]] = st
    for c in r["chains"]:
        d = st.get(c["entity"])
        if not d or not d["instances"]:
            continue
        real = [i for i in d["instances"] if i["frac"] >= MIN_OBSERVED_FRACTION]
        ghosts = [i for i in d["instances"] if i["frac"] < MIN_OBSERVED_FRACTION]
        if ghosts and c["copies"] > len(real) and len(real) >= 1:
            fixes.append((r["no"], r["id"], c["entity"], c["copies"], len(real),
                          [(g["chain"], g["observed"], g["construct_len"]) for g in ghosts]))
            print(f"#{r['no']:2d} {r['id']} e{c['entity']}  copies {c['copies']} -> {len(real)}"
                  f"   {c['desc'][:40]}")
            for g in ghosts:
                print(f"      幽灵拷贝 链{g['chain']}: 仅 {g['observed']}/{g['construct_len']} 个残基有坐标")

json.dump(report, open("perchain.json", "w"))

# ---- apply fixes and recompute ----
fixmap = {(f[1], f[2]): f[4] for f in fixes}
new = []
for r in recs:
    q = dict(r)
    chains = []
    for c in r["chains"]:
        c = dict(c)
        k = (r["id"], c["entity"])
        if k in fixmap:
            c["copies_asu_declared"] = c["copies"]
            c["copies"] = fixmap[k]
            c["copies_note"] = "装配中的空拷贝已剔除(仅少数残基有坐标)"
        d = (report.get(r["id"]) or {}).get(c["entity"])
        if d and d["instances"]:
            best = d["instances"][0]
            c["observed_best"] = best["observed"]
            c["construct_len_ref"] = best["construct_len"]
            c["frac_observed"] = round(best["frac"], 3)
            c["unobserved_ranges"] = best["ranges"]
        chains.append(c)
    q["chains"] = chains
    tok = sum(c["len"] * c["copies"] for c in chains)
    tok += sum((l["atoms"] or 0) * l["count"] for l in r["ligands"])
    tok += sum(i["count"] for i in r["ions"])
    q["tokens_before_ghostfix"] = r["tokens"]
    q["tokens"] = tok
    new.append(q)

print(f"\ncopy-number fixes applied: {len(fixes)}")
for q in new:
    if q["tokens"] != q["tokens_before_ghostfix"]:
        print(f"   #{q['no']} {q['id']}: token {q['tokens_before_ghostfix']} -> {q['tokens']}")

print("\n=== 观测比例最低的链(按最完整的那份拷贝算) ===")
rows = []
for q in new:
    for c in q["chains"]:
        if c["type"] == "Protein" and c.get("frac_observed") is not None:
            rows.append((c["frac_observed"], q["no"], q["id"], c["entity"], c))
rows.sort()
for frac, no, pid, _ent, c in rows[:14]:
    print(f"#{no:2d} {pid} e{c['entity']}  有坐标 {c['observed_best']}/{c['construct_len_ref']}"
          f" = {frac*100:5.1f}%   {c['desc'][:40]}")

json.dump(new, open("targets41_final.json", "w"), ensure_ascii=False, indent=1)
print(f"\ntoken range: {min(q['tokens'] for q in new)}-{max(q['tokens'] for q in new)}")
print("over 5000:", sum(1 for q in new if q["tokens"] > 5000))
print("wrote targets41_final.json")
