"""For each of the 41 targets, quantify the three-way sequence discrepancy.

Three distinct sequences get conflated when people say "the PDB sequence":

  A. 全长天然蛋白  -- UniProt canonical sequence (the whole gene product)
  B. 构建体 construct -- entity_poly.pdbx_seq_one_letter_code_can: what the
     experimentalist actually put in the tube. Deliberately truncated to a
     domain, may carry tags, mutations, or be a fusion chimera.
  C. 观测到的残基 -- residues that actually have coordinates. B minus disordered
     loops and flexible termini, which are involuntary losses from resolution.

A -> B is a deliberate choice by the depositor.
B -> C is an involuntary loss.

The job files currently feed B. This script measures how far B sits from A, and
how far C sits from B, so the handout can state both gaps per target.
"""
import json
import os
import re
import time
import urllib.request

OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
G = "https://data.rcsb.org/graphql"
OUT = "/Users/lilindu/alphafold-homework-2026"

Q = """
query($ids:[String!]!){
 entries(entry_ids:$ids){
  rcsb_id
  polymer_entities{
   rcsb_id
   entity_poly{pdbx_seq_one_letter_code_can rcsb_entity_polymer_type}
   rcsb_polymer_entity{pdbx_description}
   rcsb_polymer_entity_align{
    reference_database_name
    reference_database_accession
    aligned_regions{entity_beg_seq_id ref_beg_seq_id length}}
   uniprots{
    rcsb_id
    rcsb_uniprot_protein{name{value} sequence}}
   polymer_entity_instances{
    rcsb_id
    rcsb_polymer_entity_instance_container_identifiers{auth_asym_id}
    rcsb_polymer_instance_feature{
     type
     feature_positions{beg_seq_id end_seq_id}}}}}}
"""


def post(p, t=300):
    r = urllib.request.Request(G, data=json.dumps(p).encode(),
                               headers={"Content-Type": "application/json"})
    return json.load(OP.open(r, timeout=t))


recs = json.load(open(os.path.join(OUT, "targets41.json")))
ids = sorted({r["id"] for r in recs})

raw = {}
if os.path.exists("seqgap_raw.json"):
    raw = json.load(open("seqgap_raw.json"))
todo = [i for i in ids if i not in raw]
print("fetching", len(todo), flush=True)
for i in range(0, len(todo), 12):
    ch = todo[i:i + 12]
    for a in range(4):
        try:
            d = post({"query": Q, "variables": {"ids": ch}})
            if "data" not in d:
                raise RuntimeError(str(d.get("errors"))[:200])
            for e in d["data"]["entries"] or []:
                raw[e["rcsb_id"]] = e
            break
        except Exception as ex:
            print("retry", str(ex)[:100], flush=True)
            time.sleep(4)
json.dump(raw, open("seqgap_raw.json", "w"))
print("fetched", len(raw), flush=True)

TAGS = [
    ("His-tag", re.compile(r"H{6,}")),
    ("Strep-tag II", re.compile(r"WSHPQFEK")),
    ("FLAG", re.compile(r"DYKDDDDK")),
    ("Myc", re.compile(r"EQKLISEEDL")),
    ("TEV site", re.compile(r"ENLYFQ[GS]")),
    ("Thrombin site", re.compile(r"LVPR[GS]S")),
    ("HRV3C site", re.compile(r"LEVLFQGP")),
    ("Avi-tag", re.compile(r"GLNDIFEAQKIEWHE")),
    ("SUMO/Ubl 尾", re.compile(r"GGSGGS{2,}")),
]

report = {}
for r in recs:
    eid = r["id"]
    e = raw.get(eid)
    if not e:
        continue
    ents = {}
    for pe in e.get("polymer_entities") or []:
        ent = pe["rcsb_id"].split("_")[-1]
        ep = pe.get("entity_poly") or {}
        seq = (ep.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").strip().upper()
        ptype = ep.get("rcsb_entity_polymer_type") or "?"
        desc = ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()

        # --- A: full-length reference ---
        ups = []
        for u in pe.get("uniprots") or []:
            s = ((u.get("rcsb_uniprot_protein") or {}).get("sequence") or "")
            nm = (((u.get("rcsb_uniprot_protein") or {}).get("name") or {}).get("value") or "")
            ups.append({"acc": u.get("rcsb_id"), "len": len(s) if s else None, "name": nm})

        # coverage of the reference, from the deposited alignment
        cover, ref_span = None, None
        aligns = pe.get("rcsb_polymer_entity_align") or []
        best = None
        for al in aligns:
            if (al.get("reference_database_name") or "").upper().startswith("UNIPROT"):
                regs = al.get("aligned_regions") or []
                tot = sum(x.get("length") or 0 for x in regs)
                if regs and (best is None or tot > best[0]):
                    lo = min(x["ref_beg_seq_id"] for x in regs)
                    hi = max(x["ref_beg_seq_id"] + (x["length"] or 0) - 1 for x in regs)
                    best = (tot, al.get("reference_database_accession"), lo, hi)
        if best:
            cover = best[0]
            ref_span = (best[1], best[2], best[3])

        # --- C: unobserved residues, per chain instance ---
        unobs_union, inst_detail = set(), []
        for pi in pe.get("polymer_entity_instances") or []:
            ch = ((pi.get("rcsb_polymer_entity_instance_container_identifiers") or {})
                  .get("auth_asym_id"))
            miss = []
            for f in pi.get("rcsb_polymer_instance_feature") or []:
                if f.get("type") == "UNOBSERVED_RESIDUE_XYZ":
                    for fp in f.get("feature_positions") or []:
                        b, en = fp.get("beg_seq_id"), fp.get("end_seq_id") or fp.get("beg_seq_id")
                        if b:
                            miss.append((b, en))
                            unobs_union.update(range(b, en + 1))
            inst_detail.append({"chain": ch, "unobserved": miss,
                                "n_unobs": sum(e2 - b2 + 1 for b2, e2 in miss)})

        tags = [nm for nm, rx in TAGS if seq and rx.search(seq)]
        ents[ent] = {
            "type": ptype, "desc": desc, "construct_len": len(seq),
            "uniprots": ups, "ref_coverage": cover, "ref_span": ref_span,
            "n_unobserved_union": len(unobs_union),
            "instances": inst_detail, "tags": tags,
            "is_fusion": len(ups) >= 2,
        }
    report[eid] = ents

json.dump(report, open("seqgap.json", "w"), ensure_ascii=False, indent=1)

# ---------------- summary ----------------
print("\n" + "=" * 78)
print("A→B  构建体 vs 全长(蛋白链,有 UniProt 参考的)")
print("=" * 78)
rows = []
for r in recs:
    for ent, d in (report.get(r["id"]) or {}).items():
        if d["type"] != "Protein":
            continue
        up = d["uniprots"][0] if d["uniprots"] else None
        if not up or not up.get("len"):
            continue
        frac = (d["ref_coverage"] or 0) / up["len"] if up["len"] else None
        rows.append((frac, r["no"], r["id"], ent, d, up))
rows.sort()
for frac, no, pid, ent, d, up in rows:
    span = d["ref_span"]
    sp = f"{span[1]}-{span[2]}" if span else "?"
    print(f"#{no:2d} {pid} e{ent}  构建体 {d['construct_len']:4d} aa / 全长 {up['len']:5d} aa"
          f"  = {frac*100:5.1f}%  覆盖 {sp:>12s}  {d['desc'][:38]}")

print("\n" + "=" * 78)
print("B→C  未观测到坐标的残基(构建体里有、但结构里看不到)")
print("=" * 78)
for r in recs:
    for ent, d in (report.get(r["id"]) or {}).items():
        if d["type"] != "Protein" or not d["n_unobserved_union"]:
            continue
        pct = 100.0 * d["n_unobserved_union"] / max(1, d["construct_len"])
        print(f"#{r['no']:2d} {r['id']} e{ent}  未观测 {d['n_unobserved_union']:4d}"
              f" / {d['construct_len']:4d} aa = {pct:5.1f}%   {d['desc'][:40]}")

print("\n" + "=" * 78)
print("需要特别处理:融合构建体 / 含表达标签")
print("=" * 78)
for r in recs:
    for ent, d in (report.get(r["id"]) or {}).items():
        if d["is_fusion"] or d["tags"]:
            what = []
            if d["is_fusion"]:
                what.append("融合(" + "+".join(u["acc"] or "?" for u in d["uniprots"]) + ")")
            if d["tags"]:
                what.append("标签: " + ", ".join(d["tags"]))
            print(f"#{r['no']:2d} {r['id']} e{ent}  {'; '.join(what)}   {d['desc'][:44]}")
