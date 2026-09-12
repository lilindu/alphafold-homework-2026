"""Switch protein inputs from the construct sequence (B) to the full-length
UniProt sequence (A), which is what a researcher actually has in practice.

Four things must be recomputed, not just substituted:

1. PTM POSITIONS. Modification positions are currently in construct numbering
   (entity label_seq_id). Under full-length input they must be remapped to
   UniProt numbering via the deposited alignment, or the modification lands on
   the wrong residue. This is the highest-risk part of the change.

2. TOKENS. Full-length is longer -- sometimes far longer (mesothelin 17 -> 622
   aa; West Nile envelope 102 -> 3433 aa). The 5,000-token ceiling may now be
   exceeded, which would make a target unusable.

3. WHICH CHAINS CAN BE CONVERTED. A chain gets a full-length sequence only if it
   maps to exactly one UniProt entry. Synthetic peptides, de novo designed
   proteins and fusion chimeras have no single natural full-length form, so they
   keep the deposited sequence -- for them the construct IS the whole molecule.

4. STANDARD RESIDUES. UniProt canonical sequences can contain U (selenocysteine)
   or other non-standard letters the Server rejects; must re-check.

Expression tags disappear as a side effect: they are not part of the UniProt
sequence. DNA/RNA chains are unchanged -- an oligo used in crystallography has no
"full-length" form.
"""
import json
import os

OUT = "/Users/lilindu/alphafold-homework-2026"

recs = json.load(open(os.path.join(OUT, "targets41.json")))
raw = json.load(open("seqgap_raw.json"))

STD = set("ACDEFGHIKLMNPQRSTVWY")
TOKEN_LIMIT = 5000


def entity_info(eid):
    """entity_id -> alignment + uniprot payload."""
    out = {}
    e = raw.get(eid)
    if not e:
        return out
    for pe in e.get("polymer_entities") or []:
        ent = pe["rcsb_id"].split("_")[-1]
        ups = []
        for u in pe.get("uniprots") or []:
            seq = ((u.get("rcsb_uniprot_protein") or {}).get("sequence") or "")
            ups.append({"acc": u.get("rcsb_id"), "seq": seq})
        regions = []
        acc = None
        for al in pe.get("rcsb_polymer_entity_align") or []:
            if (al.get("reference_database_name") or "").upper().startswith("UNIPROT"):
                acc = al.get("reference_database_accession")
                for x in al.get("aligned_regions") or []:
                    if all(x.get(k) is not None
                           for k in ("entity_beg_seq_id", "ref_beg_seq_id", "length")):
                        regions.append((x["entity_beg_seq_id"], x["ref_beg_seq_id"],
                                        x["length"]))
        out[ent] = {"uniprots": ups, "align_acc": acc, "regions": regions}
    return out


def remap(pos, regions):
    """construct position -> reference position."""
    for ebeg, rbeg, length in regions:
        if ebeg <= pos < ebeg + length:
            return rbeg + (pos - ebeg)
    return None


converted, kept, problems = [], [], []
new_recs = []

for r in recs:
    info = entity_info(r["id"])
    chains, notes = [], []
    ok = True

    for c in r["chains"]:
        c = dict(c)
        if c["type"] != "Protein":
            chains.append(c)
            continue

        d = info.get(c["entity"]) or {}
        ups = [u for u in d.get("uniprots") or [] if u.get("seq")]
        regions = d.get("regions") or []

        if len(ups) == 1 and regions:
            full = ups[0]["seq"].upper()
            bad = set(full) - STD
            if bad:
                problems.append((r["no"], r["id"], c["entity"],
                                 f"UniProt 序列含非标准残基 {''.join(sorted(bad))},保留构建体"))
                c["source"] = "construct"
                c["reason"] = "全长序列含 Server 不支持的残基"
                chains.append(c)
                kept.append((r["id"], c["entity"]))
                continue
            c["construct_seq"] = c["seq"]
            c["construct_len"] = c["len"]
            c["seq"] = full
            c["len"] = len(full)
            c["source"] = "uniprot"
            c["uniprot"] = ups[0]["acc"]
            c["regions"] = regions
            converted.append((r["id"], c["entity"], c["construct_len"], c["len"]))
        else:
            why = ("无 UniProt 对应(合成肽 / 人工设计蛋白)" if not ups
                   else f"融合构建体({len(ups)} 个 UniProt),无单一全长形式"
                   if len(ups) > 1 else "缺少比对信息")
            c["source"] = "construct"
            c["reason"] = why
            kept.append((r["id"], c["entity"]))
        chains.append(c)

    # ---- remap PTM positions ----
    new_ptms = []
    for p in r["ptms"]:
        ch = next((c for c in chains if c["entity"] == p["entity"]), None)
        if ch is None:
            continue
        if ch.get("source") == "uniprot":
            np_ = remap(p["pos"], ch["regions"])
            if np_ is None or not (1 <= np_ <= ch["len"]):
                problems.append((r["no"], r["id"], p["entity"],
                                 f"修饰 {p['code']}@{p['pos']} 无法映射到全长编号"))
                ok = False
                continue
            parent = {"SEP": "S", "TPO": "T", "PTR": "Y", "ALY": "K", "MLY": "K",
                      "M3L": "K", "MLZ": "K", "CIR": "R", "HIP": "H", "NEP": "H",
                      "HYP": "P", "KCR": "K"}.get(p["code"])
            aa = ch["seq"][np_ - 1]
            if parent and aa != parent:
                problems.append((r["no"], r["id"], p["entity"],
                                 f"修饰 {p['code']} 映射到全长第 {np_} 位是 {aa},"
                                 f"应为 {parent} —— 编号对不上"))
                ok = False
                continue
            q = dict(p)
            q["pos_construct"] = p["pos"]
            q["pos"] = np_
            q["parent_aa"] = aa
            new_ptms.append(q)
        else:
            new_ptms.append(dict(p))

    tokens = sum(c["len"] * c["copies"] for c in chains)
    tokens += sum((l["atoms"] or 0) * l["count"] for l in r["ligands"])
    tokens += sum(i["count"] for i in r["ions"])

    q = dict(r)
    q["chains"] = chains
    q["ptms"] = new_ptms
    q["tokens_construct"] = r["tokens"]
    q["tokens"] = tokens
    q["input_mode"] = "full-length"
    q["ptm_ok"] = ok
    new_recs.append(q)

print(f"protein chains converted to full length: {len(converted)}")
print(f"chains kept as deposited:               {len(kept)}")
print()
print("=== token change ===")
over = []
for q in sorted(new_recs, key=lambda x: -x["tokens"]):
    d = q["tokens"] - q["tokens_construct"]
    flag = ""
    if q["tokens"] > TOKEN_LIMIT:
        flag = "   *** 超过 5000 上限 ***"
        over.append(q)
    if d or flag:
        print(f"#{q['no']:2d} {q['id']}  {q['tokens_construct']:5d} -> {q['tokens']:5d}"
              f"  ({d:+6d}){flag}")

print(f"\ntargets now over the 5,000-token limit: {len(over)}")
for q in over:
    print(f"   #{q['no']} {q['id']} = {q['tokens']}")
    for c in q["chains"]:
        if c["type"] == "Protein" and c.get("source") == "uniprot":
            print(f"       {c['uniprot']} {c['construct_len']} -> {c['len']} aa"
                  f" x{c['copies']}   {c['desc'][:44]}")

if problems:
    print(f"\n=== problems ({len(problems)}) ===")
    for no, pid, ent, msg in problems:
        print(f"   #{no} {pid} e{ent}: {msg}")

json.dump(new_recs, open("targets41_full.json", "w"), ensure_ascii=False, indent=1)
print("\nwrote targets41_full.json")
