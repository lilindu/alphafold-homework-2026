"""Does the modified residue actually sit at a protein-protein interface?

Titles are suggestive but not evidence, so this measures contacts from the
deposited coordinates:

  * download the mmCIF, parse the _atom_site loop directly (Biopython's parser
    is fussy about modified residues appearing as HETATM inside a polymer)
  * locate every Server-supported modified residue
  * compute the minimum heavy-atom distance from that residue to any atom of a
    DIFFERENT polymer chain
  * < 4.0 A  -> the PTM is in direct contact across the interface
    < 5.0 A  -> peripheral
    otherwise the PTM is buried inside its own chain and is NOT an interface case

Also records how many of the PTM's own contacts are to the partner chain, which
separates "the modification IS the binding determinant" (many contacts, e.g. a
phosphate held by an arginine cage) from "it happens to be nearby".
"""
import gzip
import json
import os
import urllib.request
import numpy as np

OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
CACHE = "/tmp/afhw/cif"
os.makedirs(CACHE, exist_ok=True)

SERVER_PTM = {"SEP", "TPO", "PTR", "NEP", "HIP", "ALY", "MLY", "M3L", "MLZ",
              "2MR", "AGM", "MCS", "HYP", "HY3", "LYZ", "AHB", "P1L", "SNN",
              "SNC", "TRF", "KCR", "CIR", "YHA"}
AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE",
} | SERVER_PTM
NUC3 = {"DA", "DT", "DG", "DC", "A", "U", "G", "C", "DU", "PSU", "5MC"}


def fetch(pdb):
    p = os.path.join(CACHE, f"{pdb}.cif.gz")
    if not os.path.exists(p):
        url = f"https://files.rcsb.org/download/{pdb}.cif.gz"
        data = OP.open(url, timeout=180).read()
        with open(p, "wb") as f:
            f.write(data)
    with gzip.open(p, "rt", errors="replace") as f:
        return f.read()


def parse_atoms(text):
    """Yield dicts for the _atom_site loop."""
    lines = text.splitlines()
    i, n = 0, len(lines)
    while i < n:
        if lines[i].strip() == "loop_":
            j, cols = i + 1, []
            while j < n and lines[j].lstrip().startswith("_"):
                cols.append(lines[j].strip())
                j += 1
            if cols and cols[0].startswith("_atom_site."):
                names = [c.split(".", 1)[1] for c in cols]
                idx = {k: t for t, k in enumerate(names)}
                rows = []
                while j < n:
                    ln = lines[j]
                    s = ln.strip()
                    if s.startswith("#") or s == "loop_" or s.startswith("_"):
                        break
                    if s:
                        parts = s.split()
                        if len(parts) >= len(names):
                            rows.append(parts)
                    j += 1
                return idx, rows
            i = j
        else:
            i += 1
    return None, []


def analyse(pdb):
    idx, rows = parse_atoms(fetch(pdb))
    if not idx:
        return None
    need = ["group_PDB", "label_atom_id", "label_comp_id", "label_asym_id",
            "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "auth_asym_id",
            "auth_seq_id", "pdbx_PDB_model_num", "type_symbol"]
    for k in need:
        if k not in idx:
            return None

    coords, comp, asym, seqid, elem, authseq = [], [], [], [], [], []
    for r in rows:
        if r[idx["pdbx_PDB_model_num"]] not in ("1", "."):
            continue
        if r[idx["type_symbol"]] == "H":
            continue
        c = r[idx["label_comp_id"]]
        if c in ("HOH", "DOD"):
            continue
        try:
            xyz = (float(r[idx["Cartn_x"]]), float(r[idx["Cartn_y"]]),
                   float(r[idx["Cartn_z"]]))
        except ValueError:
            continue
        coords.append(xyz)
        comp.append(c)
        asym.append(r[idx["label_asym_id"]])
        seqid.append(r[idx["label_seq_id"]])
        authseq.append(r[idx["auth_seq_id"]])
        elem.append(r[idx["type_symbol"]])
    if not coords:
        return None
    X = np.asarray(coords)
    comp = np.asarray(comp)
    asym = np.asarray(asym)
    authseq = np.asarray(authseq)

    # polymer chains only (a chain counts as polymer if it holds >=4 residues)
    poly = {}
    for a in set(asym.tolist()):
        m = asym == a
        kinds = set(comp[m].tolist())
        res = set(zip(comp[m].tolist(), authseq[m].tolist()))
        if len(res) >= 4 and (kinds & AA3 or kinds & NUC3):
            poly[a] = "protein" if (kinds & AA3) else "nucleic"

    out = []
    ptm_mask = np.isin(comp, list(SERVER_PTM))
    for code, aseq, a in sorted(set(zip(comp[ptm_mask].tolist(),
                                        authseq[ptm_mask].tolist(),
                                        asym[ptm_mask].tolist()))):
        if a not in poly:
            continue
        sel = (comp == code) & (authseq == aseq) & (asym == a)
        P = X[sel]
        best = {"code": code, "chain": a, "auth_seq": aseq,
                "min_dist": None, "partner": None, "n_close": 0,
                "partners": {}}
        for b, kind in poly.items():
            if b == a:
                continue
            Q = X[asym == b]
            if len(Q) == 0:
                continue
            d = np.sqrt(((P[:, None, :] - Q[None, :, :]) ** 2).sum(-1))
            mn = float(d.min())
            nclose = int((d < 4.0).sum())
            best["partners"][b] = {"kind": kind, "min": round(mn, 2), "n<4A": nclose}
            if best["min_dist"] is None or mn < best["min_dist"]:
                best["min_dist"] = round(mn, 2)
                best["partner"] = b
                best["n_close"] = nclose
        out.append(best)
    return {"pdb": pdb, "chains": poly, "ptms": out}


if __name__ == "__main__":
    valid = json.load(open("census_valid.json"))
    cands = [k for k, r in valid.items()
             if r["ptms"] and sum(1 for c in r["chains"] if c["type"] == "Protein") >= 2]
    print("PTM candidates with >=2 protein entities:", len(cands), flush=True)
    res = {}
    if os.path.exists("iface.json"):
        res = json.load(open("iface.json"))
    for k in cands:
        if k in res:
            continue
        try:
            r = analyse(k)
            res[k] = r
            if r:
                for p in r["ptms"]:
                    print(f"  {k} {p['code']}{p['auth_seq']} chain{p['chain']} "
                          f"min={p['min_dist']} n<4A={p['n_close']} -> {p['partner']}",
                          flush=True)
        except Exception as ex:
            print("  ERR", k, str(ex)[:90], flush=True)
            res[k] = None
        json.dump(res, open("iface.json", "w"))
    print("analysed:", len(res))
