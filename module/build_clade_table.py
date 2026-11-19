#!/usr/bin/env python3
"""Generate module/nlr_clades.tsv from the RefPlantNLR Table S1 subclass column.

Every anchor is a verified RefPlantNLR identifier, so the clade definitions are traceable to
published assignments rather than hand-typed. Re-run if Table S1 is updated.

    python3 module/build_clade_table.py \
        --table-s1 /path/to/Table_S1_RefPlantNLR.xlsx \
        --tnp-fasta /path/to/outgroup.fasta \
        --output module/nlr_clades.tsv

Table S1 encodes the clade in the `Subclass` column as a ';'-separated path, e.g.
`NLR;CNL;NRC;G08//CNL14`. Decoded:

    CNL -> CC-NLR    GT -> TIR-NLR    G10 -> CCG10-NLR    GR -> CCR-NLR (RPW8-type CC)

`CC-01` and `CC-02` are deliberately absent: they are defined by the tree as the two daughters
of the CC MRCA, then labelled by which one contains the NRC anchors and which the MLA1 anchors.
Hardcoding their membership would defeat the point of building a phylogeny.
"""
from __future__ import annotations

import argparse
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

# The NRC helper clade is exactly Table S1's `NRC;G08//CNL14` plus NRC0. Everything else
# carrying an NRC label is a sensor. This is a data-grounded split, not a typed list.
NRC_HELPER_TOKENS = {"G08", "CNL14"}
NRC_HELPER_EXTRA = {"NRC0"}


def read_table_s1(path: str) -> dict[str, str]:
    z = zipfile.ZipFile(path)
    shared: list[str] = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(f"{NS}si"):
            shared.append("".join(t.text or "" for t in si.iter(f"{NS}t")))

    def cidx(ref: str) -> int:
        n = 0
        for ch in re.match(r"([A-Z]+)", ref).group(1):
            n = n * 26 + ord(ch) - 64
        return n - 1

    rows = []
    for row in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).iter(f"{NS}row"):
        cells = {}
        for c in row.findall(f"{NS}c"):
            if not c.get("r"):
                continue
            v = c.find(f"{NS}v")
            if c.get("t") == "s" and v is not None:
                val = shared[int(v.text)]
            else:
                val = v.text if v is not None else ""
            cells[cidx(c.get("r"))] = (val or "").strip()
        if cells:
            rows.append([str(cells.get(i, "")) for i in range(max(cells) + 1)])

    head = rows[0]
    i_id, i_sub = head.index("ID"), head.index("Subclass")
    out = {}
    for r in rows[1:]:
        if len(r) > max(i_id, i_sub) and r[i_id]:
            out[r[i_id]] = r[i_sub]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table-s1", required=True)
    ap.add_argument("--tnp-fasta", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    subs = read_table_s1(a.table_s1)
    tnp = [l[1:].split()[0] for l in open(a.tnp_fasta) if l.startswith(">")]
    if not tnp:
        print(f"ERROR: no sequences in {a.tnp_fasta}", file=sys.stderr)
        return 1

    def tokens(s: str) -> set[str]:
        out = set()
        for part in s.split(";"):
            for q in part.split("//"):
                if q.strip():
                    out.add(q.strip())
        return out

    by_token: dict[str, list[str]] = defaultdict(list)
    for pid, sub in subs.items():
        for t in tokens(sub):
            by_token[t].append(pid)

    nrc_all = sorted(by_token.get("NRC", []) + by_token.get("NRC-S", []))
    nrc_helpers = sorted(
        [p for p in nrc_all if tokens(subs[p]) & NRC_HELPER_TOKENS] +
        [p for p in nrc_all if p in NRC_HELPER_EXTRA])
    nrc_sensors = sorted(set(nrc_all) - set(nrc_helpers))

    # Curated CC structure, per the user 2026-08-05. This supersedes an earlier guess that
    # put NRC in CC-01 and MLA in CC-02: the tree placed NRC and MLA in the SAME daughter of
    # the CC node, and that is correct -- they are both in CC clade 2. CC clade 1 is a
    # separate group (Ptr1 / alpha-WOLF / Pm2 / I2 / R3a / Xa1).
    #
    # Some entries are families in RefPlantNLR rather than single ids, so they are expanded
    # by prefix: "alpha-WOLF" is alpha-WOLF10/11/12, "Pm2" is Pm2a/2d/2e/21, and Ptr1 has
    # orthologues CaPtr1 / NbPtr1 / NiattPtr1 / NitabPtr1.
    def expand(names, prefixes=()):
        out = set()
        for n in names:
            if n in subs:
                out.add(n)
        for pre in prefixes:
            for p in subs:
                if not p.startswith(pre):
                    continue
                # An allelic series continues with a letter or a separator, not another
                # digit: "Pm2a".."Pm2j" are Pm2 alleles but "Pm21" is a different wheat gene
                # from Haynaldia villosa and must not be absorbed.
                rest = p[len(pre):]
                if rest[:1].isdigit() and pre[-1:].isdigit():
                    continue
                out.add(p)
        return sorted(out)

    cc01 = sorted(
        set(expand(["I2", "R3a", "Xa1"], ("Ptr1", "CaPtr1", "NbPtr1", "NiattPtr1",
                                          "NitabPtr1", "alpha-WOLF", "Pm2")))
        | set(by_token.get("CNL8", []))
        | set(by_token.get("CNL6", []))
        | set(by_token.get("ZfBED", []))
    )
    mla = sorted(set(by_token.get("MLA1", [])) |
                 set(expand(["Sr35", "Sr50", "Sr33", "DES1"])))
    zar1 = expand([], ("ZAR1", "NbZAR1"))
    sensor1 = expand(["Bs2", "Rx", "Gpa2", "Rx2", "Rx4", "Rpa1", "CNL-11990"],
                     ("Rpi-amr1", "Rpi-amr3", "Rpi-nig3"))
    sensor2 = expand(["CaMi", "Mi-1", "Rpi-blb2", "Hero-A", "R8", "Sw5-b",
                      "NbPRFa", "NbPRFb", "Prf", "R1", "NtTPN1"])
    rga5 = sorted(set(by_token.get("RGA5", [])) | set(by_token.get("MIC1", [])))
    cc02 = sorted(
        set(nrc_all) | set(mla) | set(zar1)
        | set(by_token.get("Pik-1", []))
        | set(by_token.get("Pik-2", []))
        | set(by_token.get("RGA4", []))
        | set(rga5)
        | set(by_token.get("RPP13", []))
        | set(by_token.get("CNL4", []))
    )

    # name, parent, role, anchors, note
    rows = [
        ("TNP", "root", "outgroup", sorted(tnp),
         "divergent TIR-NLRs supplied separately; roots the tree"),
        ("TIR-NLR", "root", "anchored", sorted(by_token.get("GT", [])), "Table S1 GT"),
        ("CCG10-NLR", "root", "anchored", sorted(by_token.get("G10", [])),
         "Table S1 G10; no CC domain is detectable in any of these"),
        ("CCR-NLR", "root", "anchored", sorted(by_token.get("GR", [])),
         "Table S1 GR, the RPW8-type CC helpers"),
        ("CCR-NRG1", "CCR-NLR", "anchored",
         sorted(p for p in by_token.get("GR", []) if "NRG1" in p or p == "DAR5"),
         "NRG1-type"),
        ("CCR-ADR1", "CCR-NLR", "anchored",
         sorted(p for p in by_token.get("GR", []) if "NRG1" not in p and p != "DAR5"),
         "ADR1-type"),
        ("CC-NLR", "root", "anchored", sorted(by_token.get("CNL", [])), "Table S1 CNL"),
        ("CC-01", "CC-NLR", "anchored", cc01,
         "CC clade 1: Ptr1 / alpha-WOLF / Pm2 / I2 / R3a / Xa1"),
        ("CC-02", "CC-NLR", "anchored", cc02,
         "CC clade 2: NRC and MLA together, ZAR1 typically basal"),
        ("CC-02-ZAR1", "CC-02", "anchored", zar1,
         "ZAR1 / NbZAR1, usually but not always basal in CC clade 2"),
        ("CC-02-MLA", "CC-02", "anchored", mla,
         "Table S1 MLA1 plus Sr35, Sr50, Sr33, DES1"),
        ("CC-02-NRCh", "CC-02", "anchored", nrc_helpers,
         "Table S1 NRC;G08//CNL14 plus NRC0 -- the canonical NRC helpers"),
        ("CC-02-NRC-sensor-1", "CC-02", "anchored", sensor1,
         "Bs2, Rx, Gpa2, Rx2, Rx4, Rpi-amr1*, Rpi-amr3*, Rpi-nig3*, Rpa1, CNL-11990"),
        ("CC-02-NRC-sensor-2", "CC-02", "anchored", sensor2,
         "CaMi, Mi-1, Rpi-blb2, Hero-A, R8, Sw5-b, NbPRFa/b, Prf, R1, NtTPN1"),
    ]

    # finer named subclades nested under CC-01 or CC-02 (MIC1 merged into RGA5; CNL16
    # omitted because five of its six Table S1 anchors lack RefPlantNLR sequences).
    for name, parent, ids, note in [
        ("Pik-1", "CC-02", sorted(by_token.get("Pik-1", [])), "Pik-1 sensor"),
        ("Pik-2", "CC-02", sorted(by_token.get("Pik-2", [])), "Pik-2 helper"),
        ("RGA4",  "CC-02", sorted(by_token.get("RGA4", [])),  "RGA4 helper"),
        ("RGA5",  "CC-02", rga5,                              "RGA5 sensor"),
        ("RPP13", "CC-02", sorted(by_token.get("RPP13", [])), "RPP13"),
        ("ZfBED", "CC-01", sorted(by_token.get("ZfBED", [])), "integrated ZfBED"),
        ("CNL8",  "CC-01", sorted(by_token.get("CNL8", [])),  "G04 / I2, L-locus"),
        ("CNL6",  "CC-01", sorted(by_token.get("CNL6", [])),  "G07"),
        ("CNL4",  "CC-02", sorted(by_token.get("CNL4", [])),  "G11"),
    ]:
        if len(ids) >= 2:
            rows.append((name, parent, "anchored", ids, note))

    with open(a.output, "w") as fh:
        fh.write("# NLRtracker clade definitions - generated by build_clade_table.py\n")
        fh.write("# Anchors are RefPlantNLR identifiers from Table S1's Subclass column.\n")
        fh.write("# role: anchored = clade is the MRCA of its anchors\n")
        fh.write("#       tree_defined = derived from the tree, see notes\n")
        fh.write("#       outgroup = used to root the tree\n")
        fh.write("clade\tparent\trole\tn_anchors\tanchors\tnotes\n")
        for name, parent, role, ids, note in rows:
            fh.write(f"{name}\t{parent}\t{role}\t{len(ids)}\t{','.join(ids)}\t{note}\n")

    print(f"wrote {a.output}")
    print(f"{'clade':20s} {'parent':10s} {'role':13s} anchors")
    for name, parent, role, ids, _ in rows:
        print(f"  {name:18s} {parent:10s} {role:13s} {len(ids)}")
    print(f"\nNRC split: {len(nrc_helpers)} helpers, {len(nrc_sensors)} sensors")
    print(f"  helpers: {', '.join(nrc_helpers)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
