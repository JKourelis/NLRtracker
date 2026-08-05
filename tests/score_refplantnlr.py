#!/usr/bin/env python3
"""Score NLRtracker.tsv against the published RefPlantNLR subclasses.

    score_refplantnlr.py <NLRtracker.tsv> <Table_S1.xlsx>

Prints three numbers: exact, ambiguous-but-contains-the-truth, wrong.

Table S1's `Subclass` column encodes the clade as a ';'-separated path, e.g.
`NLR;CNL;NRC;G08//CNL14`. The second token is the clade:

    CNL -> CC-NLR    GT -> TIR-NLR    G10 -> CCG10-NLR    GR -> CCR-NLR

`GR` maps to CCR-NLR, not to `RPW8`: NLRtracker reserves the bare name `RPW8` for
NLR-associated RPW8 proteins that lack an NB-ARC.

An ambiguous call such as "CC-NLR or CCG10-NLR" is counted separately from a wrong one when it
contains the published answer. Reporting uncertainty honestly is not the same as being wrong,
and collapsing the two would hide the distinction the support tiers exist to make.
"""
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
LABEL = {"CNL": "CC-NLR", "GT": "TIR-NLR", "G10": "CCG10-NLR", "GR": "CCR-NLR"}


def read_table_s1(path):
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(f"{NS}si"):
            shared.append("".join(t.text or "" for t in si.iter(f"{NS}t")))

    def cidx(ref):
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
            val = shared[int(v.text)] if c.get("t") == "s" and v is not None else (
                v.text if v is not None else "")
            cells[cidx(c.get("r"))] = (val or "").strip()
        if cells:
            rows.append([str(cells.get(i, "")) for i in range(max(cells) + 1)])

    head = rows[0]
    i_id, i_sub = head.index("ID"), head.index("Subclass")
    out = {}
    for r in rows[1:]:
        if len(r) > max(i_id, i_sub) and r[i_id]:
            parts = [p for p in r[i_sub].split(";") if p]
            if len(parts) > 1:
                # Normalise the identifier exactly as the pipeline does. InterProScan strips
                # '/' from sequence names, so Table S1's "Bph1/9-1" has to become "Bph19-1"
                # to match. Without this the four Bph1/9 alleles never matched anything and
                # were silently dropped from the score, so the totals covered 477 proteins
                # while appearing to cover 481.
                out[r[i_id].replace("/", "")] = parts[1]
    return out


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    truth = read_table_s1(sys.argv[2])

    pred = {}
    with open(sys.argv[1]) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_id, i_sc = head.index("seqname"), head.index("Subclass (putative)")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(i_id, i_sc):
                pred[f[i_id]] = f[i_sc]

    res = Counter()
    for pid, call in pred.items():
        clade = truth.get(pid)
        if not clade or clade not in LABEL:
            continue
        want = LABEL[clade]
        got = call.replace(" (probable)", "")
        if got == want:
            res["exact"] += 1
        elif " or " in got and want.replace("-NLR", "") in [
                q.replace("-NLR", "") for q in got.split(" or ")]:
            res["ambiguous"] += 1
        else:
            res["wrong"] += 1

    print(res["exact"], res["ambiguous"], res["wrong"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
