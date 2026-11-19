#!/usr/bin/env python3
"""Build module/nlr_signatures.tsv and module/nlr_classes.tsv.

Reads the signature -> domain-class mapping out of the legacy case_when chain in the 1.x
script, which git still holds:  git show v1.0.3:module/NLRtracker.R > /tmp/NLRtracker.legacy.R
It then adds the curated InterPro additions, enriches every accession from
the InterProScan datadir's entries.json, and writes two data files that drive the R pipeline.

Re-run this when InterPro is updated:

    python3 module/build_signature_table.py \
        --entries  ~/Bioinformatics/interproscan6/data/interpro/111.0/entries.json \
        --legacy   /tmp/NLRtracker.legacy.R \
        --outdir   module

Provenance for every curated addition is in docs/signatures.md.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import OrderedDict

# --------------------------------------------------------------------------------------------
# Class-level metadata.
#
# region        the arbitrated region a class contributes to, or "-" for none
# merge_group   classes sharing a merge_group are merged when their intervals overlap
# precedence    higher wins when a merge collapses several classes into one interval
# drop_after    the class only exists to extend a region; it is removed once merging is done
# simple        the letter used in the "Domain architecture simplified" string
# colour/shape  iTOL and figure styling
# --------------------------------------------------------------------------------------------
CLASSES = [
    # class,             region,   merge_group, prec, drop_after, simple, colour,    shape
    ("NBARC",            "NB-ARC", "nbarc",      40, False, "N", "#FED9A6", "RE"),
    ("NACHT",            "NB-ARC", "nbarc",      30, False, "(N)", "#FED9A6", "RE"),
    # PLOOP keeps its own region. A bare P-loop must not be drawn or extracted as a
    # full NB-ARC, because that overstates the evidence for a degenerate NLR. It still
    # collapses into NBARC through the shared merge group whenever the two overlap.
    ("PLOOP",            "PLOOP",  "nbarc",      20, False, "P", "#FFFFCC", "RE"),
    ("nbarc-motif",      "NB-ARC", "nbarc",      10, True,  "",  "",        ""),
    ("Winged-helix",     "NB-ARC", "nbarc",      10, True,  "",  "",        ""),
    ("nacht-winged-helix", "NB-ARC", "nbarc",    10, True,  "",  "",        ""),

    ("TIR",              "TIR",    "-",           0, False, "T", "#DECBE4", "OC"),
    ("CC",               "CC",     "-",           0, False, "C", "#B3CDE3", "OC"),
    ("R1",               "R1",     "-",           0, False, "B", "#FBB4AE", "HH"),
    ("LRR",              "LRR",    "-",           0, False, "L", "#E5D8BD", "EL"),
    ("CJID",             "C-JID",  "-",           0, False, "J", "#FDDAEC", "HH"),

    # RPW8 / MLKL / CblN share one merge group; RPW8 wins, then MLKL, then CblN.
    ("RPW8",             "RPW8",   "rpw8",       40, False, "R", "#CCEBC5", "OC"),
    ("MLKL",             "MLKL",   "rpw8",       30, False, "M", "#CCEBC5", "OC"),
    ("CblN",             "ID",     "rpw8",       20, False, "O", "#666666", "HH"),

    # Non-plant NLR discriminators. Not regions; they inform Status only.
    ("HET",              "-",      "-",           0, False, "(H)", "#666666", "HH"),
    ("NN",               "-",      "-",           0, False, "(-)", "#666666", "HH"),
    ("FISNA",            "-",      "-",           0, False, "(F)", "#666666", "HH"),
    ("CARD",             "ID",     "id",         40, False, "(C)", "#666666", "HH"),
    ("PYD",              "ID",     "id",         40, False, "(P)", "#666666", "HH"),
    ("BIR",              "ID",     "id",         40, False, "(B)", "#666666", "HH"),

    # Integrated domains. All contribute to the general ID region.
    ("PKin",             "ID",     "id",         40, False, "O", "#666666", "HH"),
    ("TRX",              "ID",     "id",         40, False, "O", "#666666", "HH"),
    ("Exo70",            "ID",     "id",         40, False, "O", "#666666", "HH"),
    ("BED",              "ID",     "id",         40, False, "O", "#666666", "HH"),
    ("TRANSPOSON",       "ID",     "id",         40, False, "O", "#666666", "HH"),
    ("WRKY",             "ID",     "id",         40, False, "O", "#666666", "HH"),
    ("HMA",              "ID",     "id",         40, False, "O", "#666666", "HH"),
    ("JAC",              "ID",     "id",         40, False, "O", "#666666", "HH"),
    ("BRX",              "ID",     "id",         40, False, "O", "#666666", "HH"),
    # Super-structure repeats lose to OTHER when they overlap it (legacy behaviour).
    ("a/b",              "ID",     "id",         10, False, "O", "#666666", "HH"),
    ("ANK",              "ID",     "id",         10, False, "O", "#666666", "HH"),
    ("TPR",              "ID",     "id",         10, False, "O", "#666666", "HH"),
    ("WD40",             "ID",     "id",         10, False, "O", "#666666", "HH"),
    ("ARM",              "ID",     "id",         10, False, "O", "#666666", "HH"),
    ("OTHER",            "ID",     "id",         20, False, "O", "#666666", "HH"),
    ("PLOOP-OTHER",      "ID",     "id",         20, False, "O", "#666666", "HH"),

    # Motif-derived and evidence-only classes. Never regions.
    ("disease",          "-",      "-",           0, True,  "",  "",        ""),
    ("cc-motif",         "-",      "-",           0, True,  "",  "",        ""),
    ("rnbs-d",           "-",      "-",           0, True,  "",  "",        ""),
    ("other-motif",      "-",      "-",           0, True,  "",  "",        ""),
    ("linker-MHD",       "-",      "-",           0, True,  "",  "",        ""),
]

# --------------------------------------------------------------------------------------------
# Curated InterPro additions.
#
# Every row here is justified by evidence recorded in docs/signatures.md:
# a hit count on RefPlantNLR (481 curated plant NLRs) and on a full Nicotiana tabacum
# proteome (75 127 proteins).
#
# subclass_evidence values:
#   CC-NLR / TIR-NLR                 diagnostic, may decide a subclass
#   CC-NLR:supporting                may corroborate, never decides
#   CCG10-NLR:gated                  only fires when no true CC signature is present
#   -                                no subclass role
# --------------------------------------------------------------------------------------------
ADDITIONS = [
    # signature, class,  subclass_evidence,   n_refplantnlr, n_ntab, note
    ("PF23282", "Winged-helix", "TIR-NLR",          99,  166,
     "NB-ARC WHD; DIAGNOSTIC TIR-lineage, purity 100% covers 94%"),
    ("PF23559", "Winged-helix", "CC-NLR",          351,  952,
     "NB-ARC WHD; DIAGNOSTIC CC-lineage, purity 98% covers 99% - best CC marker in the set"),
    ("PF25895", "Winged-helix", "-",                 0,    6,
     "NB-ARC WHD, plant disease resistance; rare in RefPlantNLR, confirmed in proteome scan"),
    ("PF23598", "LRR",          "-",               297, 1064,
     "LRR region and extraction evidence (also hits cell-surface LRR receptors such as FLS2); 87% NLR purity too low for subclass"),
    ("PF25019", "LRR",          "CC-NLR:supporting", 80,   68,
     "NLR-specific LRR; purity 100% covers 23%"),
    ("PF23247", "LRR",          "CCG10-NLR:gated",  21,   18,
     "NLR-specific LRR, named for RPS2; gated CCG10 rule, precision 100% covers 58%"),
    ("PF23286", "LRR",          "TIR-NLR:supporting", 18,  70,
     "NLR-specific LRR; purity 100% covers 17%"),
    ("cd21037", "MLKL",         "-",                 1,   38,
     "MLKL N-terminal 4H bundle; leads the MLKL group, stronger than PF06760"),
    ("SSF47668", "MLKL",        "-",                 0,    0,
     "4H bundle superfamily IPR036537, mammalian-biased; 0 in 75127 plant proteins"),
]

# Subclass rules that cannot be expressed per-signature. Consumed by 04_classify.R.
SUBCLASS_RULES = [
    # name, require_any, forbid_any, subclass, qualifier, note
    ("ccg10_gated", "PF23247", "PF18052,cd14798", "CCG10-NLR", "probable",
     "No CC domain is detectable in any CCG10-NLR; gate gives precision 100%, coverage 11/19"),
]

DB_OF_PREFIX = [
    ("G3DSA", "CATH-Gene3D"), ("SSF", "SUPERFAMILY"), ("PTHR", "PANTHER"),
    ("PF", "Pfam"), ("PS", "PROSITE"), ("SM", "SMART"), ("PR", "PRINTS"),
    ("cd", "CDD"), ("Motif", "MEME"),
]


def guess_db(acc: str) -> str:
    for pre, db in DB_OF_PREFIX:
        if acc.startswith(pre):
            return db
    return "?"


def parse_legacy(path: str) -> "OrderedDict[str, list[str]]":
    """Pull signature -> domain class out of the legacy case_when chain."""
    src = open(path).read()
    pat = re.compile(
        r'Signature\s*%in%\s*c\(\s*((?:"[^"]+"\s*,?\s*)+)\)'
        r'(?:\s*&[^~]*)?\s*~\s*"\(([^)]*)\)"',
        re.S,
    )
    out: "OrderedDict[str, list[str]]" = OrderedDict()
    for m in pat.finditer(src):
        klass = m.group(2)
        for acc in re.findall(r'"([^"]+)"', m.group(1)):
            out.setdefault(klass, [])
            if acc not in out[klass]:
                out[klass].append(acc)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--entries", required=True, help="path to interpro/<ver>/entries.json")
    ap.add_argument("--legacy", required=True, help="reference copy of the legacy NLRtracker.R")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--interpro-version", default="111.0")
    a = ap.parse_args()

    entries = json.load(open(a.entries))
    legacy = parse_legacy(a.legacy)
    known_classes = {c[0] for c in CLASSES}

    # class -> ordered signatures, legacy first then curated additions
    per_class: "OrderedDict[str, list[str]]" = OrderedDict()
    for klass, accs in legacy.items():
        if klass not in known_classes:
            print(f"WARNING: legacy class ({klass}) has no CLASSES entry - skipped",
                  file=sys.stderr)
            continue
        per_class.setdefault(klass, []).extend(accs)

    add_meta = {}
    for acc, klass, sub, n_ref, n_ntab, note in ADDITIONS:
        if klass not in known_classes:
            print(f"ERROR: addition {acc} targets unknown class ({klass})", file=sys.stderr)
            return 1
        per_class.setdefault(klass, [])
        if acc not in per_class[klass]:
            per_class[klass].append(acc)
        add_meta[acc] = (sub, n_ref, n_ntab, note)

    # de-duplicate: first class wins, matching case_when semantics
    seen: dict[str, str] = {}
    dupes = []
    sig_rows = []
    for klass in per_class:
        for acc in per_class[klass]:
            if acc in seen:
                dupes.append((acc, seen[acc], klass))
                continue
            seen[acc] = klass
            e = entries.get(acc)
            db = e["database"] if e else guess_db(acc)
            ipr = (e.get("integrated") or "-") if e else "-"
            short = (e.get("name") or "-") if e else "-"
            desc = (e.get("description") or "-") if e else "-"
            status = "ok" if e else ("internal" if acc == "CJID" else
                                     "motif" if acc.startswith("Motif") else "ABSENT")
            sub, n_ref, n_ntab, note = add_meta.get(acc, ("-", "", "", ""))
            # a fallback name is needed only when the TSV gives neither entry nor description
            needs_fallback = "yes" if (status == "ok" and ipr == "-" and desc == "-") else "no"
            sig_rows.append([acc, db, ipr, short, desc, klass, sub, status,
                             needs_fallback, str(n_ref), str(n_ntab), note])

    out_sig = f"{a.outdir}/nlr_signatures.tsv"
    with open(out_sig, "w") as fh:
        fh.write(f"# NLRtracker signature table - generated by build_signature_table.py\n")
        fh.write(f"# InterPro {a.interpro_version}\n")
        fh.write("# Provenance for curated rows: docs/signatures.md\n")
        fh.write("signature\tdatabase\tipr_entry\tipr_short_name\tipr_description"
                 "\tdomain_class\tsubclass_evidence\tstatus\tneeds_name_fallback"
                 "\tn_refplantnlr\tn_ntab_proteome\tnotes\n")
        for r in sig_rows:
            fh.write("\t".join(x.replace("\t", " ") for x in r) + "\n")

    out_cls = f"{a.outdir}/nlr_classes.tsv"
    with open(out_cls, "w") as fh:
        fh.write("# NLRtracker domain-class metadata - generated by build_signature_table.py\n")
        fh.write("# region: the arbitrated region this class feeds, or - for none\n")
        fh.write("# merge_group + precedence: overlapping classes in one group collapse,"
                 " highest precedence names the result\n")
        fh.write("# drop_after: class exists only to extend a region, removed after merging\n")
        fh.write("domain_class\tregion\tmerge_group\tprecedence\tdrop_after"
                 "\tsimple_code\tcolour\tshape\n")
        for k, region, mg, prec, drop, simple, colour, shape in CLASSES:
            fh.write(f"{k}\t{region}\t{mg}\t{prec}\t{'yes' if drop else 'no'}"
                     f"\t{simple}\t{colour}\t{shape}\n")

    out_rules = f"{a.outdir}/nlr_subclass_rules.tsv"
    with open(out_rules, "w") as fh:
        fh.write("# Subclass rules that need more than one signature - consumed by"
                 " module/R/04_classify.R\n")
        fh.write("rule\trequire_any\tforbid_any\tsubclass\tqualifier\tnotes\n")
        for r in SUBCLASS_RULES:
            fh.write("\t".join(r) + "\n")

    # ---------------------------------------------------------------- report
    print(f"wrote {out_sig}      {len(sig_rows)} signatures")
    print(f"wrote {out_cls}      {len(CLASSES)} classes")
    print(f"wrote {out_rules}   {len(SUBCLASS_RULES)} rules")

    absent = [r for r in sig_rows if r[7] == "ABSENT"]
    fallback = [r for r in sig_rows if r[8] == "yes"]
    print(f"\nstatus: ok={sum(1 for r in sig_rows if r[7]=='ok')}"
          f"  motif={sum(1 for r in sig_rows if r[7]=='motif')}"
          f"  internal={sum(1 for r in sig_rows if r[7]=='internal')}"
          f"  ABSENT={len(absent)}")
    if absent:
        print("  ABSENT from InterPro - these can never match:")
        for r in absent:
            print(f"    {r[0]:20s} ({r[5]})")
    print(f"\nneed a hardcoded name fallback ({len(fallback)}): "
          + ", ".join(r[0] for r in fallback))
    if dupes:
        print(f"\ndropped {len(dupes)} duplicate assignment(s) - first class wins:")
        for acc, first, second in dupes:
            print(f"    {acc:20s} kept ({first})  dropped ({second})")
    classes_used = {r[5] for r in sig_rows}
    unused = [c[0] for c in CLASSES if c[0] not in classes_used]
    if unused:
        print(f"\nclasses with no signatures (motif-only or evidence-only): "
              + ", ".join(unused))
    return 0


if __name__ == "__main__":
    sys.exit(main())
