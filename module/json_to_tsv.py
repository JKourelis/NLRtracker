#!/usr/bin/env python3
"""Convert InterProScan 6 JSON to the 15-column InterProScan 6 TSV.

Lets NLRtracker consume any IPS6 JSON, including output produced by other wrappers or by an
earlier InterPro release, without re-running InterProScan.

    python3 module/json_to_tsv.py --input proteins.fa.json --output proteins.tsv

Streams the file rather than loading it: real proteome JSON is 400-700 MB pretty-printed
(20M+ lines), and json.load() on that costs several GB of Python objects for no reason.

Column semantics follow the real 6.0.1 writer, which differs from the published docs in two
places (see docs/interproscan6-notes.md section 1):

  * column 10 is the match SOURCE, not the literal "T" the docs claim
  * column 15 is "-" when empty, not "0"

One row per (sequence id, match, location). A single `results` entry can carry several xrefs,
because InterProScan deduplicates identical sequences and lists every identifier that shared
the MD5 -- each one must get its own rows or downstream per-protein counts are wrong.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date

# Which quantity column 9 carries, keyed by lowercased library name. Mirrors
# ProcessOutputTSV.formatLine in the InterProScan 6 source.
SCORE_FROM_MATCH_EVALUE = {"cdd", "prints"}
SCORE_FROM_LOCATION_SCORE = {"hamap", "interpro-n", "prosite profiles", "signalp",
                             "signalp-euk", "signalp-prok"}
SCORE_NONE = {"coils", "mobidb-lite", "phobius", "prosite patterns", "deeptmhmm", "tmbed",
              "antifam"}


def stream_results(path: str):
    """Yield each object in the top-level "results" array, one at a time.

    Tracks brace depth outside string literals. Handles escapes so a brace inside a
    sequence description cannot desynchronise the parser.
    """
    with open(path, "r") as fh:
        # advance to the opening bracket of "results"
        buf = ""
        while True:
            chunk = fh.read(65536)
            if not chunk:
                raise ValueError(f'no "results" array found in {path}')
            buf += chunk
            i = buf.find('"results"')
            if i != -1:
                j = buf.find("[", i)
                if j != -1:
                    rest = buf[j + 1:]
                    break
                # keep only the tail we might still need
                buf = buf[i:]
            elif len(buf) > 1 << 20:
                buf = buf[-4096:]

        depth = 0
        in_str = False
        esc = False
        cur: list[str] = []
        pending = rest
        while True:
            if not pending:
                pending = fh.read(1 << 20)
                if not pending:
                    break
            for k, ch in enumerate(pending):
                if depth > 0:
                    cur.append(ch)
                if esc:
                    esc = False
                    continue
                if ch == "\\":
                    if in_str:
                        esc = True
                    continue
                if ch == '"':
                    in_str = not in_str
                    continue
                if in_str:
                    continue
                if ch == "{":
                    if depth == 0:
                        cur = ["{"]
                    depth += 1
                elif ch == "}":
                    depth -= 1
                    if depth == 0:
                        # parse_float/parse_int=str keeps every number as its original
                        # JSON text. InterProScan writes Java's formatting -- "1.0E-16",
                        # not "1e-16" -- and round-tripping through a Python float
                        # silently rewrites the mantissa and lowercases the exponent, so
                        # the output stops matching a real InterProScan TSV byte for byte.
                        yield json.loads("".join(cur), parse_float=str, parse_int=str)
                        cur = []
                elif ch == "]" and depth == 0:
                    return
            pending = ""


def score_for(library: str, match: dict, loc: dict):
    lib = (library or "").lower()
    if lib in SCORE_NONE:
        return "-"
    if lib in SCORE_FROM_MATCH_EVALUE:
        # The Groovy writer reads match.evalue here, but the JSON serialisation omits a
        # match-level evalue for CDD and carries it on the location instead (verified:
        # ADR1/cd21037 has no match.evalue, and its location evalue 1.02279E-7 is exactly
        # what the real TSV reports). PRINTS does have it at match level. Fall back so
        # both are correct.
        v = match.get("evalue")
        if v is None:
            v = loc.get("evalue")
    elif lib in SCORE_FROM_LOCATION_SCORE:
        v = loc.get("score")
    else:
        v = loc.get("evalue")
    return "-" if v is None else v


def rows_for_protein(entry: dict, run_date: str, goterms: bool, pathways: bool):
    seq = (entry.get("sequence") or "").strip()
    md5 = entry.get("md5", "")
    seqlen = len(seq)
    ids = [x.get("id") for x in entry.get("xref") or [] if x.get("id")]
    if not ids:
        return
    for match in entry.get("matches") or []:
        sig = match.get("signature") or {}
        lib = ((sig.get("signatureLibraryRelease") or {}).get("library")) or "-"
        acc = sig.get("accession") or "-"
        desc = sig.get("description") or "-"
        entry_obj = sig.get("entry")  # null for unintegrated signatures
        ipr_acc = (entry_obj or {}).get("accession") or "-"
        ipr_desc = (entry_obj or {}).get("description") or "-"
        src = match.get("source") or lib

        go = "-"
        if goterms:
            terms = [f"{t['id']}(InterPro)" for t in (entry_obj or {}).get("goXRefs") or []]
            terms += [f"{t['id']}(PANTHER)"
                      for t in (match.get("treegrafter") or {}).get("goXRefs") or []]
            go = "|".join(dict.fromkeys(terms)) or "-"
        pw = "-"
        if pathways:
            pws = [f"{p.get('databaseName')}:{p.get('id')}"
                   for p in (entry_obj or {}).get("pathwayXRefs") or []]
            pw = "|".join(pws) or "-"

        for loc in match.get("locations") or []:
            sc = score_for(lib, match, loc)
            for sid in ids:
                yield (sid, md5, seqlen, lib, acc, desc, loc.get("start"), loc.get("end"),
                       sc, src, run_date, ipr_acc, ipr_desc, go, pw)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="InterProScan 6 JSON")
    ap.add_argument("--output", required=True, help="TSV to write, or - for stdout")
    ap.add_argument("--fasta", default=None,
                    help="also write the matching protein FASTA. The JSON carries every "
                         "sequence, so extracting it here guarantees the identifiers match "
                         "the TSV exactly -- the original input FASTA often does not, "
                         "because wrappers rewrite headers before scanning.")
    ap.add_argument("--goterms", action="store_true")
    ap.add_argument("--pathways", action="store_true")
    ap.add_argument("--date", default=None,
                    help="dd-MM-yyyy for column 11; defaults to today")
    a = ap.parse_args()

    run_date = a.date or date.today().strftime("%d-%m-%Y")

    # Provenance: the release the annotations were actually computed against.
    ips_ver = ipr_ver = "unknown"
    with open(a.input) as fh:
        head = fh.read(4096)
    for key, name in (("interproscan-version", "ips"), ("interpro-version", "ipr")):
        i = head.find(f'"{key}"')
        if i != -1:
            seg = head[i:i + 120]
            parts = seg.split('"')
            if len(parts) > 3:
                if name == "ips":
                    ips_ver = parts[3]
                else:
                    ipr_ver = parts[3]

    out = sys.stdout if a.output == "-" else open(a.output, "w")
    fa = open(a.fasta, "w") if a.fasta else None
    n_prot = n_row = n_seq = 0
    libs: dict[str, int] = {}
    try:
        for entry in stream_results(a.input):
            n_prot += 1
            if fa is not None:
                seq = (entry.get("sequence") or "").strip()
                for x in entry.get("xref") or []:
                    if x.get("id") and seq:
                        fa.write(f">{x['id']}\n{seq}\n")
                        n_seq += 1
            for r in rows_for_protein(entry, run_date, a.goterms, a.pathways):
                out.write("\t".join("" if v is None else str(v) for v in r) + "\n")
                n_row += 1
                libs[r[3]] = libs.get(r[3], 0) + 1
            if n_prot % 20000 == 0:
                print(f"  ... {n_prot:,} proteins, {n_row:,} rows", file=sys.stderr)
    finally:
        if out is not sys.stdout:
            out.close()
        if fa is not None:
            fa.close()

    print(f"InterProScan {ips_ver} / InterPro {ipr_ver}", file=sys.stderr)
    print(f"{n_prot:,} protein records -> {n_row:,} TSV rows", file=sys.stderr)
    if fa is not None:
        print(f"{n_seq:,} sequences -> {a.fasta}", file=sys.stderr)
    for lib, n in sorted(libs.items(), key=lambda kv: -kv[1]):
        print(f"  {lib:22s} {n:,}", file=sys.stderr)
    if ipr_ver not in ("unknown", "109.0"):
        print(f"\nNOTE: annotations are InterPro {ipr_ver}, while the shipped signature table is "
              f"109.0.\n  All nine NLR signatures added in 2.0.0 exist in 108.0, so nothing "
              f"NLR-relevant is lost.\n  Record this in run_info.yaml.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
