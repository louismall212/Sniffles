#!/usr/bin/env python3
#
# Sniffles2
# Edit SV candidates in existing .snf files (remove them or mark them as QC-failed)
#
"""
Edit SV candidates in an existing Sniffles2 .snf file, selected by their Sniffles IDs.

The IDs can come from:
  - a text file: one ID per line (first whitespace-separated token; '#' comments and blank lines are ignored),
    with or without the 'Sniffles2.' prefix, e.g. 'Sniffles2.DUP.1AS3' or 'DUP.1AS3'
  - a VCF (.vcf or .vcf.gz):
      * single-sample Sniffles VCF: the ID column is used
      * multi-sample (merged) Sniffles VCF: the per-sample FORMAT/ID of --sample is used
        (merged IDs such as 'Sniffles2.DUP.3M1' do not exist in any .snf file)

Actions:
  remove    drop the candidates from the .snf
  qc-false  keep them but set qc=False (and FILTER to --filter-label); multi-sample calling only
            reports an SV if at least one sample has a candidate with qc=True and FILTER=PASS

Examples:
  python -m sniffles.snf_edit --input s1.snf --output s1.edited.snf --ids ids.txt --action remove
  python -m sniffles.snf_edit --input s1.snf --output s1.edited.snf --ids s1.vcf.gz --action qc-false
  python -m sniffles.snf_edit --input s1.snf --output s1.edited.snf --ids merged.vcf --sample s1 --action qc-false

Must be run with the same Sniffles version that wrote the .snf (candidates are pickled SVCall objects).
Only open .snf files you trust: unpickling can execute code.
"""
import argparse
import gzip
import json
import os
import sys
import types
from collections import Counter

if __package__ in (None, ""):
    # Running as a plain script (python .../src/sniffles/snf_edit.py): prefer the installed sniffles,
    # because the .snf must be read with the version that wrote it. Fall back to this repository's copy.
    _here = os.path.dirname(os.path.abspath(__file__))
    sys.path = [p for p in sys.path if os.path.abspath(p or ".") != _here]  # don't import sniffles' modules as top-level
    try:
        import sniffles.snf  # noqa: F401
    except ImportError:
        sys.path.insert(0, os.path.dirname(_here))

from sniffles import sv
from sniffles.config import BUILD
from sniffles.snf import SNFile

ACTIONS = ("remove", "qc-false")


def strip_prefix(svid: str, prefix: str) -> str:
    return svid[len(prefix):] if prefix and svid.startswith(prefix) else svid


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def default_sample_name(snf_path: str, header: dict) -> str:
    """Same rule the sniffles main script uses for the sample column name in multi-sample calling."""
    sample_id = header["config"].get("sample_id")
    if sample_id is not None:
        return sample_id
    return os.path.splitext(os.path.basename(snf_path))[0]


def read_ids_txt(path: str, prefix: str) -> dict[str, tuple | None]:
    ids = {}
    with open_text(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                ids[strip_prefix(line.split()[0], prefix)] = None
    return ids


def read_ids_vcf(path: str, prefix: str, sample: str | None, merged_tolerance: int) -> dict[str, tuple | None]:
    """
    Returns {snf_id: (contig, pos, svtype, max_pos_diff) or None}. The location is used to double-check the match:
    exact for single-sample VCFs, within merged_tolerance bp for merged VCFs (merged POS comes from one of the
    merged samples' candidates, so it can differ from this sample's candidate).
    """
    ids = {}
    sample_col = None
    n_samples = 0
    with open_text(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = fields[9:]
                n_samples = len(samples)
                if n_samples > 1 or sample is not None:
                    if sample is None:
                        sys.exit(f"{path} has {n_samples} samples ({', '.join(samples)}); give --sample")
                    if sample not in samples:
                        sys.exit(f"Sample '{sample}' not in {path} (samples: {', '.join(samples)})")
                    sample_col = 9 + samples.index(sample)
                continue
            contig, pos, info = fields[0], int(fields[1]), fields[7]
            svtype = next((kv[7:] for kv in info.split(";") if kv.startswith("SVTYPE=")), None)
            if sample_col is None:
                raw_ids = [fields[2]]
                loc = (contig, pos, svtype, 0)
            else:
                fmt = fields[8].split(":")
                if "ID" not in fmt:
                    sys.exit(f"{path}: no FORMAT/ID field, cannot map merged calls to per-sample .snf IDs")
                values = fields[sample_col].split(":")
                i = fmt.index("ID")
                if i >= len(values) or values[i] in (".", "NULL"):
                    continue  # SV not present in this sample's .snf
                raw_ids = values[i].split(",")  # intra-sample merges list several IDs
                loc = (contig, pos, svtype, merged_tolerance)
            for raw in raw_ids:
                if raw not in (".", ""):
                    ids[strip_prefix(raw, prefix)] = loc
    return ids


def edit_snf(in_path: str, out_path: str, ids: dict, action: str, filter_label: str):
    # The reading code only needs this config attribute
    reader = SNFile(types.SimpleNamespace(combine_close_handles=False), open(in_path, "rb"), in_path)
    reader.read_header()
    header = reader.header

    found = Counter()
    mismatches = []
    changed = []
    new_index = {}
    chunks = []
    offset = 0
    removed = 0

    for contig, blocks in reader.index.items():
        for block_start in blocks:
            # read_blocks returns one dict per stored chunk; several tasks can write chunks for the same block
            for block in reader.read_blocks(contig, block_start):
                for svtype in sv.TYPES:
                    kept = []
                    for cand in block.get(svtype, []):
                        if cand.id not in ids:
                            kept.append(cand)
                            continue
                        found[cand.id] += 1
                        loc = ids[cand.id]
                        if loc is not None and not location_matches(loc, cand):
                            mismatches.append((cand.id, loc[:3], (cand.contig, cand.pos, cand.svtype)))
                            kept.append(cand)
                            continue
                        changed.append(cand)
                        if action == "remove":
                            removed += 1
                            continue
                        cand.qc = False
                        cand.filter = filter_label
                        kept.append(cand)
                    block[svtype] = kept

                data = gzip.compress(pickle_block(reader, block))
                new_index.setdefault(contig, {}).setdefault(block_start, []).append((offset, len(data)))
                chunks.append(data)
                offset += len(data)
    reader.close()

    header["index"] = new_index
    header["snf_candidate_count"] = max(0, header["snf_candidate_count"] - removed)
    with open(out_path, "wb") as out:
        # same header serialization as SNFileBase.write_results
        out.write((json.dumps(header, default=lambda obj: "<Unstored_Object>") + "\n").encode())
        for data in chunks:
            out.write(data)

    return found, mismatches, changed


def location_matches(loc: tuple, cand: sv.SVCall) -> bool:
    contig, pos, svtype, max_diff = loc
    # vcf.py writes POS as call.pos, or 1 if call.pos is 0
    return contig == cand.contig and svtype == cand.svtype and abs(pos - max(cand.pos, 1)) <= max_diff


def pickle_block(snf: SNFile, block: dict) -> bytes:
    """Serialize one block with SNFileBase.serialize_block."""
    snf.blocks = {0: block}
    return snf.serialize_block(0)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="input .snf")
    p.add_argument("--output", required=True, help="output .snf (must differ from input)")
    p.add_argument("--ids", required=True, help="IDs to edit: .txt (one per line) or .vcf/.vcf.gz")
    p.add_argument("--action", required=True, choices=ACTIONS)
    p.add_argument("--sample", default=None,
                   help="sample column to take FORMAT/ID from in a multi-sample VCF "
                        "(default for multi-sample VCFs: the name sniffles gives this .snf when merging)")
    p.add_argument("--filter-label", default="MANUAL", help="FILTER value for --action qc-false (default: MANUAL)")
    p.add_argument("--no-location-check", action="store_true",
                   help="do not require CHROM/SVTYPE/POS of the VCF record to match the .snf candidate")
    p.add_argument("--merged-pos-tolerance", type=int, default=None,
                   help="max. POS difference between a merged VCF record and this sample's candidate "
                        "(default: --combine-match-max stored in the .snf header, normally 1000)")
    p.add_argument("--report", default=None, help="write a TSV of the edited candidates")
    p.add_argument("--allow-overwrite", action="store_true", help="overwrite --output if it exists")
    args = p.parse_args(argv)

    if os.path.abspath(args.input) == os.path.abspath(args.output):
        sys.exit("--output must differ from --input")
    if os.path.exists(args.output) and not args.allow_overwrite:
        sys.exit(f"{args.output} exists; use --allow-overwrite")

    with open(args.input, "rb") as f:
        header = json.loads(f.readline())
    prefix = header["config"].get("id_prefix") or "Sniffles2."
    snf_build = header["config"].get("build")
    if snf_build != BUILD:
        print(f"WARNING: {args.input} was written by Sniffles {snf_build}, this is Sniffles {BUILD} "
              f"({os.path.dirname(sv.__file__)}). Use the same version if reading or merging fails.", file=sys.stderr)

    ids_lower = args.ids.lower()
    if ids_lower.endswith((".vcf", ".vcf.gz", ".vcf.bgz")):
        sample = args.sample
        if sample is None:
            with open_text(args.ids) as f:
                for line in f:
                    if line.startswith("#CHROM"):
                        if len(line.rstrip("\n").split("\t")) > 10:
                            sample = default_sample_name(args.input, header)
                            print(f"Multi-sample VCF: using sample column '{sample}'", file=sys.stderr)
                        break
        elif sample != default_sample_name(args.input, header):
            print(f"WARNING: --sample '{sample}' differs from this .snf's sample name "
                  f"'{default_sample_name(args.input, header)}'. Sniffles IDs repeat between samples, "
                  f"so make sure the column belongs to {args.input}.", file=sys.stderr)
        tolerance = args.merged_pos_tolerance
        if tolerance is None:
            tolerance = header["config"].get("combine_match_max") or 1000
        ids = read_ids_vcf(args.ids, prefix, sample, tolerance)
    else:
        ids = read_ids_txt(args.ids, prefix)
    if args.no_location_check:
        ids = dict.fromkeys(ids)

    found, mismatches, changed = edit_snf(args.input, args.output, ids, args.action, args.filter_label)

    missing = [i for i in ids if i not in found]
    dup_ids = [i for i, n in found.items() if n > 1]
    print(f"{len(ids)} IDs requested, {len(changed)} candidates {'removed' if args.action == 'remove' else 'set to qc=False'}"
          f" -> {args.output}", file=sys.stderr)
    if missing:
        print(f"WARNING: {len(missing)} IDs not found in {args.input}: {', '.join(missing[:10])}"
              f"{' ...' if len(missing) > 10 else ''}", file=sys.stderr)
    if dup_ids:
        print(f"WARNING: {len(dup_ids)} IDs occur more than once in the .snf; all were edited: {', '.join(dup_ids[:10])}",
              file=sys.stderr)
    for svid, want, have in mismatches:
        print(f"WARNING: {svid} NOT edited, location differs (VCF {want}, .snf {have}); "
              f"the VCF may come from another run. Use --no-location-check to force.", file=sys.stderr)

    if args.report:
        with open(args.report, "w") as r:
            r.write("id\tcontig\tpos\tend\tsvtype\tsvlen\tsupport\taction\n")
            for c in changed:
                r.write(f"{prefix}{c.id}\t{c.contig}\t{c.pos}\t{c.end}\t{c.svtype}\t{c.svlen}\t{c.support}\t{args.action}\n")


if __name__ == "__main__":
    main()
