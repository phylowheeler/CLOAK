#!/usr/bin/env python3
"""
Realign one protein FASTA/MSA with MUSCLE5 while treating selected ambiguity
characters as missing for the re-alignment step.

Main behavior:
1. Read an input FASTA file that may already be aligned.
2. Remove user-specified missing/ambiguity characters from each sequence before
   sending sequences to MUSCLE.
3. Any taxon that becomes empty after that cleanup is NOT sent to MUSCLE.
4. Run MUSCLE on the remaining informative sequences.
5. Reinsert all-empty taxa as gap-only rows in the final alignment.
6. Preserve the original taxon order in the output FASTA.

Important detail:
- We STRIP X characters before MUSCLE rather than converting them to '-' first.
  That is intentional: MUSCLE discards input gaps anyway, so stripping and
  gap-conversion are equivalent for alignment input, while stripping is simpler.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class RecordMeta:
    rec_id: str
    description_suffix: str
    raw_seq: str
    cleaned_seq: str
    informative: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Realign one FASTA/MSA with MUSCLE5.")
    parser.add_argument("--input", required=True, help="Input FASTA/MSA file")
    parser.add_argument("--output", required=True, help="Output realigned FASTA file")
    parser.add_argument("--muscle", required=True, help="Path to MUSCLE5 executable")
    parser.add_argument("--threads", type=int, default=2, help="Threads for MUSCLE")
    parser.add_argument(
        "--seqtype",
        choices=["amino", "nt", "auto"],
        default="amino",
        help="Sequence type to pass to MUSCLE",
    )
    parser.add_argument(
        "--mode",
        choices=["align", "super5"],
        default="align",
        help="MUSCLE5 command to use",
    )
    parser.add_argument(
        "--missing-chars",
        default="Xx-.?*",
        help=(
            "Characters to treat as missing for re-alignment purposes. "
            "Default: Xx-.?*"
        ),
    )
    parser.add_argument(
        "--summary",
        default=None,
        help="Optional TSV file with a small per-locus summary",
    )
    return parser.parse_args()


def clean_sequence(seq: str, missing_chars: Sequence[str]) -> str:
    missing = {c.upper() for c in missing_chars}
    seq = seq.upper().replace(" ", "")
    return "".join(ch for ch in seq if ch not in missing)


def get_description_suffix(record: SeqRecord) -> str:
    desc = record.description or ""
    rec_id = record.id
    if not desc or desc == rec_id:
        return ""
    if desc.startswith(rec_id):
        return desc[len(rec_id):].lstrip()
    return desc


def load_records(infile: Path, missing_chars: Sequence[str]) -> List[RecordMeta]:
    raw_records = list(SeqIO.parse(str(infile), "fasta"))
    if not raw_records:
        raise ValueError(f"No FASTA records found in {infile}")

    ids = [r.id for r in raw_records]
    if len(ids) != len(set(ids)):
        dupes = sorted({x for x in ids if ids.count(x) > 1})
        raise ValueError(
            "Duplicate FASTA IDs detected. IDs must be unique for safe "
            f"reconstruction after MUSCLE. Duplicates: {', '.join(dupes[:10])}"
        )

    metas: List[RecordMeta] = []
    for record in raw_records:
        raw_seq = str(record.seq).upper()
        cleaned_seq = clean_sequence(raw_seq, missing_chars)
        metas.append(
            RecordMeta(
                rec_id=record.id,
                description_suffix=get_description_suffix(record),
                raw_seq=raw_seq,
                cleaned_seq=cleaned_seq,
                informative=(len(cleaned_seq) > 0),
            )
        )
    return metas


def write_temp_fasta(records: Sequence[RecordMeta], path: Path) -> None:
    seqrecords: List[SeqRecord] = []
    for meta in records:
        seqrecords.append(SeqRecord(Seq(meta.cleaned_seq), id=meta.rec_id, description=""))
    SeqIO.write(seqrecords, str(path), "fasta")


def run_muscle(
    informative_records: Sequence[RecordMeta],
    muscle_exe: Path,
    threads: int,
    seqtype: str,
    mode: str,
) -> Dict[str, str]:
    if len(informative_records) == 0:
        return {}

    if len(informative_records) == 1:
        only = informative_records[0]
        return {only.rec_id: only.cleaned_seq}

    with tempfile.TemporaryDirectory(prefix="muscle5_realign_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        temp_in = tmpdir_path / "input.fasta"
        temp_out = tmpdir_path / "output.fasta"

        write_temp_fasta(informative_records, temp_in)

        cmd = [str(muscle_exe)]
        if mode == "align":
            cmd.extend(["-align", str(temp_in)])
        else:
            cmd.extend(["-super5", str(temp_in)])
        cmd.extend(["-output", str(temp_out), "-threads", str(threads)])

        subprocess.run(cmd, check=True)

        aligned_records = list(SeqIO.parse(str(temp_out), "fasta"))
        if len(aligned_records) != len(informative_records):
            raise RuntimeError(
                "MUSCLE output record count does not match input record count: "
                f"input={len(informative_records)} output={len(aligned_records)}"
            )

        aligned_map = {rec.id: str(rec.seq).upper() for rec in aligned_records}
        expected_ids = {meta.rec_id for meta in informative_records}
        if set(aligned_map) != expected_ids:
            missing = sorted(expected_ids - set(aligned_map))
            extra = sorted(set(aligned_map) - expected_ids)
            raise RuntimeError(
                "MUSCLE output IDs do not match input IDs. "
                f"Missing: {missing[:10]} Extra: {extra[:10]}"
            )
        return aligned_map


def reconstruct_alignment(metas: Sequence[RecordMeta], aligned_map: Dict[str, str]) -> List[SeqRecord]:
    if aligned_map:
        aln_len = len(next(iter(aligned_map.values())))
    else:
        # All taxa were fully missing/ambiguous after cleanup.
        aln_len = 1

    final_records: List[SeqRecord] = []
    for meta in metas:
        if meta.rec_id in aligned_map:
            seq = aligned_map[meta.rec_id]
        else:
            seq = "-" * aln_len

        desc = meta.description_suffix
        final_records.append(SeqRecord(Seq(seq), id=meta.rec_id, description=desc))

    return final_records


def write_summary(
    summary_path: Path,
    infile: Path,
    metas: Sequence[RecordMeta],
    final_records: Sequence[SeqRecord],
    mode: str,
) -> None:
    informative_count = sum(meta.informative for meta in metas)
    missing_only_count = len(metas) - informative_count
    final_aln_len = len(final_records[0].seq) if final_records else 0

    if informative_count == 0:
        status = "all_missing_after_cleanup"
    elif informative_count == 1:
        status = "single_informative_sequence"
    else:
        status = f"muscle_{mode}"

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "input_file",
                "n_taxa_total",
                "n_informative_for_muscle",
                "n_all_missing_after_cleanup",
                "final_alignment_length",
                "status",
            ]
        )
        writer.writerow(
            [
                str(infile),
                len(metas),
                informative_count,
                missing_only_count,
                final_aln_len,
                status,
            ]
        )


def main() -> int:
    args = parse_args()

    infile = Path(args.input)
    outfile = Path(args.output)
    muscle_exe = Path(args.muscle)

    if not infile.is_file():
        raise FileNotFoundError(f"Input file not found: {infile}")
    if not muscle_exe.exists():
        raise FileNotFoundError(f"MUSCLE executable not found: {muscle_exe}")

    metas = load_records(infile, args.missing_chars)
    informative = [meta for meta in metas if meta.informative]

    aligned_map = run_muscle(
        informative_records=informative,
        muscle_exe=muscle_exe,
        threads=args.threads,
        seqtype=args.seqtype,
        mode=args.mode,
    )

    final_records = reconstruct_alignment(metas, aligned_map)

    outfile.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(final_records, str(outfile), "fasta")

    if args.summary:
        write_summary(Path(args.summary), infile, metas, final_records, args.mode)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
