"""Synthetic end-to-end checks for the compiled FASTA screen."""

import csv
import json
import os
import random
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BIN = ROOT / "bin" / ("hificheck.exe" if os.name == "nt" else "hificheck")
BIN = Path(os.environ.get("HIFICHECK_BINARY", DEFAULT_BIN))


def dna(length, seed):
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@unittest.skipUnless(BIN.is_file(), "Build HiFiCheck first: make")
class TestFastaScreen(unittest.TestCase):
    def run_screen(self, records, adapter=None, extra=(), retained=False):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fasta = tmp / "reads.fa"
            report = tmp / "report.tsv"
            fasta.write_text("".join(f">{name} extra metadata\n{seq[:777]}\n{seq[777:]}\n" for name, seq in records))
            command = [str(BIN), "-i", str(fasta), "-o", str(report), "-t", "2", *extra]
            if retained:
                command.extend(["--retained-fasta", str(tmp / "retained.fa")])
            if adapter:
                adapters = tmp / "adapter.fa"
                adapters.write_text(f">library_adapter\n{adapter}\n")
                command.extend(["--adapters", str(adapters)])
            run = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            with report.open(newline="") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            summary = json.loads((tmp / "report.tsv.summary.json").read_text())
            summary["_quarantine_ids"] = (tmp / "report.tsv.quarantine.ids").read_text().splitlines()
            summary["_direct_repeat_ids"] = (tmp / "report.tsv.direct_repeat.ids").read_text().splitlines()
            summary["_retained_fasta"] = (tmp / "retained.fa").read_text() if retained else None
            return rows, summary

    def test_clean_and_multiline_fasta(self):
        rows, summary = self.run_screen([("clean", dna(5300, 1))])
        self.assertEqual(rows, [])
        self.assertEqual(summary["reads"], 1)
        self.assertEqual(summary["bases"], 5300)
        self.assertEqual(summary["keep"], 1)
        self.assertEqual(summary["fasta_quality_and_bam_tags"], "unavailable")

    def test_internal_adapter(self):
        adapter = dna(46, 5)
        rows, summary = self.run_screen([("adapter_read", dna(3000, 1) + adapter + dna(3000, 2))], adapter)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["qc_action"], "exclude_adapter")
        self.assertEqual(rows[0]["adapter_start_0"], "3000")
        self.assertEqual(rows[0]["adapter_end_0"], "3046")
        self.assertEqual(summary["adapter_hits"], 1)
        self.assertEqual(summary["_quarantine_ids"], ["adapter_read"])

    def test_direct_repeat_is_retained(self):
        unit = dna(2900, 11)
        rows, summary = self.run_screen([("array", dna(400, 3) + unit + dna(200, 7) + unit + dna(400, 4))])
        self.assertEqual(rows[0]["qc_action"], "keep_direct_repeat")
        self.assertEqual(rows[0]["repeat_orientation"], "+")
        self.assertEqual(summary["direct_repeat_hits"], 1)
        self.assertEqual(summary["_quarantine_ids"], [])
        self.assertEqual(summary["_direct_repeat_ids"], ["array"])

    def test_similar_tandem_copies_are_retained(self):
        unit = dna(3600, 23)
        changed = list(unit)
        for pos in range(0, len(changed), 100):
            changed[pos] = "A" if changed[pos] != "A" else "C"
        rows, summary = self.run_screen([("two_operons", unit + dna(200, 31) + "".join(changed))])
        self.assertEqual(rows[0]["qc_action"], "keep_direct_repeat")
        self.assertEqual(summary["_quarantine_ids"], [])

    def test_inverted_repeat(self):
        unit = dna(2900, 12)
        rows, summary = self.run_screen([("readthrough", dna(400, 3) + unit + dna(200, 7) + rc(unit) + dna(400, 4))])
        self.assertEqual(rows[0]["qc_action"], "quarantine")
        self.assertEqual(rows[0]["repeat_orientation"], "-")
        self.assertEqual(summary["inverted_repeat_hits"], 1)
        self.assertEqual(summary["_quarantine_ids"], ["readthrough"])

    def test_near_terminal_palindrome_is_excluded(self):
        unit = dna(2900, 12)
        rows, summary = self.run_screen([("palindrome", unit + dna(200, 7) + rc(unit))])
        self.assertEqual(rows[0]["qc_action"], "exclude_palindrome", rows[0])
        self.assertEqual(summary["exclude_palindrome"], 1)
        self.assertEqual(summary["_quarantine_ids"], ["palindrome"])

    def test_all_mode_and_ambiguous_bases(self):
        rows, summary = self.run_screen([("with_ambiguity", dna(2200, 1) + "NRYWSKMBVDH" + dna(2200, 2))], extra=("--all",))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["qc_action"], "keep")
        self.assertEqual(summary["reads"], 1)

    def test_retained_fasta_keeps_direct_arrays_and_preserves_headers(self):
        unit = dna(2900, 12)
        clean = dna(5000, 30)
        rows, summary = self.run_screen([
            ("clean", clean),
            ("array", unit + dna(200, 9) + unit),
            ("palindrome", unit + dna(200, 7) + rc(unit)),
        ], retained=True)
        self.assertIn(">clean extra metadata\n" + clean + "\n", summary["_retained_fasta"])
        self.assertIn(">array extra metadata\n", summary["_retained_fasta"])
        self.assertNotIn(">palindrome", summary["_retained_fasta"])
        self.assertEqual(summary["_quarantine_ids"], ["palindrome"])

    def test_existing_output_requires_force(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            fasta = folder / "input.fa"
            report = folder / "report.tsv"
            fasta.write_text(">r1\n" + dna(4500, 1) + "\n")
            report.write_text("previous report\n")
            command = [str(BIN), "-i", str(fasta), "-o", str(report), "-t", "1"]
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
            self.assertEqual(report.read_text(), "previous report\n")
            self.assertEqual(subprocess.run(command + ["--force"], capture_output=True).returncode, 0)
            self.assertTrue((folder / "report.tsv.summary.json").is_file())


if __name__ == "__main__":
    unittest.main()
