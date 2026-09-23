"""Generate temporary synthetic FASTA and measure screening throughput."""

import argparse
import os
import random
import subprocess
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", type=int, default=20000)
    parser.add_argument("--length", type=int, default=5000)
    parser.add_argument("--unique", type=int, default=1024,
                        help="Number of distinct synthetic sequences to cycle through")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--retained", action="store_true", help="Also write retained FASTA")
    parser.add_argument("--adapter", action="store_true", help="Also scan a synthetic 46 bp adapter")
    default_binary = Path(__file__).resolve().parents[1] / "bin" / ("hificheck.exe" if os.name == "nt" else "hificheck")
    parser.add_argument("--binary", type=Path, default=default_binary)
    args = parser.parse_args()
    rng = random.Random(17)
    sequences = ["".join(rng.choices("ACGT", k=args.length)) for _ in range(args.unique)]
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        fasta = folder / "synthetic.fa"
        with fasta.open("w", buffering=1024 * 1024) as out:
            for i in range(args.reads):
                out.write(f">r{i}\n{sequences[i % args.unique]}\n")
        report = folder / "flags.tsv"
        command = [str(args.binary), "-i", str(fasta), "-o", str(report),
                   "-t", str(args.threads)]
        if args.retained:
            command.extend(["--retained-fasta", str(folder / "retained.fa")])
        if args.adapter:
            adapter_path = folder / "adapter.fa"
            adapter_path.write_text(">adapter\n" + "".join(rng.choices("ACGT", k=46)) + "\n")
            command.extend(["--adapters", str(adapter_path)])
        subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
