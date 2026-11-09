# Status

Updated: 2026-11-09

Local algorithm and robustness review completed for this revision.
Fixed inverted-repeat priority, IUPAC preservation, numeric narrowing and span
multiplication overflow, input/output path aliases, shared stdin rejection,
worker exception propagation and thread-start cleanup, summary close checking,
and required build flags under Conda. Exact internal adapter hits return early.
Validation: make test (12 tests); GCC UBSan build (12 tests, no diagnostics).
Synthetic benchmark: 20,000 x 5,000 bp, adapter scan and retained output:
1 worker 9,177.0 reads/s; 8 workers 38,270.6 reads/s (single runs, not a speedup
comparison against the previous code).
Remaining limitations: 100 bp repeat bins aggregate sparse seed evidence without
alignment; high-copy seed cap can miss repeats; only one prioritized repeat is
reported; batch limits do not cap a single read's memory. Real-data sensitivity,
specificity and fault-injection tests for allocation/thread failures not run.
No remote runs submitted.
