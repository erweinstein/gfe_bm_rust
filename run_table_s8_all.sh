#!/bin/bash
# run_table_s8_all.sh <dir with the replication files> [first] [last] [out] -- the Table S.VIII search in one call:
# every (sample, G) job goes to one queue (largest G first) shared by all threads, rows are appended to <out> as
# they finish, and jobs already in <out> are skipped, so the call can be stopped and restarted. Uses the
# incremental local search; the results are the same as run_table_s8.sh's (checked on all 1,300 jobs).
# Settings as run_table_s8.sh. RAYON_NUM_THREADS sets the threads (default: all).
# Added 24 November 2026.
cd "$(dirname "$0")"; DIR=$1; FIRST=${2:-1}; LAST=${3:-100}; OUT=${4:-output/table_s8_search_incremental.csv}
cargo build --release -q
./target/release/gfe_bm_rust table_s8 "$DIR" $FIRST $LAST 1000 20 10 10 20261116 --incremental --out "$OUT"
