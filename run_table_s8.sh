#!/bin/bash
# run_table_s8.sh <dir with the replication files> [first] [last] -- runs `gfe_bm_rust table_s8` sample by sample,
# appending to output/table_s8_search.csv and skipping samples already there, so it can be stopped and restarted.
# Settings: Algorithm 1 with 1,000 starts; VNS with 20 Algorithm 1 starts, neighbourhoods up to 10 units, 10 rounds
# without improvement; seed 20261116.
cd "$(dirname "$0")"; DIR=$1; FIRST=${2:-1}; LAST=${3:-100}; OUT=output/table_s8_search.csv
cargo build --release -q
for s in $(seq $FIRST $LAST); do
  if [ -f $OUT ] && awk -F, -v s=$s 'NR>1 && $1==s' $OUT | grep -q .; then n=$(awk -F, -v s=$s 'NR>1 && $1==s' $OUT | wc -l); [ $n -eq 13 ] && continue; fi
  ./target/release/gfe_bm_rust table_s8 $DIR $s $s 1000 20 10 10 20261116 > output/sample_$s.tmp && \
    { [ -f $OUT ] || echo "sample,G,Q_authors,Q_local,Q_vns_authors,Q_alg1,Q_vns,Q_best,secs" > $OUT; grep -v '^sample' output/sample_$s.tmp >> $OUT; rm output/sample_$s.tmp; }
  echo "sample $s done $(date)"
done
