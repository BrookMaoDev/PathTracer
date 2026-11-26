#!/bin/sh
# Runs every render script one after another (each render already uses all cores). Each script's
# output goes to logs/NAME.log; the images end up in renderings/.
cd "$(dirname "$0")"
mkdir -p ../logs

status=0
for s in cornell_noIS_noES cornell_IS_noES cornell_IS_ES; do
  if ./$s.sh > ../logs/$s.log 2>&1; then
    tail -n 1 "../logs/$s.log"
  else
    echo "$s: FAILED, see logs/$s.log"
    status=1
  fi
done
exit $status
