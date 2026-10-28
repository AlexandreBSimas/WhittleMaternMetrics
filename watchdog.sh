#!/bin/bash
# Memory guard for long runs: every 20 s, if the R jobs matching PATTERN use
# more than MAX_GB in total, or the free memory drops below MIN_FREE_GB,
# stop them (the per-split results already saved are kept).
# Usage: bash watchdog.sh "<pattern>" [MAX_GB] [MIN_FREE_GB] [OWNER_PID]
# With OWNER_PID (e.g. the PID of run_all.sh), the guard keeps running while
# that process is alive, also between phases when no R job is running.
# Without it, the guard stops as soon as no matching job is running.
PAT="$1"; MAX_GB="${2:-48}"; MIN_FREE_GB="${3:-16}"; OWNER="${4:-}"
alive() {
  if [ -n "$OWNER" ]; then kill -0 "$OWNER" 2> /dev/null
  else ps -axo command= | grep -E "$PAT" | grep -qv grep; fi
}
while alive; do
  used=$(ps -axo rss=,command= | grep -E "$PAT" | grep -v grep | awk '{s+=$1} END {print int(s/1048576)}')
  free=$(vm_stat | awk '/free|inactive|speculative/ {gsub("\\.","",$NF); s+=$NF} END {print int(s*16384/1073741824)}')
  echo "$(date '+%H:%M:%S') R jobs ${used} GB, free ${free} GB"
  if [ "$used" -gt "$MAX_GB" ] || [ "$free" -lt "$MIN_FREE_GB" ]; then
    echo "memory limit reached: stopping jobs"
    pkill -TERM -f "$PAT"; sleep 5; pkill -KILL -f "$PAT"; exit 1
  fi
  sleep 20
done
