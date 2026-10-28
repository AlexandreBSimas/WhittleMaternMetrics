#!/bin/sh
# Reproduces Table 1 and Figures 1, 3 and 4 of the paper.
# Usage (from anywhere):
#   sh run_all.sh [cores]
# For each data set, the full-data fits and the cross-validation run at the
# same time, the cross-validation with `cores` worker processes (default 4,
# i.e. about 5 R processes and 10 GB of memory). A memory watchdog
# (watchdog.sh) stops the runs if they use more than 48 GB or free memory
# falls below 16 GB; finished splits are kept and a rerun resumes.
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
mkdir -p output/logs
C="${1:-4}"
bash watchdog.sh "exec/R.*(run_cv|fit_fulldata).R" 48 16 $$ > output/logs/watchdog.log 2>&1 &
for d in paris_k london; do
  # Wait only for this data set's R jobs, not for the watchdog.
  Rscript fit_fulldata.R "$d" > "output/logs/${d}_fulldata.log" 2>&1 & p1=$!
  Rscript run_cv.R "$d" lpo "$C" > "output/logs/${d}_lpo.log" 2>&1 & p2=$!
  wait $p1 $p2
done
Rscript summarise.R paris_k london > output/logs/summarise.log 2>&1
Rscript make_illustrations.R > output/logs/make_illustrations.log 2>&1
# The maps (Figure 4) need internet access for the basemap tiles and Google Chrome.
Rscript make_maps.R > output/logs/make_maps.log 2>&1 ||
  echo "make_maps.R failed (see output/logs); the other results are complete"
echo "Done: see output/"
