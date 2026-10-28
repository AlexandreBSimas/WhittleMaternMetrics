# Whittle–Matérn metrics and isotropic Gaussian processes on compact metric graphs

Code and data to reproduce the numerical results of

> David Bolin and Alexandre B. Simas (2026). *Whittle–Matérn metrics and
> isotropic Gaussian processes on compact metric graphs.* arXiv preprint
> [arXiv:2610.35732](https://arxiv.org/abs/2610.35732).

The Whittle–Matérn metric `d_kappa`, `kappa >= 0`, is the variogram of a
first-order Whittle–Matérn field on a metric graph; `d_0` is its limit as
`kappa -> 0`. For every completely monotone function `r`, `r(d_kappa)` is a
valid isotropic covariance function on every compact metric graph. The code
fits such models to traffic data from Paris and London, and compares them
with other Gaussian process models by likelihood and cross-validation.

| Paper | Produced by | Output |
|---|---|---|
| Figure 1 | `make_illustrations.R` | `output/metric_multigraph.png` |
| Figure 3 | `make_illustrations.R` | `output/samples_isotropic.png` |
| Figure 4 | `make_maps.R` | `output/paris_occupancy.png`, `output/london_speed.png` |
| Table 1, `p` and `-log L` | `fit_fulldata.R`, `summarise.R` | `output/<dataset>_fulldata_fit.csv`, `output/tex/fulldata_rows.tex` |
| Table 1, RMSE and CRPS | `run_cv.R`, `summarise.R` | `output/<dataset>_lpo_summary.csv`, `output/tex/<dataset>_lpo_rows.tex` |

Figure 2 is a diagram drawn in the LaTeX source. All results are already in
`output/`.

## Requirements

- **R** (tested with R 4.6.1) with `MetricGraph` (1.6.0.9000), `rSPDE`
  (2.6.0.9001), `Matrix` and `sf`.
- **Figures 1 and 3:** also `MASS`, `ggplot2`, `viridis` and `patchwork`.
- **Figure 4:** also `leaflet`, `webshot2`, `htmlwidgets` and `viridis`, with Google
  Chrome and internet access for the basemap tiles.
- **Data preparation** (optional, `data_prep/`): curl, Python 3 with `pyproj`,
  R with `jsonlite`, and internet access.

## Running

From the repository root:

```sh
sh run_all.sh 4      # everything, with 4 cross-validation workers per data set
```

This takes about 1.5 hours on a laptop. Step by step:

```sh
Rscript fit_fulldata.R paris_k        # full-data fits (also FEM on a 0.1 km mesh)
Rscript fit_fulldata.R london
Rscript run_cv.R paris_k lpo 4        # leave 75% and 90% of the sites out, 25 draws each
Rscript run_cv.R london lpo 4
Rscript summarise.R paris_k london    # summary CSVs and the LaTeX rows of Table 1
Rscript make_illustrations.R          # Figures 1 and 3
Rscript make_maps.R                   # Figure 4
```

Each cross-validation split is saved under `output/<dataset>_lpo/splits/` as
it finishes. A rerun reuses the saved splits, so an interrupted run resumes
where it stopped, and with the splits in this repository `run_cv.R` only
recomputes the score tables. Delete `output/<dataset>_lpo/` to recompute
the cross-validation from scratch.

**Memory.** Each R process needs about 2 GB. `run_all.sh` runs the full-data
fits and the cross-validation of a data set at the same time, and starts
`watchdog.sh`, which stops the R jobs if together they use more than 48 GB
or free memory drops below 16 GB (saved splits are kept). Set
`VECLIB_MAXIMUM_THREADS=1` (macOS) or `OPENBLAS_NUM_THREADS=1` when running
several R processes in parallel, as `run_all.sh` does. To run a long job
detached from the terminal, use `python3 detach.py <log> Rscript run_cv.R ...`.

## Data (`data/`)

| Data set | Response | Locations × days | Network |
|---|---|---|---|
| Paris (`paris_k`) | log occupancy rate (%), 08:00–09:00, working days, 4 Sep 2025 – 16 Oct 2026 | 314 × 144 | counted road segments of the City of Paris, 478 km |
| London (`london`) | mean speed (mph), 17:00–18:00, every day, 1 Jan – 30 Jun 2025 | 256 × 88 | OpenStreetMap motorways and trunk roads with their links, 2488 km |

**Paris** (`data/paris/`). Permanent traffic counters of the City of Paris
(opendata.paris.fr, datasets `comptages-routiers-permanents` and
`referentiel-comptages-routiers`), downloaded on 23 October 2026. The
network is built from the City's reference geometry of the counted road
segments, not from OpenStreetMap.

| File | Content |
|---|---|
| `network_edges.csv` | One edge per pair of junctions joined by a counted segment (either direction); polyline in Lambert-93 km, with the end points at the junctions |
| `sites.csv` | One site per edge, from the direction with the better coverage, at the edge midpoint |
| `panel_k.csv` | Occupancy (%), days × sites, complete |
| `panel_q.csv` | Flow (vehicles per hour), same days and sites; not used in the paper |

The panels keep the sites with values on at least 90% of the working days,
and then remove days and sites, one at a time, until no values are missing.

**London** (`data/london/`). National Highways WebTRIS speeds and the
OpenStreetMap road network (© OpenStreetMap contributors, ODbL).

| File | Content |
|---|---|
| `panel_1700.csv` | Mean speed, 17:00–18:00, days × sites, complete. A site-day needs at least 3 of the four 15-minute values; sites need values on at least 95% of the days |
| `sites_selected.csv` | Active MIDAS sites in the box 0.65°W–0.40°E, 51.22°N–51.78°N, at least 0.3 km apart |
| `osm_london.rds` | OpenStreetMap motorway, motorway_link, trunk and trunk_link lines (sf object) |

**Rebuilding the data.** `sh data_prep/paris/prepare.sh` and
`sh data_prep/london/prepare.sh` download the raw data, rebuild the panels
and copy them to `data/`. The Paris source keeps only about the last 13
months of hourly data, so a rerun of the Paris preparation covers a later
period and gives different panels from those in `data/paris/`.

## Files

| File | Content |
|---|---|
| `R/data.R` | Loaders for both data sets: the metric graph with the observation locations as vertices, the responses, mean covariates (easting, northing), coordinates and geodesic distances |
| `R/models.R` | The metrics `d_0`, `d_kappa` and `d_res`; the nine covariance models; profile likelihood; starting values; exact Gaussian prediction |
| `R/design.R` | Cross-validation splits; RMSE, CRPS, log scores and coverage |
| `fit_fulldata.R` | Full-data maximum likelihood fits |
| `run_cv.R` | Cross-validation for one data set |
| `summarise.R` | Summary CSVs and the LaTeX rows of Table 1 |
| `make_illustrations.R` | Figures 1 and 3 |
| `make_maps.R` | Figure 4 |
| `run_all.sh`, `watchdog.sh`, `detach.py` | Running everything, the memory guard, detached runs |
| `data_prep/` | Download and panel scripts for both data sets |

## Method in brief

- **Models.** The mean is intercept + easting + northing, shared over the
  days, which are treated as independent replicates. All models have
  independent Gaussian noise. The covariances are:
  - exponential, generalised Cauchy and Matérn (`nu <= 1/2`), each in `d_0`
    and in `d_kappa`;
  - the Whittle–Matérn field with `alpha = 1`;
  - the Whittle–Matérn field with `alpha` estimated in (1/2, 5/2], using the
    rSPDE rational approximation of order 2 on a 0.25 km finite element mesh;
  - exponential in the resistance distance of the graph as given (`d_res`).
- **Estimation.** Maximum likelihood, with the mean coefficients profiled out
  by generalised least squares. `kappa` is estimated without bounds. The
  Cauchy and Matérn fits are also started from the fitted exponential model,
  so their likelihoods are at least as high.
- **Computation of `d_kappa`.** The precision matrix `Q^kappa` of the field at
  the vertices is assembled in closed form from the edge lengths (it agrees
  with `MetricGraph::spde_precision` to about 1e-11), and its sparse Cholesky
  factorisation is reused across values of `kappa`.
- **London graph.** The degree-2 vertices of the OpenStreetMap lines are
  merged after the sites are inserted as vertices (15,325 to 3,990 vertices),
  which reduces memory use from about 24 GB to about 2 GB per process.
- **Cross-validation.** 75% or 90% of the locations are held out, 25 random
  draws each. All parameters are re-estimated on each training set, the
  held-out values are predicted from the exact conditional distribution, and
  the scores are averaged over the draws.

`run_cv.R` also implements a buffered leave-region-out design
(`Rscript run_cv.R <dataset> lro <cores>`), which is not used in the paper.

## Citation

```bibtex
@article{bolin_simas2026wmmetrics,
  title   = {Whittle--Mat\'ern metrics and isotropic Gaussian processes on compact metric graphs},
  author  = {Bolin, David and Simas, Alexandre B.},
  journal = {arXiv preprint arXiv:2610.35732},
  year    = {2026},
  doi     = {10.48550/arXiv.2610.35732}
}
```
