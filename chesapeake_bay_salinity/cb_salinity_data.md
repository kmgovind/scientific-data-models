# Chesapeake Bay Salinity Data

## Overview

This directory contains Chesapeake Bay salinity data derived from [NOAA's Chesapeake Bay Operational Forecast System (CBOFS)](https://tidesandcurrents.noaa.gov/ofs/cbofs/cbofs.html), along with a Julia download script and a Jupyter notebook for exploration and visualization.

Salinity is typically reported in Practical Salinity Units (PSU). In the Chesapeake Bay, it varies in space and time because of freshwater inflow, ocean exchange, tides, and mixing.

This makes the dataset useful for environmental field estimation, forecasting, and adaptive sampling experiments.

---

## Data Source

The source is NOAA's Chesapeake Bay Operational Forecast System (CBOFS), a hydrodynamic model that produces forecast and nowcast products for the bay and surrounding waters. The data in this directory focus on salinity outputs stored as NetCDF files.

---

## Data Layout

The salinity fields are organized as spatiotemporal model outputs that can be thought of conceptually as:

```text
S(latitude, longitude, depth, time)
```

where the grid varies across location, depth, and forecast time.

The downloaded files in this directory are stored under `datafiles/` and use names like:

```text
chesapeake_salinity_YYYYMMDD_tHHzz_nNNN.nc
```

For example, a file such as `chesapeake_salinity_20260616_t00z_n001.nc` represents one model output for a specific date, cycle, and forecast hour.

---

## Files In This Directory

### `cb_salinity_download.jl`

Downloads CBOFS salinity NetCDF files into `datafiles/`.

### `cb_salinity_data_viz.ipynb`

Notebook for exploratory analysis and visualization of the salinity data.

### `chesapeake_salinity_animation.mp4`

Rendered animation showing salinity evolution through time.

### `datafiles/`

Folder containing the downloaded NetCDF files used by the notebook and animation.

### `cb_salinity_data.md`

This document.

---

## Suggested Analyses

The notebook in this directory can be used to explore:

* Surface salinity maps at a fixed time and depth
* Time series at a fixed location
* Vertical salinity structure across depth
* Animations of salinity changes through time

These views are useful for understanding the salt-freshwater transition, spatial gradients, and short-term variability in the bay.

---

## Notes

The download script targets NOAA CBOFS catalog and file-server URLs and writes only the files that are available for the requested date range. Because the data are large, the repository may contain only a subset of outputs in `datafiles/`.

---

## Reference

NOAA Chesapeake Bay Operational Forecast System (CBOFS)

For model documentation and data access details, see NOAA's official CBOFS documentation.
