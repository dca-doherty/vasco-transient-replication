# vasco-transient-replication
Publicly Sourced Data Replication of VASCO Palomar Transients Script 

# Reproducing the VASCO Transient Catalog from POSS-I Plates

## What This Is

This document describes an open-source pipeline that independently reproduces the transient catalog published by Solano, Villarroel, and Rodrigo (2022) in Monthly Notices of the Royal Astronomical Society. Their catalog identified 5,399 "vanishing" sources on photographic plates from the First Palomar Observatory Sky Survey (POSS-I), conducted between 1949 and 1958. These are point sources that appeared on the original red-sensitive plates but have no counterpart in any modern astronomical survey.

The pipeline downloads the original digitized plate scans from NASA's Infrared Science Archive (IRSA), extracts sources using SExtractor and PSFEx, applies morphological quality filters, performs a red-versus-blue plate comparison to isolate true transients, and then crossmatches surviving candidates against modern catalogs (Gaia, PanSTARRS, SuperCOSMOS) to remove persistent objects. Infrared catalog crossmatching is available as an option but disabled by default, since the published 107,000-source VASCO catalog was never crossmatched against IR (confirmed by Villarroel). The pipeline also applies a southern hemisphere declination cutoff and supports cross-plate deduplication for multi-plate runs. The output is a list of transient candidates that can be compared directly against the published VASCO catalog.

The goal is full reproducibility. Anyone with a Linux machine, the required software, and enough disk space should be able to run this pipeline and arrive at the same set of transients that VASCO published.

## Why This Matters

The VASCO catalog is the foundation for several lines of active research, including a 2025 paper in Scientific Reports that found statistical correlations between transient detections and nuclear weapons testing. Independent reproduction of the catalog strengthens all downstream analyses. It also makes the catalog auditable. Until now, the extraction pipeline was described in the literature but the code was never released. This pipeline fills that gap.

## What You Need

### Hardware

You need a machine with at least 16 GB of RAM. PSFEx builds a point spread function model from tens of thousands of source stamps and will use around 800 MB of memory per plate. SExtractor itself is lightweight, but the intermediate catalogs for a full plate can reach 400 MB when VIGNET stamps are included.

Disk space is the real constraint. Each full POSS-I plate scan is roughly 400 MB. There are 936 red plates and 936 blue plates in the survey. If you download everything, that is approximately 750 GB of raw plate data. The intermediate SExtractor catalogs, PSFEx models, and cached catalog query results will roughly double that. Plan for 1.5 TB of free disk space if you intend to process the entire survey.

A solid-state drive will make a noticeable difference. SExtractor reads the entire plate into memory, but the catalog queries write and read many small cache files. On a spinning disk, the I/O overhead adds up across hundreds of plates.

### Operating System

The pipeline runs on Linux. It was developed and tested on Ubuntu 24.04 under Windows Subsystem for Linux (WSL2). It should work on any recent Debian-based distribution. macOS may work with Homebrew-installed dependencies, but this has not been tested. Windows native is not supported because SExtractor and PSFEx do not have Windows builds.

### Software Dependencies

Install the following before running the pipeline:

**SExtractor** (also called source-extractor in newer packages):
```
sudo apt install source-extractor
```
The pipeline checks for the command names `source-extractor`, `sextractor`, and `sex` in that order. Any of these will work.

**PSFEx**
```
sudo apt install psfex
```
PSFEx builds a spatially varying point spread function model from the extracted sources. This is used in the second SExtractor pass to compute SPREAD_MODEL, which separates real astronomical sources from cosmic ray hits and plate defects.

**Python 3.10 or later** with the following packages
```
pip install numpy pandas astropy astroquery matplotlib scipy
```

Specific versions are not critical, but the pipeline was tested with numpy 1.24, pandas 2.0, astropy 6.0, and astroquery 0.4.7. Earlier versions may work but have not been verified.

**A Python virtual environment is recommended** to avoid conflicts with system packages
```
python3 -m venv venv
source venv/bin/activate
pip install numpy pandas astropy astroquery matplotlib scipy
```

### Network Access

The pipeline downloads plate scans from IRSA and queries VizieR for catalog crossmatching. You need a reliable internet connection. Catalog queries for a single plate pull several million rows across three primary catalogs (or seven with `--ir` enabled) and can take 5 to 15 minutes depending on your connection and VizieR server load. The pipeline caches all query results to disk, so subsequent runs of the same plate skip the network entirely.

If you are behind a corporate firewall or proxy, you may need to configure the `HTTP_PROXY` and `HTTPS_PROXY` environment variables for both the plate downloads and the astroquery VizieR calls.


## How the Pipeline Works

### Step 1: Download the Plates

Each POSS-I field was photographed twice: once through a red-sensitive emulsion (plate prefix XE) and once through a blue-sensitive emulsion (plate prefix XO). The field numbers match. For example, field 582 has red plate XE582 and blue plate XO582.

The digitized scans are hosted by IRSA at
- Red: `https://irsa.ipac.caltech.edu/data/DSS/images/dss1red/dss1red_XE{NUM}.fits`
- Blue: `https://irsa.ipac.caltech.edu/data/DSS/images/dss1blue/dss1blue_XO{NUM}.fits`

Each plate is a single FITS file, roughly 14,000 by 14,000 pixels, about 400 MB uncompressed. The pipeline checks if the file already exists before downloading, so interrupted runs can be resumed without re-downloading.

### Step 2: Source Extraction (SExtractor + PSFEx)

Source extraction happens in two passes per plate.

**Pass 1** runs SExtractor without a PSF model. This extracts all detectable sources and includes 35x35 pixel image stamps (VIGNET) around each source. These stamps are needed by PSFEx to build the PSF model. The pass 1 catalog is large (300 to 400 MB) because of the stamps, but it is only used as input to PSFEx and is not carried forward.

**PSFEx** reads the pass 1 catalog and builds a spatially varying PSF model across the plate. Photographic plates have PSFs that change across the field due to optical distortions and emulsion variations. PSFEx fits a polynomial model (degree 2 in X and Y) to capture this. This step takes 1 to 12 minutes per plate and uses about 800 MB of RAM. The output is a `.psf` file used in pass 2. Diagnostic check images and XML output are disabled to prevent I/O hangs on large catalogs.

**Pass 2** runs SExtractor again, this time with the PSF model loaded. This enables the SPREAD_MODEL parameter, which measures how much each source's light profile deviates from the local PSF. Real stars have SPREAD_MODEL near zero. Cosmic rays, hot pixels, and plate defects tend to have strongly negative values. Pass 2 does not include VIGNET stamps, so the catalog is much smaller (under 10 MB) and runs in about 30 to 40 seconds.

**A note on WCS coordinates:** The POSS-I plate scans use a polynomial World Coordinate System (WCS) stored as PLT* keywords in the FITS header. SExtractor cannot read this format and will output all-zero coordinates. The pipeline works around this by using Astropy's WCS module to convert pixel coordinates to sky coordinates after extraction. Astropy reads the PLT* polynomial correctly.

### Step 3: Morphological Filters

The pipeline applies the same quality filters described in Solano et al. (2022)

- **FLAGS = 0**: No SExtractor warnings (blending, saturation, edge proximity)
- **SPREAD_MODEL > -0.01**: Rejects cosmic rays and sharp defects. The published threshold is -0.002, but our PSFEx models on photographic plates produce slightly different SPREAD_MODEL distributions than the original pipeline. At -0.002, lose roughly half of the known VASCO sources. At -0.01, retain 97% of them while still removing the worst artifacts.
- **FWHM between 2 and 7 pixels**: Rejects unresolved noise spikes (too small) and extended objects (too large)
- **ELONGATION < 1.3**: Rejects streaks and satellite trails
- **Symmetry < 2 pixels**: Rejects sources where the bounding box is not roughly square
- **SNR >= 30**: Rejects faint or noisy detections

One filter that is deliberately omitted is the Median Absolute Deviation (MAD) clipping on FWHM and elongation. Solano describes computing the median and MAD within 30-arcmin tessellated regions and removing sources deviating more than 2 sigma from the local median. Testing showed that this removes real transients whose morphology legitimately differs from the stellar population. Transients are inherently morphological outliers, and MAD clipping penalises exactly that. Applying MAD clipping globally reduced recall by 30 percent and locally by 26 percent with no meaningful improvement in F1 score. The hard cuts on FWHM, elongation, and SPREAD_MODEL handle the worst artifacts without this step.

### Step 4: Red Versus Blue Comparison

This is the critical step that separates real transients from persistent sources. A source that appears on the red plate but not on the blue plate is a transient candidate, because it was only present during one of the two exposures (taken at different times). A source that appears on both plates is a persistent astronomical object and gets removed.

The pipeline crossmatches the filtered red and blue source lists at a 5 arcsecond radius. Any red source with a blue counterpart is marked as persistent and excluded. What remains are the transient candidates.

This step typically removes thousands of sources and is the primary reason the final count is much lower than the post-filter count. Without it, you will have roughly 50 to 100 percent more candidates than the published catalog.

Not all POSS-I fields have blue plate scans available on IRSA. If the blue plate cannot be downloaded, the pipeline falls back to red-only mode and skips this step. Use the `--red-only` flag to force this behavior intentionally.

### Step 5: Catalog Crossmatch

Surviving candidates are crossmatched against modern astronomical catalogs to remove any source that has a known counterpart today. The catalogs are queried via VizieR with a 5 degree search radius centered on the plate:

**Primary catalogs** (always run, queried in parallel)
- Gaia EDR3 (300,000+ sources per plate). Proper motions are corrected from the Gaia epoch (2016) back to the POSS-I epoch (approximately 1951) to account for stellar motion over 65 years.
- PanSTARRS DR2 (1.5 to 2 million sources per plate)

**Infrared catalogs** (optional, disabled by default, enable with `--ir`)
- AllWISE
- 2MASS
- CatWISE2020
- UKIDSS-LAS (northern sky only)

Villarroel confirmed that the published 107,000-source VASCO catalog was never crossmatched against infrared catalogs. The IR step is kept as an option for users who want stricter filtering, but it's off by default to match the actual methodology that produced the published catalog.

**Second-epoch optical**:
- SuperCOSMOS (digitization of POSS-II plates from the 1980s and 1990s)

All query results are cached to disk in the `cache/` subdirectory. If you rerun the pipeline on the same plate, the cached results are loaded instantly. Clearing the cache forces fresh queries.

### Step 6: Southern Hemisphere Exclusion

POSS-I included some southern fringe plates that extend below declination -30 degrees. These were excluded from the published VASCO catalog. The pipeline applies a declination cutoff (default -30 degrees) and removes any candidate below that limit. This can be adjusted with the `--dec-min` flag if needed, but the default matches the published boundary.

### Step 7: Cross-Plate Deduplication

POSS-I plates overlap by about 1 degree at their edges. The same transient source can appear on two to four adjacent plates. For multi-plate runs, the pipeline can deduplicate by checking new candidates against previously processed results. Any candidate within 5 arcseconds of an existing detection in a previous plate's output is flagged as a duplicate and removed. This is only active when you provide a `--dedup-dir` pointing to the directory containing your previous `transients_*.csv` files.

### Step 8: VASCO Comparison

The final candidate list is crossmatched against the published VASCO catalog at 5 arcseconds. The pipeline reports recall (what fraction of VASCO entries recovered), precision (what fraction of our candidates match VASCO), and F1 score.

## Running the Pipeline

### Single Plate

```
python3 poss_transients.py --plate 582
```
This processes plate XE582 (red) and XO582 (blue), applies all filters, runs the crossmatch, and saves the results to `./output/transients_582.csv`.

### Validating Against VASCO

If you have a copy of the VASCO catalog (see the section on obtaining it above), you can compare your results:

```
python3 poss_transients.py --plate 582 --vasco vasco_catalog.csv
```
The `--vasco` flag is optional. Without it, the pipeline runs the full extraction and crossmatch but skips the comparison step.

### Red Only Mode

```
python3 poss_transients.py --plate 582 --red-only
```
Skips the blue plate download and red-versus-blue comparison. Useful for quick testing or when the blue plate is not available.

### Custom Output Directory

```
python3 poss_transients.py --plate 582 --output-dir /data/poss_i/582
```

### With IR Crossmatch

```
python3 poss_transients.py --plate 582 --vasco vasco_catalog.csv --ir
```
Enables the infrared catalog crossmatch (AllWISE, 2MASS, CatWISE2020, UKIDSS). This is off by default because the published VASCO catalog was not IR-crossmatched. Turning it on will remove a small number of additional candidates (typically 1-2% of the total) that have infrared counterparts.

### Processing All 936 Plates

There is no built-in batch mode yet. To process every plate, you would write a wrapper script that loops over plate numbers. Using `--dedup-dir` pointed at a shared output directory handles the overlapping-plate problem automatically:

```bash
for num in $(seq 1 936); do
    python3 poss_transients.py --plate $num --output-dir ./output --dedup-dir ./output
done
```

Each plate's results are saved as `transients_{num}.csv` in the output directory. When the next plate runs, it reads all existing result files and removes any duplicate detections from overlapping plate edges.

**Be careful with this.** Processing all 936 plates will
- Download approximately 750 GB of plate scans (1,872 FITS files)
- Take roughly 35 to 45 minutes per plate, or about 550 to 700 hours total
- Generate hundreds of gigabytes of intermediate files and cached catalog queries
- Make thousands of VizieR API calls (VizieR does not have a formal rate limit, but sustained high-volume querying may result in temporary blocks)

A more practical approach is to process plates in batches of 10 to 20, verify the results look reasonable, and then scale up. Run multiple instances in parallel if your hardware supports it, but be mindful of RAM (each instance needs about 1 GB) and network bandwidth.

If you have access to a computing cluster or cloud instance, the pipeline is embarrassingly parallel across plates. Each plate is fully independent. You could process all 936 plates in a day with 30 to 40 concurrent workers on a machine with 64 GB of RAM and a fast network connection.

## Expected Performance

Based on testing with plate XE582 (166 published VASCO transients):

| Metric | Value |
|--------|-------|
| Sources extracted (red) | 74,330 |
| After morphological filters | 44,424 |
| After red-vs-blue comparison | 12,599 |
| After catalog crossmatch | 241 |
| VASCO detection rate | 99.4% (165/166 extracted) |
| VASCO recall after all filters | 86.7% (144/166) |
| Precision | 59.8% (144/241) |
| F1 score | 0.708 |

### Timing (single plate, SSD, 16 GB RAM)

| Step | Time |
|------|------|
| Plate download (400 MB) | 10 to 15 seconds |
| SExtractor pass 1 | 30 to 40 seconds |
| PSFEx | 1 to 12 minutes |
| SExtractor pass 2 | 10 to 20 minutes |
| Gaia query (or cache load) | 2 to 4 minutes |
| PanSTARRS query | 2 to 4 minutes |
| SuperCOSMOS query | 1 to 2 minutes |
| IR queries (if --ir enabled) | 3 to 6 minutes |
| Total per plate | 35 to 55 minutes |

## Known Limits

**SPREAD_MODEL calibration.** The PSFEx models built from photographic plate scans produce SPREAD_MODEL distributions that are shifted compared to models built from CCD data. The published Solano threshold of -0.002 rejects roughly half of known VASCO sources when applied to our PSF models. Used -0.01 as a compromise. This retains 97% of VASCO sources but allows more false positives through. The red-versus-blue comparison handles most of these.

**No asteroid removal.** Solano used the SkyBoT service (IMCCE) to identify known asteroids by checking if any minor planet was at each source's coordinates at the time of observation. This pipeline does not currently include SkyBoT queries. The published catalog flagged 189 sources as asteroids out of 298,165 initial detections, so the impact is small but not zero.

**No visual inspection.** Solano's team visually inspected a subset of candidates and reclassified 178 as artifacts. This pipeline does not include a visual inspection step. Automated morphological filters handle the majority of artifacts, but some will slip through.

**Blue plate availability.** Not every POSS-I field has a blue plate scan on IRSA. When the blue plate is missing, the pipeline falls back to red-only mode, which produces more false positives. The red-versus-blue comparison is the single most important filter for precision.

**Single-threaded extraction.** SExtractor and PSFEx are both single-threaded. The dominant time cost is SExtractor pass 2 with SPREAD_MODEL evaluation (10 to 20 minutes per plate). There is no way to speed this up on a single plate, but multiple plates can be processed in parallel.

**VizieR query volume.** Each plate queries roughly 2 to 4 million catalog rows across the primary catalogs (Gaia, PanSTARRS, SuperCOSMOS). With `--ir` enabled, this increases to 5 to 10 million rows across seven catalogs. VizieR is a shared public resource. If you are processing many plates in rapid succession, space out your runs or use the cached results to minimize repeated queries.

## File Structure

After processing a plate, the output directory looks like this:

```
output/
    dss1red_XE582.fits          # Downloaded red plate (400 MB)
    dss1blue_XO582.fits         # Downloaded blue plate (400+ MB)
    transients_582.csv          # Final transient candidates
    work_red/                   # SExtractor/PSFEx working files (red)
        sex.conf
        extract.param
        default.conv
        psfex.conf
        output.cat              # SExtractor catalog (FITS LDAC)
        output.psf              # PSFEx PSF model
    work_blue/                  # Same structure for blue plate
    cache/                      # Cached VizieR query results
        Gaia_EDR3.pkl
        Gaia_EDR3_full.pkl
        PanSTARRS_DR2.pkl
        SuperCOSMOS.pkl
        AllWISE.pkl             # Only present if --ir was used
        2MASS.pkl               # Only present if --ir was used
        CatWISE2020.pkl         # Only present if --ir was used
        UKIDSS-LAS.pkl          # Only present if --ir was used
```
Delete the `cache/` directory to force fresh catalog queries on the next run. Delete the FITS files to force fresh downloads. The working directories can be deleted after processing; they are not needed for the final results.

## Obtaining VASCO Catalog

The published VASCO transient catalog is freely available from the Spanish Virtual Observatory (SVO) archive, hosted by the Centro de Astrobiologia (CAB) at INTA in Madrid

http://svocats.cab.inta-csic.es/vanish/

This archive contains the full output of the Solano et al. (2022) pipeline

298,165 total sources detected only on POSS-I red plates (not in any modern optical survey)
Of those, 288,770 had crossmatches in other archives (mainly infrared)
5,399 remained as unidentified transient candidates after all filtering
172,163 were detected in the infrared but not in optical or ultraviolet

The archive provides a web-based query interface and supports Virtual Observatory (VO) standard access protocols. You can download the full catalog or query by position, plate ID, or other parameters.

For comparison against this pipeline, you need a CSV file with at minimum these columns

`Name`: Plate identifier ("XE582")
`RA`: Right ascension in decimal degrees (J2000)
`Dec`: Declination in decimal degrees (J2000)

The SVO archive lets you export in several formats including CSV and VOTable. If you are working with the broader set of approximately 107,000 transient candidates that includes sources with infrared crossmatches (used in some VASCO analyses beyond the Solano 2022 paper), those are also available through the same archive.

Not all 936 POSS-I survey plates have transients in the catalog. The published 5,399 unidentified transients span a subset of plates.

## References

Solano, E., Villarroel, B., and Rodrigo, C. (2022). "Discovering vanishing objects in POSS I red images using the Virtual Observatory." Monthly Notices of the Royal Astronomical Society, 515(1), 1380-1391.
Bertin, E. and Arnouts, S. (1996). "SExtractor: Software for source extraction." Astronomy and Astrophysics Supplement Series, 117, 393-404.
Bertin, E. (2011). "Automated Morphometry with SExtractor and PSFEx." Astronomical Data Analysis Software and Systems XX, ASP Conference Series, 442, 435.
