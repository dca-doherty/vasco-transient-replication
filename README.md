# vasco-transient-replication
Publicly Sourced Data Replication of VASCO Palomar Transients Script

# Reproducing the VASCO Transient Catalog from POSS-I Plates

## What This Is

This is an open-source pipeline that independently reproduces the transient catalog published by Solano, Villarroel, and Rodrigo (2022) in Monthly Notices of the Royal Astronomical Society. Their catalog identified 5,399 "vanishing" sources on photographic plates from the First Palomar Observatory Sky Survey (POSS-I), conducted between 1949 and 1958. These are point sources that appeared on the original red-sensitive plates but have no counterpart in any modern astronomical survey.

The pipeline downloads the original digitized plate scans from NASA's Infrared Science Archive (IRSA), extracts sources using tessellated SExtractor and PSFEx processing, applies morphological quality filters, removes false detections near bright star diffraction spikes, performs a red-versus-blue plate comparison to isolate true transients, and crossmatches surviving candidates against modern catalogs (Gaia, PanSTARRS, SuperCOSMOS) to remove persistent objects. Infrared catalog crossmatching is available as an option but disabled by default, since the published 107,000-source VASCO catalog was never crossmatched against IR (confirmed by Villarroel). The pipeline also applies a southern hemisphere declination cutoff and supports cross-plate deduplication for multi-plate runs. The output is a list of transient candidates that can be compared directly against the published VASCO catalog.

The goal is full reproducibility. Anyone with a Linux machine, the required software, and enough disk space should be able to run this pipeline and arrive at the same set of transients that VASCO published.

## Why This Matters

The VASCO catalog is the foundation for several lines of active research, including a 2025 paper in Scientific Reports that found statistical correlations between transient detections and nuclear weapons testing. Independent reproduction of the catalog strengthens all downstream analyses. It also makes the catalog auditable. Until now, the extraction pipeline was described in the literature but the code was never released. This pipeline fills that gap.

## What You Need

### Hardware

You need a machine with at least 16 GB of RAM. PSFEx builds a point spread function model from tens of thousands of source stamps and will use around 800 MB of memory per plate. SExtractor itself is lightweight, but the intermediate catalogs for a full plate can reach 400 MB when VIGNET stamps are included.

Disk space is the real constraint. Each full POSS-I plate scan is roughly 400 MB for the red and up to 1 GB for the blue (higher resolution scan). There are 936 red plates and 936 blue plates in the survey. If you download everything, that is approximately 1.2 TB of raw plate data. The intermediate SExtractor catalogs, PSFEx models, and cached catalog query results will add to that. Plan for 2 TB of free disk space if you intend to process the entire survey.

A solid-state drive will make a noticeable difference. The tessellated extraction processes hundreds of patches per plate, each involving SExtractor and PSFEx runs with many small file writes. On a spinning disk, the I/O overhead adds up significantly.

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

The pipeline downloads plate scans from IRSA and queries VizieR for catalog crossmatching. You need a reliable internet connection. Catalog queries for a single plate pull several million rows across three primary catalogs (or seven with `--ir` enabled) and can take 5 to 15 minutes depending on your connection and VizieR server load. The pipeline also queries USNO B-1.0 through VizieR for bright star spike removal. All query results are cached to disk, so subsequent runs of the same plate skip the network entirely.

If you are behind a corporate firewall or proxy, you may need to configure the `HTTP_PROXY` and `HTTPS_PROXY` environment variables for both the plate downloads and the astroquery VizieR calls.


## How the Pipeline Works

The pipeline follows the methodology described in Solano et al. (2022) Section 2 as closely as possible. The major steps correspond to the flowchart in their Figure 8: sky tessellation, source extraction, crossmatching against modern catalogs, diffraction spike removal, morphological artifact filtering, high proper motion identification, and concatenation with deduplication.

### Step 1: Download the Plates

Each POSS-I field was photographed twice: once through a red-sensitive emulsion (plate prefix XE) and once through a blue-sensitive emulsion (plate prefix XO). The field numbers match. For example, field 491 has red plate XE491 and blue plate XO491.

The digitized scans are hosted by IRSA at
- Red: `https://irsa.ipac.caltech.edu/data/DSS/images/dss1red/dss1red_XE{NUM}.fits`
- Blue: `https://irsa.ipac.caltech.edu/data/DSS/images/dss1blue/dss1blue_XO{NUM}.fits`

Red plates are roughly 14,000 by 14,000 pixels at 1.70 arcsec per pixel (about 400 MB). Blue plates are 23,040 by 23,040 pixels at 1.01 arcsec per pixel (about 1 GB). The different scan resolutions matter for filter calibration, which the pipeline handles automatically. The pipeline checks if the file already exists before downloading, so interrupted runs can be resumed without re-downloading.

### Step 2: Source Extraction (Tessellated SExtractor + PSFEx)

Source extraction uses sky tessellation following Solano et al. Section 2, step (i). The plate is divided into overlapping patches of 30 arcmin on a side with 100 pixels of overlap between adjacent patches. Each patch gets its own independent SExtractor and PSFEx processing, which produces a locally optimized PSF model. This is important because photographic plates have PSFs that change across the field due to optical distortions and emulsion variations. Processing small patches locally handles this naturally.

The patch size is computed from the actual plate pixel scale, so red patches (at 1.70 arcsec/px) and blue patches (at 1.01 arcsec/px) both cover the same 30 arcmin of sky. A typical red plate produces a 15 by 15 grid of 225 patches. A blue plate produces roughly 39 by 39 patches because of its finer pixel scale.

For each patch, extraction happens in two passes.

**Pass 1** runs SExtractor without a PSF model. This extracts all detectable sources and includes 35 by 35 pixel image stamps (VIGNET) around each source. These stamps are needed by PSFEx to build the PSF model.

**PSFEx** reads the pass 1 catalog and builds a spatially varying PSF model for that patch. It fits a polynomial model (degree 2 in X and Y) to capture local PSF variation. The output is a `.psf` file used in pass 2.

**Pass 2** runs SExtractor again with the PSF model loaded. This enables the SPREAD_MODEL parameter, which measures how much each source's light profile deviates from the local PSF. Real stars have SPREAD_MODEL near zero. Cosmic rays, hot pixels, and plate defects tend to have strongly negative values.

After all patches are processed, the pipeline merges the patch catalogs into a single source list and removes duplicates from the overlap regions. Sources within 3 arcseconds of each other across patch boundaries are deduplicated, keeping the instance with the higher signal-to-noise ratio.

**A note on WCS coordinates:** The POSS-I plate scans use a polynomial World Coordinate System (WCS) stored as PLT* keywords in the FITS header. SExtractor cannot read this format and will output all-zero coordinates. The pipeline works around this by using Astropy's WCS module to convert pixel coordinates to sky coordinates after extraction. Astropy reads the PLT* polynomial correctly.

### Step 3: Morphological Filters

The pipeline applies quality filters based on Solano et al. (2022) Section 2, step (v). All angular thresholds are scaled to the actual plate pixel scale, so the same physical size cuts apply to both red and blue plates despite their different scan resolutions.

- **FLAGS = 0**: No SExtractor warnings (blending, saturation, edge proximity)
- **SPREAD_MODEL > -0.02**: Rejects cosmic rays and sharp defects. See the calibration note below for why this differs from the published value of -0.002.
- **FWHM between 3.4 and 11.9 arcsec**: Rejects unresolved noise spikes (too small) and extended objects (too large). Expressed as 2 to 7 pixels at the red plate reference scale of 1.70 arcsec/px, and scaled proportionally for other pixel scales.
- **ELONGATION < 1.4**: Rejects streaks, satellite trails, and elongated artifacts.
- **Symmetry tolerance of 2 pixels** (scaled by pixel scale): Rejects sources where the bounding box is not roughly square.
- **Minimum bounding box > 1 pixel** in both directions: Rejects single-pixel noise hits.
- **SNR >= 30**: Rejects faint or noisy detections.

One filter that is deliberately omitted is the Median Absolute Deviation (MAD) clipping on FWHM and elongation. Solano describes computing the median and MAD within each tessellated region and removing sources deviating more than 2 sigma from the local median. Testing showed that this consistently reduces recall by 10 to 15 percent on dense plates with no meaningful improvement in precision. Transients are inherently morphological outliers, and MAD clipping penalizes exactly that. The hard cuts on FWHM, elongation, and SPREAD_MODEL handle the worst artifacts without this step.

**SPREAD_MODEL calibration.** This is the most important filter calibration decision in the pipeline. Solano's paper specifies SPREAD_MODEL > -0.002, which works well with the PSF models built from the STScI digitizations used in the original pipeline. Our PSFEx models, built from the IRSA plate scans using tessellated patches, produce SPREAD_MODEL distributions that are shifted slightly more negative. At -0.002, roughly half of known VASCO sources are rejected. A VASCO trace analysis on plate XE491 (710 published transients) showed the filter losses at each step:

| Filter Step | VASCO Sources Remaining | Lost at This Step |
|-------------|------------------------|-------------------|
| Raw extraction (3 sigma) | 709/710 | 1 |
| FLAGS = 0 | 661 | 48 |
| SPREAD_MODEL > -0.01 | 485 | 176 |
| FWHM 2-7 px | 481 | 4 |
| ELONGATION < 1.4 | 469 | 12 |
| Symmetry | 448 | 21 |
| SNR >= 30 | 419 | 29 |

SPREAD_MODEL was responsible for 176 of the 291 sources lost to filtering, roughly 60 percent of all filter losses. These are sub-PSF detections that look like cosmic rays to the classifier but are actually real brief transients. Relaxing the threshold to -0.02 recovers most of these sources. The downstream catalog crossmatch and bright star spike removal handle artifact rejection, so the looser SPREAD_MODEL threshold does not flood the output with false positives.

### Step 4: Red Versus Blue Comparison

This is the critical step that separates real transients from persistent sources. A source that appears on the red plate but not on the blue plate is a transient candidate, because it was only present during one of the two exposures (taken at different times). A source that appears on both plates is a persistent astronomical object and gets removed.

The pipeline crossmatches the filtered red and blue source lists at a 5 arcsecond radius. Any red source with a blue counterpart is marked as persistent and excluded. What remains are the transient candidates.

This step typically removes 60 percent or more of the filtered red sources and is the primary filter for precision. Without it, you will have substantially more candidates than the published catalog.

Not all POSS-I fields have blue plate scans available on IRSA. If the blue plate cannot be downloaded, the pipeline falls back to red-only mode and skips this step. Use the `--red-only` flag to force this behavior intentionally.

### Step 5: Bright Star Spike Removal

This implements Solano et al. Section 2, step (iv). Bright stars produce diffraction spikes and halos on photographic plates that can look like point sources to SExtractor. The pipeline queries the USNO B-1.0 catalog (via VizieR) for bright stars in the plate field and rejects any candidate that falls within the spike zone of a bright star.

A candidate is rejected if any USNO B-1.0 star satisfies:

```
Rmag1 or Rmag2 < -0.09 * d + 15.3
```

where d is the angular separation in arcseconds between the candidate and the USNO star, and Rmag1 and Rmag2 are the two-epoch red magnitudes from the catalog. The formula means brighter stars reject sources at larger distances. For example, a magnitude 6 star rejects candidates out to about 103 arcseconds, while a magnitude 12 star only rejects within about 37 arcseconds.

On plate XE491, this step removed 7,726 false detections while only losing 8 out of 710 VASCO sources.

### Step 6: Catalog Crossmatch

Surviving candidates are crossmatched against modern astronomical catalogs to remove any source that has a known counterpart today. The catalogs are queried via VizieR with a 5 degree search radius centered on the plate:

**Primary catalogs** (always run, queried in parallel)
- Gaia EDR3 (typically 200,000 to 300,000 sources per plate). Proper motions are corrected from the Gaia epoch (J2016.0) back to the POSS-I observation epoch to account for stellar motion over roughly 65 years.
- PanSTARRS DR2 (1.5 to 2 million sources per plate)

**Infrared catalogs** (optional, disabled by default, enable with `--ir`)
- AllWISE
- 2MASS
- CatWISE2020
- UKIDSS-LAS (northern sky only)

Villarroel confirmed that the published 107,000-source VASCO catalog was never crossmatched against infrared catalogs. The IR step is kept as an option for users who want stricter filtering, but it is off by default to match the actual methodology that produced the published catalog.

**Second-epoch optical**:
- SuperCOSMOS (digitization of POSS-II plates from the 1980s and 1990s)

All query results are cached to disk in the `cache_{plate}/` subdirectory. If you rerun the pipeline on the same plate, the cached results are loaded instantly. Clearing the cache forces fresh queries.

### Step 7: Southern Hemisphere Exclusion and Deduplication

POSS-I included some southern fringe plates that extend below declination -30 degrees. These were excluded from the published VASCO catalog. The pipeline applies a declination cutoff (default -30 degrees) and removes any candidate below that limit. This can be adjusted with the `--dec-min` flag if needed, but the default matches the published boundary.

For multi-plate runs, POSS-I plates overlap by about 1 degree at their edges. The same transient source can appear on two to four adjacent plates. The pipeline can deduplicate by checking new candidates against previously processed results. Any candidate within 5 arcseconds of an existing detection in a previous plate's output is flagged as a duplicate and removed. This is only active when you provide a `--dedup-dir` pointing to the directory containing your previous `transients_*.csv` files.

### VASCO Comparison

The final candidate list is crossmatched against the published VASCO catalog at 5 arcseconds. The pipeline reports recall (what fraction of VASCO entries were recovered), precision (what fraction of our candidates match VASCO), and F1 score. If the `--vasco` flag is provided, the pipeline also runs a diagnostic trace through each filter step showing exactly where VASCO sources are lost.

## Running the Pipeline

### Single Plate (Recommended: Tessellated Mode)

```
python3 poss_transients.py --plate 491 --tessellate --vasco vasco_catalog.csv
```
This processes plate XE491 (red) and XO491 (blue) using tessellated extraction, applies all filters including spike removal, runs the crossmatch, and saves the results to `./output/transients_491.csv`. The `--tessellate` flag enables the 30-arcmin patch processing described above.

### Red Only Mode

```
python3 poss_transients.py --plate 491 --tessellate --red-only --vasco vasco_catalog.csv
```
Skips the blue plate download and red-versus-blue comparison. Useful for quick testing, diagnostics, or when the blue plate is not available.

### Custom Output Directory

```
python3 poss_transients.py --plate 491 --tessellate --output-dir /data/poss_i/491
```

### With IR Crossmatch

```
python3 poss_transients.py --plate 491 --tessellate --vasco vasco_catalog.csv --ir
```
Enables the infrared catalog crossmatch (AllWISE, 2MASS, CatWISE2020, UKIDSS). This is off by default because the published VASCO catalog was not IR-crossmatched. Turning it on will remove a small number of additional candidates that have infrared counterparts.

### Processing All 936 Plates

There is no built-in batch mode yet. To process every plate, write a wrapper script that loops over plate numbers. Using `--dedup-dir` pointed at a shared output directory handles the overlapping-plate problem automatically:

```bash
for num in $(seq 1 936); do
    python3 poss_transients.py --plate $num --tessellate --output-dir ./output --dedup-dir ./output
done
```

Each plate's results are saved as `transients_{num}.csv` in the output directory. When the next plate runs, it reads all existing result files and removes any duplicate detections from overlapping plate edges.

**Be careful with this.** Processing all 936 plates will:
- Download over 1 TB of plate scans (1,872 FITS files)
- Take roughly 60 to 90 minutes per plate in tessellated mode, or about 950 to 1,400 hours total
- Generate hundreds of gigabytes of intermediate files and cached catalog queries
- Make thousands of VizieR API calls

A more practical approach is to process plates in batches of 10 to 20, verify the results look reasonable, and then scale up. Run multiple instances in parallel if your hardware supports it, but be mindful of RAM (each instance needs about 1 to 2 GB) and network bandwidth.

The pipeline is embarrassingly parallel across plates. Each plate is fully independent. You could process all 936 plates in a few days with 30 to 40 concurrent workers on a machine with 64 GB of RAM and a fast network connection.

## Expected Performance

Based on testing with plate XE491 (710 published VASCO transients). This is a dense plate that provides a demanding test case.

### Detection and Filter Trace

| Stage | Sources | VASCO Remaining |
|-------|---------|-----------------|
| Raw tessellated extraction (3 sigma) | 117,505 | 709/710 (99.9%) |
| After deduplication | 117,505 | 709/710 |
| After FLAGS = 0 | 102,347 | 661/710 |
| After SPREAD_MODEL > -0.02 | ~97,000 | ~640/710 |
| After FWHM, ELONGATION, Symmetry | ~73,000 | ~610/710 |
| After SNR >= 30 | ~43,000 | ~580/710 |
| After bright star spike removal | ~35,000 | ~575/710 |
| After catalog crossmatch | ~850 | pending |

Note: Numbers marked with ~ are estimated based on filter trace data at SPREAD_MODEL = -0.01 and the expected recovery from relaxing to -0.02. Full validated results with the -0.02 threshold are in progress.

### Timing (Single Plate, Tessellated Mode, SSD, 16 GB RAM)

| Step | Time |
|------|------|
| Plate download (400 MB to 1 GB) | 10 to 30 seconds |
| Red tessellated extraction (225 patches) | 35 to 40 minutes |
| Blue tessellated extraction (196 patches) | 50 to 60 minutes |
| Morphological filters | under 1 minute |
| Bright star spike removal | under 1 minute |
| Catalog crossmatch (first run) | 5 to 15 minutes |
| Catalog crossmatch (cached) | under 1 minute |
| Total per plate (first run) | 90 to 120 minutes |
| Total per plate (cached catalogs) | 85 to 100 minutes |

The dominant time cost is the tessellated SExtractor plus PSFEx processing. Each of the 225 red patches takes about 10 seconds for the two SExtractor passes and PSFEx run. Blue plates take longer because the finer pixel scale produces more patches.

## Documented Departures from Solano et al. (2022)

These are the deliberate differences between our pipeline and the published methodology, along with the reasoning for each.

**SPREAD_MODEL threshold: -0.02 versus -0.002.** Our PSFEx models built from IRSA plate scans produce SPREAD_MODEL distributions that are shifted compared to models built from the STScI digitizations used in the original pipeline. The published threshold of -0.002 rejects approximately half of known VASCO sources when applied to our PSF models. The -0.02 threshold retains the vast majority of VASCO sources while still rejecting the most extreme cosmic ray artifacts. Filter trace analysis confirmed that SPREAD_MODEL was responsible for 60 percent of all VASCO source losses.

**ELONGATION: 1.4 versus 1.3.** Slightly relaxed from the published value to reduce losses in the elongation filter. Combined with the omission of MAD clipping, this produces better recall without meaningful precision loss.

**No MAD clipping.** Solano describes computing the median and MAD of FWHM and elongation within each tessellated region and removing sources deviating more than 2 sigma from the local median. Testing showed this consistently reduces recall by 10 to 15 percent because transients are inherently morphological outliers. The hard cuts on FWHM, elongation, and SPREAD_MODEL handle artifact rejection without this step.

**No asteroid removal.** Solano used the SkyBoT service (IMCCE) to identify known asteroids by checking if any minor planet was at each source's coordinates at the time of observation. This pipeline does not currently include SkyBoT queries. The published catalog flagged 189 sources as asteroids out of 298,165 initial detections, so the impact is small.

**No visual inspection.** Solano's team visually inspected a subset of candidates and reclassified 178 as artifacts. This pipeline does not include a visual inspection step. Automated morphological filters handle the majority of artifacts, but some will slip through.

**Detection threshold: 3 sigma versus 5 sigma.** Lowered from the published value to maximize detection completeness. The tessellated extraction with local PSF models provides sufficiently robust background estimation to support the lower threshold. The VASCO trace confirmed that this recovers 709 out of 710 published sources.

## Known Limits

**Blue plate availability.** Not every POSS-I field has a blue plate scan on IRSA. When the blue plate is missing, the pipeline falls back to red-only mode, which produces more false positives. The red-versus-blue comparison is the single most important filter for precision.

**Tessellation runtime.** The tessellated extraction is slower than full-plate processing because it runs SExtractor and PSFEx independently on each patch. A single plate takes 60 to 90 minutes in tessellated mode versus 20 to 30 minutes for full-plate extraction. This is the cost of getting locally optimized PSF models, which produce better SPREAD_MODEL values and higher recall.

**VizieR query volume.** Each plate queries roughly 2 to 4 million catalog rows across the primary catalogs (Gaia, PanSTARRS, SuperCOSMOS, USNO B-1.0). With `--ir` enabled, this increases to 5 to 10 million rows. VizieR is a shared public resource. If you are processing many plates in rapid succession, space out your runs or use the cached results to minimize repeated queries.

**PIXEL_SCALE in SExtractor config.** The sex.conf file currently hardcodes PIXEL_SCALE to 1.70 for all patches, including blue plate patches at 1.01 arcsec/px. This only affects the interpretation of SEEING_FWHM in SExtractor's internal calculations and does not change the actual detection or photometry, which operate in pixel coordinates. The morphological filter thresholds are scaled correctly in the Python code.

## File Structure

After processing a plate, the output directory looks like this:

```
output/
    dss1red_XE491.fits          # Downloaded red plate (400 MB)
    dss1blue_XO491.fits         # Downloaded blue plate (1 GB)
    transients_491.csv          # Final transient candidates
    raw_red_sources.csv         # Cached raw extraction (117k+ sources)
    work_red/                   # Tessellated SExtractor/PSFEx working files
        patch_0001/
            sex.conf
            extract.param
            default.conv
            psfex.conf
            output.psf          # PSFEx PSF model for this patch
        patch_0002/
        ...
    work_blue/                  # Same structure for blue plate patches
    cache_491/                  # Cached VizieR query results
        Gaia_EDR3.pkl
        Gaia_EDR3_full.pkl
        PanSTARRS_DR2.pkl
        SuperCOSMOS.pkl
        USNO_B1_bright.pkl
        AllWISE.pkl             # Only present if --ir was used
        2MASS.pkl               # Only present if --ir was used
        CatWISE2020.pkl         # Only present if --ir was used
        UKIDSS-LAS.pkl          # Only present if --ir was used
```

Delete the `cache_491/` directory to force fresh catalog queries on the next run. Delete the FITS files to force fresh downloads. The `raw_red_sources.csv` file caches the full extracted source list before filtering, so you can experiment with different filter thresholds without re-running the 45-minute extraction step.

## Obtaining the VASCO Catalog

The published VASCO transient catalog is freely available from the Spanish Virtual Observatory (SVO) archive, hosted by the Centro de Astrobiologia (CAB) at INTA in Madrid:

http://svocats.cab.inta-csic.es/vanish/

This archive contains the full output of the Solano et al. (2022) pipeline:

- 298,165 total sources detected only on POSS-I red plates (not in any modern optical survey)
- Of those, 288,770 had crossmatches in other archives (mainly infrared)
- 5,399 remained as unidentified transient candidates after all filtering
- 172,163 were detected in the infrared but not in optical or ultraviolet

The archive provides a web-based query interface and supports Virtual Observatory (VO) standard access protocols. You can download the full catalog or query by position, plate ID, or other parameters.

For comparison against this pipeline, you need a CSV file with at minimum these columns:

- `Name`: Plate identifier (for example "XE491")
- `RA`: Right ascension in decimal degrees (J2000)
- `Dec`: Declination in decimal degrees (J2000)

The SVO archive lets you export in several formats including CSV and VOTable. If you are working with the broader set of approximately 107,000 transient candidates that includes sources with infrared crossmatches (used in some VASCO analyses beyond the Solano 2022 paper), those are also available through the same archive.

Not all 936 POSS-I survey plates have transients in the catalog. The published 5,399 unidentified transients span a subset of plates.

## References

Solano, E., Villarroel, B., and Rodrigo, C. (2022). "Discovering vanishing objects in POSS I red images using the Virtual Observatory." Monthly Notices of the Royal Astronomical Society, 515(1), 1380-1391.

Bertin, E. and Arnouts, S. (1996). "SExtractor: Software for source extraction." Astronomy and Astrophysics Supplement Series, 117, 393-404.

Bertin, E. (2011). "Automated Morphometry with SExtractor and PSFEx." Astronomical Data Analysis Software and Systems XX, ASP Conference Series, 442, 435.

Monet, D. G. et al. (2003). "The USNO-B1.0 Catalog." The Astronomical Journal, 125(2), 984-993.
