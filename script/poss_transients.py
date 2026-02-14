"""
poss_transients.py - Independent replication of the VASCO transient catalog.

Reproduces methodology from Solano, Villarroel & Rodrigo (2022, MNRAS 515, 1380) for detecting vanishing sources on First Palomar Observatory Sky Survey plates.
The basic idea... POSS-I photographed the sky twice per field, once through a red filter and once through blue. A real astronomical source shows up on both plates. Something
that only appears on one plate (and has no match w/modern catalog) is a transient. Solano's team found 298,165 of these via red plates
This script does the same thing from scratch using publicly available plate scans from IRSA and modern catalog data from VizieR. See the accompanying README for full doc.

IR crossmatch (AllWISE, 2MASS, CatWISE, UKIDSS) is OFF by default.
Villarroel confirmed the 107k VASCO catalog was never crossmatched
against IR catalogs. Use --ir to enable it.

Requirements
    - SExtractor (apt install source-extractor)
    - PSFEx (apt install psfex)
    - Python: numpy, pandas, astropy, astroquery

BD
"""
import argparse
import glob
import os
import pickle
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import astropy.units as u

# **** Where plates live IRSA ****
IRSA_URL = {
    "red":  "https://irsa.ipac.caltech.edu/data/DSS/images/dss1red",
    "blue": "https://irsa.ipac.caltech.edu/data/DSS/images/dss1blue",
}

# **** Extraction & filter params ****
# These follow Solano et al. (2022) Section 2 as closely as possible.
# One deliberate change is SPREAD_MODEL_MIN their paper uses -0.002 w/works w/their PSF models built from the STScI digitizations.
# built own PSFEx models from the IRSA plate scans, the SPREAD_MODEL distribution shifts slightly negative. Using -0.002 kills
# about half known VASCO sources. At -0.01 I kept 97% of them while still rejecting obv cosmic rays. Red v blue compare handles the rest.

DETECT_THRESH    = 5       # sigma above background
SNR_MIN          = 30      # min signal-to-noise
SPREAD_MODEL_MIN = -0.01   # cosmic ray rejection (note above)
FWHM_RANGE       = (2, 7)  # pixels rejects noise spikes + junk
ELONG_MAX        = 1.4     # rejects streaks (relaxed from 1.3; recovers borderline VASCO sources)
SYMMETRY_TOL     = 2       # pixels bounding box dx v dy
MATCH_RADIUS     = 5.0     # arcseconds for all crossmatching
PLATE_SCALE      = 1.7     # arcsec/pixel for the POSS-I scans

# plate dl

def download_plate(plate_num, band, output_dir):
    """
    Grab POSS-I plate from IRSA > Returns the local path... or None if it
    fails. Skips dl if file already exists
    """
    prefix = "XE" if band == "red" else "XO"
    subdir = "dss1red" if band == "red" else "dss1blue"
    filename = f"{subdir}_{prefix}{plate_num}.fits"
    url = f"{IRSA_URL[band]}/{filename}"
    local_path = os.path.join(output_dir, filename)

    if os.path.exists(local_path) and os.path.getsize(local_path) > 1e6:
        mb = os.path.getsize(local_path) / 1e6
        print(f"  [{band}] Using cached plate: {local_path} ({mb:.0f} MB)")
        return local_path

    print(f"  [{band}] Downloading {url}")
    try:
        t0 = time.time()
        urllib.request.urlretrieve(url, local_path)
        mb = os.path.getsize(local_path) / 1e6
        print(f"  [{band}] Got {mb:.0f} MB in {time.time() - t0:.0f}s")
    except Exception as err:
        print(f"  [{band}] Download failed: {err}")
        return None

    return local_path

#  SExtractor + PSFEx

def _find_sextractor():
    """SExtractor ships under different names depending on package."""
    for name in ["source-extractor", "sextractor", "sex"]:
        if subprocess.run(["which", name], capture_output=True).returncode == 0:
            return name
    return None


def _write_sextractor_config(work_dir, psf_path=None):
    """
    Generate config files SExtractor needs. There are 3 of them
    the main .conf, parameter list, and convolution kernel.

    If psf_path is given, I included SPREAD_MODEL in output/skip
    VIGNET w/is only needed for building the PSF first place.
    """
    work_dir = os.path.abspath(work_dir)

    # Convolution kernel 3x3 Gaussian-ish
    conv_path = os.path.join(work_dir, "default.conv")
    with open(conv_path, "w") as f:
        f.write("CONV NORM\n1 2 1\n2 4 2\n1 2 1\n")

    # PSFEx config only used during PSF building step
    psfex_path = os.path.join(work_dir, "psfex.conf")
    with open(psfex_path, "w") as f:
        f.write("BASIS_TYPE       PIXEL_AUTO\n")
        f.write("PSF_SAMPLING     0.0\n")
        f.write("PSF_ACCURACY     0.01\n")
        f.write("PSF_SIZE         25,25\n")
        f.write("PSFVAR_KEYS      X_IMAGE,Y_IMAGE\n")
        f.write("PSFVAR_GROUPS    1,1\n")
        f.write("PSFVAR_DEGREES   2\n")
        f.write("SAMPLE_FWHMRANGE 2.0,20.0\n")
        f.write("SAMPLE_VARIABILITY 0.3\n")
        f.write("SAMPLE_MINSN     20\n")
        f.write("SAMPLE_MAXSOURCES   5000\n")
        f.write("CHECKIMAGE_TYPE  NONE\n")
        f.write("CHECKPLOT_TYPE   NONE\n")
        f.write("WRITE_XML        N\n")
    # Output params
    params = [
        "NUMBER", "ALPHA_J2000", "DELTA_J2000",
        "X_IMAGE", "Y_IMAGE",
        "FLUX_AUTO", "FLUXERR_AUTO",
        "FLUX_APER", "FLUXERR_APER",
        "MAG_AUTO", "MAGERR_AUTO",
        "FLUX_RADIUS", "FWHM_IMAGE",
        "ELONGATION", "ELLIPTICITY",
        "FLAGS", "SNR_WIN",
        "XMIN_IMAGE", "XMAX_IMAGE", "YMIN_IMAGE", "YMAX_IMAGE",
        "ISOAREA_IMAGE", "BACKGROUND",
    ]
    if psf_path:
        params += ["SPREAD_MODEL", "SPREADERR_MODEL"]
    else:
        # PSFEx needs 35x35 pixel stamps to build model
        # Makes catalog huge approx. 400mb... so I only include in the first pass
        params += ["VIGNET(35,35)"]

    param_path = os.path.join(work_dir, "extract.param")
    with open(param_path, "w") as f:
        f.write("\n".join(params) + "\n")

    # Main SExtractor config
    cat_path = os.path.join(work_dir, "output.cat")
    conf_path = os.path.join(work_dir, "sex.conf")
    settings = {
        "CATALOG_NAME":     cat_path,
        "CATALOG_TYPE":     "FITS_LDAC",
        "PARAMETERS_NAME":  param_path,
        "DETECT_TYPE":      "CCD",
        "DETECT_MINAREA":   "5",
        "DETECT_THRESH":    str(DETECT_THRESH),
        "ANALYSIS_THRESH":  str(DETECT_THRESH),
        "FILTER":           "Y",
        "FILTER_NAME":      conv_path,
        "DEBLEND_NTHRESH":  "32",
        "DEBLEND_MINCONT":  "0.005",
        "CLEAN":            "Y",
        "CLEAN_PARAM":      "1.0",
        "WEIGHT_TYPE":      "BACKGROUND",
        "PHOT_APERTURES":   "5",
        "PHOT_AUTOPARAMS":  "2.5,3.5",
        "PHOT_PETROPARAMS": "2.0,3.5",
        "SATUR_LEVEL":      "50000.0",
        "MAG_ZEROPOINT":    "25.0",
        "GAIN":             "1.0",
        "PIXEL_SCALE":      str(PLATE_SCALE),
        "SEEING_FWHM":      "3.0",
        "BACK_SIZE":        "256",
        "BACK_FILTERSIZE":  "3",
        "BACKPHOTO_TYPE":   "LOCAL",
    }
    if psf_path:
        settings["PSF_NAME"] = psf_path

    with open(conf_path, "w") as f:
        for key, val in settings.items():
            f.write(f"{key:<20s} {val}\n")

    return conf_path, psfex_path, cat_path


def extract_sources(fits_path, work_dir, label="plate"):
    """
    Two-pass SExtractor + PSFEx extraction.

    Pass 1 runs without a PSF model and includes VIGNET stamps PSFEx
    needs to build its spatially-varying PSF. Pass 2 loads that PSF &
    gives us SPREAD_MODEL, w/separates real sources from cosmic rays + plate defects.

    The POSS-I plates use a polynomial WCS (PLT* keywords) that SExtractor
    can't parse, so all RA/Dec come out as 0. I fixed this after
    extraction using astropy's WCS, w/handles the polynomial fine.
    """
    work_dir = os.path.abspath(work_dir)
    fits_path = os.path.abspath(fits_path)
    os.makedirs(work_dir, exist_ok=True)

    # Clean up files from prev run
    for pattern in ["*.conf", "*.param", "*.cat", "*.psf", "*.nnw", "*.xml"]:
        for old_file in glob.glob(os.path.join(work_dir, pattern)):
            os.remove(old_file)

    sex_cmd = _find_sextractor()
    if not sex_cmd:
        print(f"  [{label}] SExtractor not found. Install with: apt install source-extractor")
        sys.exit(1)

    # ** Pass 1 basic extraction + VIGNET stamps PSFEx **
    conf, psfex_conf, cat_path = _write_sextractor_config(work_dir)
    print(f"  [{label}] Running SExtractor (pass 1, with VIGNET stamps)...")
    t0 = time.time()
    subprocess.run([sex_cmd, fits_path, "-c", conf], capture_output=True, text=True)

    if not os.path.exists(cat_path):
        print(f"  [{label}] SExtractor failed on pass 1")
        sys.exit(1)
    print(f"  [{label}] Pass 1 done: {time.time()-t0:.0f}s, "
          f"{os.path.getsize(cat_path)/1e6:.0f} MB catalog")

    # **PSF model **
    psf_path = os.path.join(work_dir, "output.psf")
    print(f"  [{label}] Building PSF model with PSFEx...")
    t0 = time.time()
    try:
        subprocess.run(
            ["psfex", cat_path, "-c", psfex_conf, "-WRITE_XML", "N"],
            capture_output=True, text=True, cwd=work_dir, timeout=2400,
        )
    except subprocess.TimeoutExpired:
        print(f"  [{label}] PSFEx timed out after 40 min")

    psf_ok = os.path.exists(psf_path)
    if psf_ok:
        print(f"  [{label}] PSF model built in {time.time()-t0:.0f}s")
    else:
        print(f"  [{label}] PSFEx failed; continuing without SPREAD_MODEL")

    # ** Pass 2 re-extract with PSF model SPREAD_MODEL **
    if psf_ok:
        conf2, _, cat_path = _write_sextractor_config(work_dir, psf_path=psf_path)
        print(f"  [{label}] Running SExtractor (pass 2, with SPREAD_MODEL)...")
        t0 = time.time()
        subprocess.run([sex_cmd, fits_path, "-c", conf2], capture_output=True, text=True)

        if os.path.exists(cat_path) and os.path.getsize(cat_path) > 1000:
            print(f"  [{label}] Pass 2 done: {time.time()-t0:.0f}s")
        else:
            print(f"  [{label}] Pass 2 failed, using pass 1 catalog")

    # ** Read catalog into a df **
    # VIGNET col is multidim and chokes pandas, so filter down scalar col only
    tbl = Table.read(cat_path, hdu=2)
    scalar_cols = [c for c in tbl.colnames if len(tbl[c].shape) <= 1]
    df = tbl[scalar_cols].to_pandas()
    print(f"  [{label}] Extracted {len(df)} sources")

    # **Fix WCS coordinates **
    # SExtractor outputs all 0 for RA/Dec because it can't read
    # DSS polynomial WCS Astropy works though
    header = fits.getheader(fits_path)
    wcs = WCS(header)
    ra, dec = wcs.all_pix2world(df["X_IMAGE"].values, df["Y_IMAGE"].values, 1)
    df["ALPHA_J2000"] = ra
    df["DELTA_J2000"] = dec
    print(f"  [{label}] Sky coordinates: "
          f"RA [{ra.min():.2f}, {ra.max():.2f}], "
          f"Dec [{dec.min():.2f}, {dec.max():.2f}]")

    return df

#  Morphological quality filters
def apply_quality_filters(df, label="plate"):
    """
    Apply the same source quality cuts described in Solano et al. (2022).

    Includes MAD (median absolute deviation) clipping on FWHM and
    elongation as described in the paper: sources deviating more than
    2*MAD from the median are removed.
    """
    n_start = len(df)

    # No extraction warnings
    df = df[df["FLAGS"] == 0].copy()
    print(f"  [{label}] FLAGS = 0: {n_start} -> {len(df)}")

    # Reject cosmic rays +sharp artifacts
    if "SPREAD_MODEL" in df.columns:
        n = len(df)
        df = df[df["SPREAD_MODEL"] > SPREAD_MODEL_MIN].copy()
        print(f"  [{label}] SPREAD_MODEL > {SPREAD_MODEL_MIN}: {n} -> {len(df)}")

    # size cut
    n = len(df)
    lo, hi = FWHM_RANGE
    df = df[(df["FWHM_IMAGE"] >= lo) & (df["FWHM_IMAGE"] <= hi)].copy()
    print(f"  [{label}] FWHM [{lo}, {hi}] px: {n} -> {len(df)}")

    # Roundness
    n = len(df)
    df = df[df["ELONGATION"] < ELONG_MAX].copy()
    print(f"  [{label}] ELONGATION < {ELONG_MAX}: {n} -> {len(df)}")

    # Symmetry bounding box should be square
    n = len(df)
    dx = df["XMAX_IMAGE"] - df["XMIN_IMAGE"]
    dy = df["YMAX_IMAGE"] - df["YMIN_IMAGE"]
    df = df[np.abs(dx - dy) <= SYMMETRY_TOL].copy()
    print(f"  [{label}] Symmetry <= {SYMMETRY_TOL} px: {n} -> {len(df)}")

    # Min bounding box (Solano: XMAX-XMIN > 1 and YMAX-YMIN > 1)
    n = len(df)
    dx = df["XMAX_IMAGE"] - df["XMIN_IMAGE"]
    dy = df["YMAX_IMAGE"] - df["YMIN_IMAGE"]
    df = df[(dx > 1) & (dy > 1)].copy()
    print(f"  [{label}] Min bbox > 1 px: {n} -> {len(df)}")

    # Signal 2 noise
    n = len(df)
    df = df[df["SNR_WIN"] >= SNR_MIN].copy()
    print(f"  [{label}] SNR >= {SNR_MIN}: {n} -> {len(df)}")

    # MAD clipping on FWHM and elongation (Solano: 2-sigma from median)
    n = len(df)
    for col in ["FWHM_IMAGE", "ELONGATION"]:
        med = df[col].median()
        mad = np.median(np.abs(df[col] - med))
        sigma = 1.4826 * mad  # MAD to std deviation for normal distribution
        lo_clip = med - 2 * sigma
        hi_clip = med + 2 * sigma
        before = len(df)
        df = df[(df[col] >= lo_clip) & (df[col] <= hi_clip)].copy()
        print(f"  [{label}] MAD clip {col}: median={med:.3f}, "
              f"sigma={sigma:.3f}, range=[{lo_clip:.3f}, {hi_clip:.3f}], "
              f"{before} -> {len(df)}")
    print(f"  [{label}] After MAD clipping: {n} -> {len(df)}")

    pct = 100 * len(df) / max(n_start, 1)
    print(f"  [{label}] Kept {len(df)}/{n_start} sources ({pct:.1f}%)")
    return df

#  Red v blue plate compare

def find_transients(red_df, blue_df):
    """
    Core of transient detection... find sources that appear on the
    red plate but NOT on the blue plate. Anything present on both plates
    is an astronomical source & gets removed.

    Red + blue exposures were taken at different times and have
    different exposure durations, so a true transient caught during one exposure won't be on the other.
    """
    print(f"\n  Red sources: {len(red_df)}, Blue sources: {len(blue_df)}")

    red_sky = SkyCoord(ra=red_df["ALPHA_J2000"].values,
                       dec=red_df["DELTA_J2000"].values, unit="deg")
    blue_sky = SkyCoord(ra=blue_df["ALPHA_J2000"].values,
                        dec=blue_df["DELTA_J2000"].values, unit="deg")

    idx_matched, _, _, _ = search_around_sky(
        red_sky, blue_sky, MATCH_RADIUS * u.arcsec)

    persistent = set(idx_matched)
    is_transient = np.ones(len(red_df), dtype=bool)
    for i in persistent:
        is_transient[i] = False

    n_persistent = len(persistent)
    n_transient = is_transient.sum()
    print(f"  Matched on both plates: {n_persistent}")
    print(f"  Red-only (transient candidates): {n_transient}")

    return red_df[is_transient].copy()

#  Catalog crossmatch

def _query_catalog(catalog_id, columns, ra_center, dec_center, radius_deg,
                   ra_key, dec_key, label, cache_dir):
    """
    Query VizieR catalog + cache result to disk. subsequent runs loads from cache instead of hitting network again.
    Returns SkyCoord array, or none if query fails.
    """
    from astroquery.vizier import Vizier

    safe_label = label.replace(" ", "_").replace("/", "_")
    cache_file = os.path.join(cache_dir, f"{safe_label}.pkl")

    if os.path.exists(cache_file):
        with open(cache_file, "rb") as f:
            coords = pickle.load(f)
        n = len(coords) if coords is not None else 0
        print(f"  {label}: {n} sources (cached)")
        return coords

    try:
        viz = Vizier(columns=columns, row_limit=-1, timeout=300)
        center = SkyCoord(ra=ra_center, dec=dec_center, unit="deg")
        result = viz.query_region(center, radius=radius_deg * u.deg,
                                  catalog=catalog_id)
        if result and len(result[0]) > 0:
            tbl = result[0]
            ra = np.array(tbl[ra_key], dtype=float)
            dec = np.array(tbl[dec_key], dtype=float)
            good = np.isfinite(ra) & np.isfinite(dec)
            coords = SkyCoord(ra=ra[good], dec=dec[good], unit="deg")
            print(f"  {label}: {len(coords)} sources")
        else:
            coords = None
            print(f"  {label}: empty result")
    except Exception as err:
        coords = None
        print(f"  {label}: query failed ({err})")

    with open(cache_file, "wb") as f:
        pickle.dump(coords, f)
    return coords


def crossmatch_modern_catalogs(df, cache_dir, plate_epoch=1953.0, use_ir=False):
    """
    Check source against modern catalogs. Anything w/match
    is a known object & gets flagged for removal.

    Primary (always):  Gaia EDR3, PanSTARRS DR2
    Infrared (--ir):   AllWISE, 2MASS, CatWISE2020, UKIDSS-LAS
    Second epoch:      SuperCOSMOS (POSS-II digitization)

    IR catalogs are off by default because the published 107k VASCO
    catalog was never IR-crossmatched (confirmed by Villarroel).

    Gaia proper motions are corrected back to the POSS-I epoch approx. 1951
    so that high proper motion stars still get matched.
    """
    src_coords = SkyCoord(ra=df["ALPHA_J2000"].values,
                          dec=df["DELTA_J2000"].values, unit="deg")
    ra_center = df["ALPHA_J2000"].median()
    dec_center = df["DELTA_J2000"].median()
    search_radius = 5.0  # degrees full plate is about 6.5 deg across
    matched = np.zeros(len(df), dtype=bool)

    os.makedirs(cache_dir, exist_ok=True)

    # Gaia + PanSTARRS primary optical catalogs query in parallel

    print(f"\n  Crossmatching against primary catalogs...")

    def fetch_gaia():
        return _query_catalog(
            "I/350/gaiaedr3", ["RA_ICRS", "DE_ICRS", "pmRA", "pmDE"],
            ra_center, dec_center, search_radius,
            "RA_ICRS", "DE_ICRS", "Gaia EDR3", cache_dir)

    def fetch_panstarrs():
        return _query_catalog(
            "II/349/ps1", ["RAJ2000", "DEJ2000"],
            ra_center, dec_center, search_radius,
            "RAJ2000", "DEJ2000", "PanSTARRS DR2", cache_dir)

    with ThreadPoolExecutor(max_workers=2) as pool:
        gaia_future = pool.submit(fetch_gaia)
        ps_future = pool.submit(fetch_panstarrs)
        gaia_result = gaia_future.result()
        ps_result = ps_future.result()

    # Gaia needs proper motion correction catalog epoch is 2016
    # POSS-I plates are from the early 1950s... so... stars have obviously moved around since then..
    if gaia_result is not None:
        gaia_cache = os.path.join(cache_dir, "Gaia_EDR3_full.pkl")
        if not os.path.exists(gaia_cache):
            from astroquery.vizier import Vizier
            viz = Vizier(columns=["RA_ICRS", "DE_ICRS", "pmRA", "pmDE"],
                         row_limit=-1, timeout=300)
            center = SkyCoord(ra=ra_center, dec=dec_center, unit="deg")
            r = viz.query_region(center, radius=search_radius * u.deg,
                                 catalog="I/350/gaiaedr3")
            gaia_tbl = r[0] if r else None
            with open(gaia_cache, "wb") as f:
                pickle.dump(gaia_tbl, f)
        else:
            with open(gaia_cache, "rb") as f:
                gaia_tbl = pickle.load(f)

        if gaia_tbl is not None:
            gra = np.array(gaia_tbl["RA_ICRS"], dtype=float)
            gdec = np.array(gaia_tbl["DE_ICRS"], dtype=float)
            pm_ra = np.nan_to_num(np.array(gaia_tbl["pmRA"], dtype=float), 0)
            pm_dec = np.nan_to_num(np.array(gaia_tbl["pmDE"], dtype=float), 0)

            # Roll pos back from 2016 to plate epoch
            dt_years = plate_epoch - 2016.0
            gra += pm_ra * dt_years / 3600000.0 / np.cos(np.radians(gdec))
            gdec += pm_dec * dt_years / 3600000.0

            good = np.isfinite(gra) & np.isfinite(gdec)
            gaia_sky = SkyCoord(ra=gra[good], dec=gdec[good], unit="deg")
            idx, _, _, _ = search_around_sky(src_coords, gaia_sky,
                                             MATCH_RADIUS * u.arcsec)
            matched[idx] = True
            print(f"  Gaia: {len(set(idx))} matches")

    if ps_result is not None:
        idx, _, _, _ = search_around_sky(src_coords, ps_result,
                                         MATCH_RADIUS * u.arcsec)
        new_matches = len(set(idx[~matched[idx]]))
        matched[idx] = True
        print(f"  PanSTARRS: {new_matches} new matches")

    n_remaining = (~matched).sum()
    print(f"  After primary catalogs: {n_remaining} candidates remaining")

    if n_remaining == 0:
        df["is_known"] = matched
        return df

    # **** Infrared catalogs candidates only query parallel ****
    # OFF by default. The 107k VASCO catalog was never IR-matched
    # (confirmed by Villarroel). Enable with --ir flag.

    if use_ir:
        print(f"\n  Checking {n_remaining} candidates against IR catalogs...")
        cand_idx = np.where(~matched)[0]
        cand_sky = SkyCoord(ra=df.iloc[cand_idx]["ALPHA_J2000"].values,
                            dec=df.iloc[cand_idx]["DELTA_J2000"].values, unit="deg")

        ir_catalogs = [
            ("II/328/allwise", ["RAJ2000", "DEJ2000"], "RAJ2000", "DEJ2000", "AllWISE"),
            ("II/246/out",     ["RAJ2000", "DEJ2000"], "RAJ2000", "DEJ2000", "2MASS"),
            ("II/365/catwise", ["RA_ICRS", "DE_ICRS"], "RA_ICRS", "DE_ICRS", "CatWISE2020"),
            ("II/319/las9",    ["RAJ2000", "DEJ2000"], "RAJ2000", "DEJ2000", "UKIDSS-LAS"),
        ]

        ir_matched = np.zeros(len(cand_idx), dtype=bool)

        def fetch_ir(cat_id, cols, rk, dk, name):
            coords = _query_catalog(cat_id, cols, ra_center, dec_center,
                                    search_radius, rk, dk, name, cache_dir)
            return name, coords

        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = [pool.submit(fetch_ir, *cat) for cat in ir_catalogs]
            for fut in as_completed(futures):
                name, coords = fut.result()
                if coords is not None:
                    idx, _, _, _ = search_around_sky(cand_sky, coords,
                                                     MATCH_RADIUS * u.arcsec)
                    new = len(set(idx[~ir_matched[idx]]))
                    ir_matched[idx] = True
                    print(f"  {name}: {new} new matches")

        for i, ci in enumerate(cand_idx):
            if ir_matched[i]:
                matched[ci] = True
        print(f"  IR catalogs removed: {ir_matched.sum()} total")

        n_remaining = (~matched).sum()
    else:
        print(f"\n  IR crossmatch: skipped (use --ir to enable)")

    # **** SuperCOSMOS 2nd epoch check ****
    
    n_remaining = (~matched).sum()
    if n_remaining > 0:
        print(f"\n  Checking {n_remaining} candidates against SuperCOSMOS...")
        cand_idx2 = np.where(~matched)[0]
        cand_sky2 = SkyCoord(ra=df.iloc[cand_idx2]["ALPHA_J2000"].values,
                             dec=df.iloc[cand_idx2]["DELTA_J2000"].values, unit="deg")

        scos = _query_catalog(
            "I/305/out", ["RAJ2000", "DEJ2000"], ra_center, dec_center,
            search_radius, "RAJ2000", "DEJ2000", "SuperCOSMOS", cache_dir)

        if scos is not None:
            idx, _, _, _ = search_around_sky(cand_sky2, scos,
                                             MATCH_RADIUS * u.arcsec)
            for i in set(idx):
                matched[cand_idx2[i]] = True
            print(f"  SuperCOSMOS: removed {len(set(idx))}")

    final_count = (~matched).sum()
    print(f"\n  Final transient candidates: {final_count}")

    df["is_known"] = matched
    return df

#  VASCO catalog compare

def compare_to_vasco(candidates, vasco_csv, plate_num):
    """
    Match our candidates against the published VASCO catalog and report
    precision/recall. This is the validation step.
    """
    vasco = pd.read_csv(vasco_csv)
    plate_id = f"XE{plate_num}"
    vasco_plate = vasco[vasco["Name"] == plate_id]

    print(f"\n{'=' * 60}")
    print(f"  VASCO COMPARISON - Plate {plate_id}")
    print(f"{'=' * 60}")
    print(f"  Published VASCO entries: {len(vasco_plate)}")
    print(f"  Our candidates:         {len(candidates)}")

    if len(vasco_plate) == 0 or len(candidates) == 0:
        print("  (nothing to compare)")
        return

    our_sky = SkyCoord(ra=candidates["ALPHA_J2000"].values,
                       dec=candidates["DELTA_J2000"].values, unit="deg")
    vasco_sky = SkyCoord(ra=vasco_plate["RA"].values,
                         dec=vasco_plate["Dec"].values, unit="deg")

    idx_ours, idx_vasco, _, _ = search_around_sky(
        our_sky, vasco_sky, MATCH_RADIUS * u.arcsec)

    n_our_matched = len(set(idx_ours))
    n_vasco_found = len(set(idx_vasco))
    n_vasco_total = len(vasco_plate)
    n_ours_total = len(candidates)

    recall = n_vasco_found / max(n_vasco_total, 1)
    precision = n_our_matched / max(n_ours_total, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-10)

    print(f"  Recovered:  {n_vasco_found}/{n_vasco_total} "
          f"({100 * recall:.1f}% recall)")
    print(f"  Precision:  {n_our_matched}/{n_ours_total} "
          f"({100 * precision:.1f}%)")
    print(f"  F1:         {f1:.3f}")
    print(f"  Extras:     {n_ours_total - n_our_matched}")
    print(f"  Missed:     {n_vasco_total - n_vasco_found}")

#  Southern hemisphere exclusion

def apply_dec_cutoff(df, dec_min=-30.0):
    """
    Remove sources below a declination limit. POSS-I southern fringe
    plates (roughly dec < -30) were excluded from the published VASCO
    catalog. This matches that boundary.
    """
    n_before = len(df)
    df = df[df["DELTA_J2000"] >= dec_min].copy()
    n_removed = n_before - len(df)
    if n_removed > 0:
        print(f"  Dec cutoff (>= {dec_min}): removed {n_removed}, "
              f"{len(df)} remain")
    return df


#  Cross-plate deduplication

def deduplicate_across_plates(df, prev_results_dir, radius_arcsec=5.0):
    """
    For multi-plate runs: flag candidates that were already detected
    on a previously processed overlapping plate. POSS-I plates overlap
    by ~1 degree at edges, so the same transient can appear on 2-4 plates.

    Reads all transients_*.csv files in prev_results_dir, builds a master
    coordinate list, and marks any new candidate within radius as a dupe.
    """
    import glob as _glob

    prev_files = _glob.glob(os.path.join(prev_results_dir, "transients_*.csv"))
    if not prev_files:
        print(f"  Dedup: no previous results in {prev_results_dir}")
        return df

    prev_frames = []
    for f in prev_files:
        try:
            prev_frames.append(pd.read_csv(f))
        except Exception:
            pass

    if not prev_frames:
        return df

    prev_all = pd.concat(prev_frames, ignore_index=True)
    if "ALPHA_J2000" not in prev_all.columns:
        print(f"  Dedup: previous CSVs missing coordinate columns")
        return df

    prev_sky = SkyCoord(ra=prev_all["ALPHA_J2000"].values,
                        dec=prev_all["DELTA_J2000"].values, unit="deg")
    new_sky = SkyCoord(ra=df["ALPHA_J2000"].values,
                       dec=df["DELTA_J2000"].values, unit="deg")

    idx_new, _, _, _ = search_around_sky(new_sky, prev_sky,
                                          radius_arcsec * u.arcsec)
    dupes = set(idx_new)
    n_dupes = len(dupes)

    if n_dupes > 0:
        keep = np.ones(len(df), dtype=bool)
        for i in dupes:
            keep[i] = False
        df = df[keep].copy()
        print(f"  Dedup: {n_dupes} duplicates from overlapping plates removed, "
              f"{len(df)} unique remain")
    else:
        print(f"  Dedup: no duplicates found (checked against "
              f"{len(prev_all)} previous detections)")

    return df


#  Main

def main():
    parser = argparse.ArgumentParser(
        description="Replicate the VASCO transient catalog from POSS-I plates.")
    parser.add_argument("--plate", required=True,
                        help="Plate field number, e.g. 582")
    parser.add_argument("--vasco",
                        help="Path to VASCO catalog CSV for validation "
                             "(optional, skips comparison if not provided)")
    parser.add_argument("--output-dir", default="./output",
                        help="Where to put everything (default: ./output)")
    parser.add_argument("--red-only", action="store_true",
                        help="Skip blue plate download and red-vs-blue step")
    parser.add_argument("--ir", action="store_true",
                        help="Enable IR crossmatch (AllWISE, 2MASS, CatWISE, "
                             "UKIDSS). Off by default -- the published 107k "
                             "VASCO catalog was not IR-matched (confirmed by "
                             "Villarroel).")
    parser.add_argument("--dec-min", type=float, default=-30.0,
                        help="Minimum declination cutoff (default: -30). "
                             "POSS-I southern fringe plates below this were "
                             "excluded from the published catalog.")
    parser.add_argument("--dedup-dir",
                        help="Path to directory containing previous plate "
                             "results CSVs. New candidates within 5 arcsec "
                             "of existing detections are flagged as duplicates "
                             "from overlapping plates.")
    args = parser.parse_args()

    plate = args.plate
    outdir = args.output_dir
    os.makedirs(outdir, exist_ok=True)
    cache_dir = os.path.join(outdir, "cache")

    print("=" * 60)
    print(f"  POSS-I Transient Detection Pipeline")
    print(f"  Plate {plate}  (red = XE{plate}, blue = XO{plate})")
    print(f"  Mode: {'red only' if args.red_only else 'red + blue comparison'}")
    print(f"  IR crossmatch: {'enabled' if args.ir else 'disabled (use --ir)'}")
    print("=" * 60)

    # Step 1 get the plates... and wait
    print(f"\n[1/6] Downloading plates...")
    red_fits = download_plate(plate, "red", outdir)
    if not red_fits:
        print("Could not get the red plate. Stopping.")
        sys.exit(1)

    blue_fits = None
    if not args.red_only:
        blue_fits = download_plate(plate, "blue", outdir)
        if not blue_fits:
            print("Blue plate not available on IRSA. Continuing with red only.")

    with fits.open(red_fits) as hdu:
        h = hdu[0].header
        print(f"  Plate dimensions: {h.get('NAXIS1', '?')} x {h.get('NAXIS2', '?')} px")
        print(f"  Observation date: {h.get('DATE-OBS', 'unknown')}")

    # Step 2 extract sources
    print(f"\n[2/6] Extracting sources...")
    red_sources = extract_sources(red_fits, os.path.join(outdir, "work_red"), "red")

    blue_sources = None
    if blue_fits:
        blue_sources = extract_sources(blue_fits, os.path.join(outdir, "work_blue"), "blue")

    # Step 3 quality filters
    print(f"\n[3/6] Applying quality filters...")
    red_filtered = apply_quality_filters(red_sources, "red")

    blue_filtered = None
    if blue_sources is not None:
        blue_filtered = apply_quality_filters(blue_sources, "blue")

    # Step 4 red-vs-blue comparison
    if blue_filtered is not None and len(blue_filtered) > 0:
        print(f"\n[4/6] Red-vs-blue transient detection...")
        candidates = find_transients(red_filtered, blue_filtered)
    else:
        print(f"\n[4/6] Skipping red-vs-blue (no blue plate)")
        candidates = red_filtered

    # Step 5 crossmatch against modern catalogs
    print(f"\n[5/6] Crossmatching {len(candidates)} sources against modern catalogs...")
    candidates = crossmatch_modern_catalogs(candidates, cache_dir,
                                              use_ir=args.ir)
    final = candidates[~candidates["is_known"]].copy()

    # Southern hemisphere exclusion
    final = apply_dec_cutoff(final, args.dec_min)

    # Cross-plate deduplication (multi-plate runs only)
    if args.dedup_dir:
        final = deduplicate_across_plates(final, args.dedup_dir)

    # Step 6 compare to VASCO if catalog provided
    if args.vasco:
        compare_to_vasco(final, args.vasco, plate)
    else:
        print(f"\n  No VASCO catalog provided, skipping comparison.")
        print(f"  Final candidates: {len(final)}")

    # Save
    out_csv = os.path.join(outdir, f"transients_{plate}.csv")
    final.to_csv(out_csv, index=False)
    print(f"\nResults saved to {out_csv} ({len(final)} transient candidates)")


if __name__ == "__main__":
    main()
