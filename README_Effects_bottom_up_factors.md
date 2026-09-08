# Effects of bottom-up factors on growth and toxin content of *Alexandrium pseudogonyaulax*

R code and working data accompanying Möller et al. (2024), *Limnology and
Oceanography*. Archived, not maintained.

Two laboratory experiments on three strains of the harmful dinoflagellate
*Alexandrium pseudogonyaulax* isolated from the Danish Limfjord. One tests
nitrogen source, the other light intensity, both measuring growth, goniodomin
content, pigments and photophysiology.

## Citation

> Möller, K., Thoms, S., Tillmann, U., Krock, B., Koch, F., Peeken, I.,
> Meunier, C.L. (2024). Effects of bottom-up factors on growth and toxin content
> of a harmful algae bloom dinoflagellate. *Limnology and Oceanography* 69(6),
> 1335-1349. https://doi.org/10.1002/lno.12576

## Data

The curated, citable datasets are archived on PANGAEA:

| Dataset | DOI |
|---|---|
| Bundle | [10.1594/PANGAEA.965195](https://doi.org/10.1594/PANGAEA.965195) |
| Nitrogen source experiment | [10.1594/PANGAEA.965196](https://doi.org/10.1594/PANGAEA.965196) |
| Light intensity experiment | [10.1594/PANGAEA.965197](https://doi.org/10.1594/PANGAEA.965197) |

**On the two versions of the data.** The files on PANGAEA were restructured for
archiving: long format, one row per observation, explicit column names, added
metadata and taxonomic identifiers. The scripts in this repository were written
against the original working files, which are included here under
`Light_intensities/raw/` and `N-sources/R_txt_files/` so the analysis runs as it
did for the paper. **PANGAEA is the authoritative version and the one to cite.**
The files in this repository exist only to reproduce these scripts.

Strains: L2-D2, L4-B1 and L4-B9. The light experiment uses all three, the
nitrogen experiment uses L2-D2 and L4-B1.

## Contents

### `Light_intensities/`

| Script | Purpose |
|---|---|
| `Growth_curve_light.R` | Cell density time series per strain and light level |
| `Growth_rate_light.R` | Specific growth rates from the exponential phase |
| `Goniodomins_light.R` | Goniodomin cell quotas, normalised to cell size and POC/PON |
| `Photo_irradiance_curves_light.R` | Photosynthesis-irradiance curve fitting from FRRF data |
| `Dark_parameters_photo_irradiance_curves_light.R` | Dark-adapted FRRF parameters |
| `Pigments_light.R` | Pigment composition, cell quotas and pigment ratios |

Working data in `raw/`.

### `N-sources/`

| Script | Purpose |
|---|---|
| `Growth_curves_N_sources.R` | Cell density time series per strain and nitrogen source |
| `Growth_rate_N_sources.R` | Specific growth rates and maximum densities |
| `Goniodomins_N_sources.R` | Goniodomin cell quotas |
| `Cell_sizes_N_sources.R` | Cell size distributions |
| `Particular_nutrients_N_sources.R` | Particulate organic carbon and nitrogen |

Working data in `R_txt_files/`.

Each script is standalone. Every script header names the paper and the PANGAEA
DOI.

## Running it

Open the repository as an RStudio project so that paths resolve from the
repository root, then run any script top to bottom.

```r
source("Light_intensities/R-files/Growth_rate_light.R")
```

Packages are loaded on demand via `pacman::p_load`.

## Known limitations

The scripts were archived as they were used and have not been refactored.

- Four referenced files are not included: `PI_all.txt`, `cell_counts_L2D2.txt`,
  `cell_counts_L4B1.txt`, `cell_counts_L4B9.txt`. Retrieve the equivalent data
  from the PANGAEA archive.
- `ETR_max.txt` is included but not read by any script. It holds maximum
  electron transport rates derived from the PI curves.
- Photosynthesis-irradiance curve fitting is done per strain in long, repeated
  blocks rather than a loop. `Photo_irradiance_curves_light.R` is roughly 3,000
  lines for this reason.
- Package loading is scattered through the scripts rather than collected at the
  top, so a partial run may hit a missing package mid-script.

## Status

Archived on publication. Kept for transparency and reproducibility. Issues are
not monitored.

## License

Code: MIT. Data: CC-BY 4.0, as archived on PANGAEA.
