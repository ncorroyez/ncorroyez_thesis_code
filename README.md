# Repository of Nathan Corroyez's PhD Thesis Code at INRAE / Maison de la Télédétection, 10/2023 to at least 09/2026 (planned)

## Description

Repository of code developed and written articles during Nathan Corroyez's PhD Thesis: "Remotely-sensed vegetation properties to improve microclimate models under forest canopy".

## Table of Contents

- [Setup](#Setup)
- [Preprocessing](#Preprocessing)
- [Fonctionnalities](#Fonctionnalities)
- [License](#License)
- [Contact](#Contact)
- [Acknowledgments](#Acknowledgments)
- [References](#References)

## Setup

At the place you desire, execute:
```bash
# Clone repo
git clone git@github.com:ncorroyez/ncorroyez_thesis_code.git

# Access to the repo
cd ncorroyez_thesis_code

# Setup directories architecture (01_DATA, 03_RESULTS) & install dependencies 
Rscript setup.R
```

#### Download LiDAR data and shapefiles

LiDAR acquisitions of Aigoual, Blois, and Mormal study sites, along with their related shapefiles are available. 
If you are from Maison de la Télédétection, directly copy the data from the _mo-pulse_ server (_PROJETS/2023_2026_These_Nathan_Corroyez/LiDAR). Otherwise, please contact me.

LiDAR data are provided in Lambert-93 format (1-las_l93), UTM 31N format (2-las_utm), and normalized UTM 31N format (Z -> H, 3-las_normalized_utm). The L93 to UTM 31N conversion can be done using the '0_convert_l93_into_utm.R' file and the UTM 31N to normalized UTM 31N with the '0_normalized_height.R' file. After the conversions, Lambert-93 data are not further used. 

Shapefiles are declined in three files for each site, including the site delimitation in both Lambert-93 and UTM 31N coordinates, and the tree species (BDForêt V2, in Lambert-93).

Once data are downloaded, put them in their associated directories (e.g. for Mormal data, put LiDAR data in 01_DATA/Mormal/LiDAR and shape files in 01_DATA/Mormal/Shape).

## Preprocessing

### Canopy Height Models

Canopy Height Models (CHMs) are deduced from the subtraction of the Digital Terrain Model (DTM) from the Digital Surface Model (DSM). The DTM is calculated via the _TIN_ algorithm and the DSM is obtained using a _pitfree_ algorithm. In these calculations, a LAS catalog of the LiDAR data is created thanks to the _lidR_ package (Roussel et al., 2024). 

```bash
Rscript LiDAR/1_create_chm.R
```

### LiDAR LAI and Other Metrics

LiDAR metrics are calculated via the _lidR_ package. LiDAR LAI is deduced thanks to a gap fraction method.

LiDAR LAI, along with CHM-related and point clouds-related metrics, are obtained via _1_calculate_lidar_metrics.R_ script:
```bash
Rscript LiDAR/2_calculate_lidar_metrics.R
```

### Sentinel-2 LAI

Sentinel-2 LAI is obtained using the PROSAIL model. The R package _prosail_ (Féret et al., 2024) is employed.

```bash
Rscript Sentinel_2/2_train_predict_prosail.R
```

### Masks

Several masks are created. The final masks that are used in further analysis are:
- Mask 1: masks low vegetation (<2m) areas, routes, and clearings. All tree species and lands are included.
- Mask 2: masks low vegetation (<2m) areas, routes, and clearings. Deciduous tree species and deciduous-coniferous mix (mainly deciduous) species are included.
- Mask 3: masks low vegetation (<2m) areas, routes, and clearings. Only deciduous tree species are included.

```bash
Rscript LiDAR/3_create_masks.R
```

## Fonctionnalities

### Heterogeneity and Depth Analysis

In this step, preliminary hypotheses about heterogeneity and depth are verified.

```bash
Rscript LiDAR/4_heterogeneity_analysis.R
```

### Correct Sentinel-2 LAI with LiDAR Information: The Machine Learning Way

The feature selection is done by a Sequential Features Selector.
Several models are tested: Random Forest, and Partial Least Square Regression.

LAI correction can be done either by training in full areas, in mixed deciduous-coniferous areas, or in deciduous-only areas.

```bash
Rscript LiDAR/5_s2_lai_correction_via_ml.R
```

#### On full areas

#### On deciduous-flexible areas

#### On deciduous-only areas

## License

This work is intended to be in the public domain (GPL-3).

## Contact

Please don't hesitate to initiate contact with me or one of my supervisors for any questions, remarks, or advice about this work.

- nathan.corroyez@inrae.fr (1st year PhD Student)
- nathan.corroyez14@gmail.com (non-academic email address)
- sylvie.durrieu@inrae.fr (director of research)
- jean-baptiste.feret@inrae.fr (co-director of research)
- jerome.ogee@inrae.fr (co-director of research)

## Acknowledgments

We would like to thank the CNES agency for the TOSCA Grant N°00007689 along with INRAE's MathNum Department and Agence Nationale de la Recherche via the MaCCMic ANR Project, (ANR-21-CE32-0012) for the cofounding of the PhD Thesis of Nathan Corroyez, and the IMPRINT ANR Project (ANR-19-CE32-0005) that funded the LiDAR acquisitions of Mormal, Blois, and Aigoual sites.

## References

Feret J, de Boissieu F (2024)._prosail: PROSAIL leaf and canopy radiative transfer model and inversion routines_. R package version 2.4.1, [https://gitlab.com/jbferet/prosail](https://gitlab.com/jbferet/prosail)

Roussel, J.R., Auty, D., Coops, N. C., Tompalski, P., Goodbody, T. R. H., Sánchez Meador, A., Bourdon, J.F., De Boissieu, F., Achim, A. (2021). _lidR: An R package for analysis of Airborne Laser Scanning (ALS) data_. Remote Sensing of Environment, 251 (August), 112061. [doi:10.1016/j.rse.2020.112061](https://doi.org/10.1016/j.rse.2020.112061).
