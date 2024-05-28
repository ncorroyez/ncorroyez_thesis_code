# Repository of Nathan Corroyez's PhD Thesis Code at INRAE / Maison de la Télédétection, 10/2023 to at least 09/2026 (planned)

## Description

Repository of code developed during Nathan Corroyez's PhD Thesis: "Remotely-sensed vegetation properties to improve microclimate models under forest canopy".

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

#### Download LiDAR data and shape files

LiDAR acquisitions of Aigoual, Blois, and Mormal study sites are available. 
If you are from Maison de la Télédétection, directly copy the data from mo-pulse server (_PROJETS/2023_2026_These_Nathan_Corroyez/LiDAR). Otherwise, please contact me.

LiDAR data are provided in Lambert 93 format (1-las_l93), UTM 31N format (2-las_utm), and normalized UTM 31N format (Z -> H, 3-las_normalized_utm). The L93 to UTM 31N conversion can be done using 'convert_l93_into_utm.R' file and the UTM 31N to normalized UTM 31N with 'normalized_height.R' file. After the conversions, L93 data are not further used. 

Once data are downloaded, put them in their associated directories (e.g. for Mormal data, put LiDAR data in 01_DATA/Mormal/LiDAR and shape files in 01_DATA/Mormal/Shape).

### Preprocessing

#### Canopy Height Models

```bash
Rscript LiDAR/0.create_mne.R
```

#### LiDAR LAI (at the moment, PAI)

a

#### Sentinel-2 LAI

a

#### Masks

a

### Fonctionnalities

#### Heterogeneity and Depth Analysis

a

#### Correct Sentinel-2 LAI with LiDAR Information: The Machine Learning Way

##### On full areas

##### On deciduous-flexible areas

##### On deciduous-only areas

### License

This work is intended to be in the public domain (GPL-3).

### Contact

Please don't hesitate to initiate contact with me or one of my supervisors for any questions, remarks or advice about this work.

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
