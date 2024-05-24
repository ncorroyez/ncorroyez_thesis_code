# Repository of Nathan Corroyez's PhD Thesis Code at INRAE / Maison de la Télédétection, 10/2023 to 10/2026 (forecast)

## Description

Repository of code developed during Nathan Corroyez's PhD Thesis: "Remotely-sensed vegetation properties to improve microclimate models under forest canopy".

## Table of Contents

- [Setup](##setup)
- [Preprocessing](##preprocessing)
- [Fonctionnalities](##fonctionnalities)
- [License](##license)
- [Contact](##contact)
- [Acknowledgments](##acknowledgments)

## Setup

At the place you desire, execute:
```bash
# Clone repo
git clone git@github.com:ncorroyez/ncorroyez_thesis_code.git

# Access to the repo
cd ncorroyez_thesis_code

# Setup directories architecture & install dependencies 
Rscript setup.R
```

#### Download LiDAR data and shape files

LiDAR acquisitions of Aigoual, Blois, and Mormal study sites are available. 
If you are from Maison de la Télédétection, directly copy the data from mo-pulse server (_PROJETS/2023_2026_These_Nathan_Corroyez/LiDAR). Otherwise, please contact me.

Once data are downloaded, put them in their associated directories (e.g. for Mormal data, put LiDAR data in 01_DATA/Mormal/LiDAR and shape files in 01_DATA/Mormal/Shape).

### Preprocessing

#### Canopy Height Models

a

#### LiDAR LAI (at the moment, PAI)

a

#### Sentinel-2 LAI

a

#### Masks

a

### Fonctionnalities

####

### License

This work is intended to be in the public domain.

### Contact

Please don't hesitate to initiate contact with me or one of my supervisors for any questions, remarks or advice about this work.

- nathan.corroyez@inrae.fr (1st year PhD Student)
- nathan.corroyez14@gmail.com (non-academic email address)
- sylvie.durrieu@inrae.fr (director of research)
- jean-baptiste.feret@inrae.fr (co-director of research)
- jerome.ogee@inrae.fr (co-director of research)

## Acknowledgments
We would like to thank the CNES agency for the TOSCA Grant N°00007689 along with INRAE's MathNum Department and Agence Nationale de la Recherche via the MaCCMic ANR Project, (ANR-21-CE32-0012) for the cofounding of the PhD Thesis of Nathan Corroyez, and the IMPRINT ANR Project (ANR-19-CE32-0005) that funded the LiDAR acquisitions of Mormal, Blois, and Aigoual sites.