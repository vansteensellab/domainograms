# Domainograms

> **Perturbations in 3D genome organization can promote acquired drug resistance**
> 
> _Anna G Manjon; Stefano Giustino Manzo; Stefan Prekovic; Leon Potgeter; Tom van Schaik; Ning Qing Liu; Koen Flach; Daniel Peric-Hupkes; Stacey Joosten; Hans Teunissen; Anoek Friskes; Mila Ilic; Dorine Hintzen; Vinícius H Franceschini-Santos; Wilbert Zwart; Elzo de Wit; Bas van Steensel; René H Medema_

---

This repository contains a wrapper script to plot domainograms from GATC fragments files from pA-DamID experiments.
The code was used for the "**­­Perturbations in 3D genome organization can promote acquired drug resistance**" manuscript. 

With this Rscript file, you only need to provide the path to your `*-gatc.counts.txt.gz` files for Dam-only and antibody files.
Besides the domainogram, this script also produces in the output folder the correlation between the experimental and control data before and after LOESS correction.

# Installation

To run this script, first clone this repo. 
Then, if you need, you can create a conda environment with the provided `domainograms_conda_env.yml` file to install all its dependencies automatically:

```sh
conda env create -n plot_domainograms -f domainograms_conda_env.yml
```

# Usage

You can run this script from the terminal with

```sh
Rscript ./plot_domainograms.R [params]
```

