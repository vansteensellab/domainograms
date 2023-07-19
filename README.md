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

# Reproducing manuscript domainograms

To get the domainograms for the manuscript, the `plot_domainograms_manjon_etal.R` script was called like:

```sh
DATA_DIR="/DATA/usr/s.manzo/Projects/SGM20220704_pA_DamID_seq_AnnaB5_NewABCB1_clones/results/counts/bin-gatc"
all_smp=( pADamID-TXR3 pADamID-TXR4 pADamID-TXR5 pADamID-TXR6 pADamID-iCA3P_TXR )
all_ctr=( pADamID-DPURO3 pADamID-DPURO3 pADamID-RPE0 pADamID-RPE0 pADamID-iCA3P )

for i in {0..4}; do

        SMP="${all_smp[${i}]}"
        CTR="${all_ctr[${i}]}"

        Rscript plot_domainograms_manjon_etal.R \
        --recombination FALSE \
        --exp_dam_files ${DATA_DIR}/${SMP}_r1_Dam-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_r2_Dam-gatc.counts.txt.gz \
        --exp_antibody_files ${DATA_DIR}/${SMP}_r1_Lmnb2-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_r2_Lmnb2-gatc.counts.txt.gz \
        --ctrl_dam_files ${DATA_DIR}/${CTR}_r1_Dam-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_r2_Dam-gatc.counts.txt.gz \
        --ctrl_antibody_files ${DATA_DIR}/${CTR}_r1_Lmnb2-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_r2_Lmnb2-gatc.counts.txt.gz \
        --annotation /DATA/usr/b.v.steensel/Projects/LADrewiring/data/ucsc/hg38/ncbiRefSeq.txt \
        --output_dir VF230419_${SMP/pADamID-/}_vs_${CTR/pADamID-/} \
        --expression_tbl FALSE \
        --organism human \
        --smooth_window_size 301 \
        --my_genes "ABCB1" \
        --recombination_centered FALSE \
        --plot_title ${SMP/pADamID-/}_vs_${CTR/pADamID-/}
done

################################################################################################################################################

SMP="pADamID-iCA3_c7_TXR"
CTR="pADamID-iCA3_c7"
        Rscript plot_domainograms_manjon_etal.R \
                --recombination FALSE \
                --exp_dam_files ${DATA_DIR}/${SMP}_r1_Dam-gatc.counts.txt.gz\
                --exp_antibody_files ${DATA_DIR}/${SMP}_r1_Lmnb2-gatc.counts.txt.gz\
                --ctrl_dam_files ${DATA_DIR}/${CTR}_r1_Dam-gatc.counts.txt.gz\
                --ctrl_antibody_files ${DATA_DIR}/${CTR}_r1_Lmnb2-gatc.counts.txt.gz\
                --annotation /DATA/usr/b.v.steensel/Projects/LADrewiring/data/ucsc/hg38/ncbiRefSeq.txt \
                --output_dir VF230419_${SMP/pADamID-/}_vs_${CTR/pADamID-/} \
                --expression_tbl FALSE \
                --organism human \
                --smooth_window_size 301 \
                --my_genes "ABCB1" \
                --recombination_centered FALSE \
                --plot_title ${SMP/pADamID-/}_vs_${CTR/pADamID-/}

################################################################################################################################################

DATA_DIR="/DATA/scratch/usr/t.v.schaik/proj/sManzo_pADamID/ts210903_samples_Anna_pADamID_RPE/results/counts/bin-gatc"
all_smp=( RPE_iCut_BIX_2 RPE_iCut_GSK126_500 RPE_iCut_5AZA_62_5)
all_ctr=( RPE_iCut_DMSO RPE_iCut_DMSO RPE_iCut_DMSO)

for i in {0..2}; do
        SMP="${all_smp[${i}]}"
        CTR="${all_ctr[${i}]}"

        Rscript plot_domainograms_manjon_etal.R \
                --recombination FALSE \
                --exp_dam_files ${DATA_DIR}/${SMP}_Dam_R1-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_Dam_R2-gatc.counts.txt.gz \
                --exp_antibody_files ${DATA_DIR}/${SMP}_LaminB2_R1-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_LaminB2_R2-gatc.counts.txt.gz \
                --ctrl_dam_files ${DATA_DIR}/${CTR}_Dam_R1-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_Dam_R2-gatc.counts.txt.gz \
                --ctrl_antibody_files ${DATA_DIR}/${CTR}_LaminB2_R1-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_LaminB2_R2-gatc.counts.txt.gz \
                --annotation /DATA/usr/b.v.steensel/Projects/LADrewiring/data/ucsc/hg38/ncbiRefSeq.txt \
                --output_dir VF230419_${SMP/RPE_iCut_/}_vs_${CTR/RPE_iCUT_/} \
                --expression_tbl FALSE \
                --organism human \
                --smooth_window_size 301 \
                --my_genes "ABCB1" \
                --recombination_centered FALSE \
                --plot_title ${SMP/RPE_iCut_/}_vs_${CTR/RPE_iCut_/}
done

################################################################################################################################################

DATA_DIR="/DATA/scratch/usr/t.v.schaik/proj/sManzo_pADamID/ts210903_samples_Anna_pADamID_RPE/results/counts/bin-gatc"
all_smp=( RPE_iCut_BIX_2 RPE_iCut_GSK126_500 RPE_iCut_5AZA_62_5)
all_ctr=( RPE_iCut_DMSO RPE_iCut_DMSO RPE_iCut_DMSO)
for i in {0..2}; do 
        SMP="${all_smp[${i}]}"
        CTR="${all_ctr[${i}]}"

        Rscript plot_domainograms_manjon_etal.R \
                --recombination FALSE \
                --exp_dam_files ${DATA_DIR}/${SMP}_Dam_R1-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_Dam_R2-gatc.counts.txt.gz \
                --exp_antibody_files ${DATA_DIR}/${SMP}_LaminB2_R1-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_LaminB2_R2-gatc.counts.txt.gz \
                --ctrl_dam_files ${DATA_DIR}/${CTR}_Dam_R1-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_Dam_R2-gatc.counts.txt.gz \
                --ctrl_antibody_files ${DATA_DIR}/${CTR}_LaminB2_R1-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_LaminB2_R2-gatc.counts.txt.gz \
                --annotation /DATA/usr/b.v.steensel/Projects/LADrewiring/data/ucsc/hg38/ncbiRefSeq.txt \
                --output_dir VF230419_${SMP/RPE_iCut_/}_vs_${CTR/RPE_iCUT_/} \
                --expression_tbl FALSE \
                --organism human \
                --smooth_window_size 301 \
                --my_genes "ABCB1" \
                --recombination_centered FALSE \
                --plot_title ${SMP/RPE_iCut_/}_vs_${CTR/RPE_iCut_/}
                done

################################################################################################################################################

DATA_DIR="/DATA/scratch/usr/t.v.schaik/proj/sManzo_pADamID/ts210903_samples_Anna_pADamID_RPE/results/counts/bin-gatc"
all_smp=( RPE_CRISPRa_ABCB1 RPE_CRISPRa_ABCB4 RPE_CRISPRa_combo RPE_CRISPRa_RUNCD3B)
all_ctr=( RPE_CRISPRa_WT RPE_CRISPRa_WT RPE_CRISPRa_WT RPE_CRISPRa_WT)

for i in {0..3}; do 
        SMP="${all_smp[${i}]}"
        CTR="${all_ctr[${i}]}"

        Rscript plot_domainograms_manjon_etal.R \
                --recombination FALSE \
                --exp_dam_files ${DATA_DIR}/${SMP}_Dam_R1-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_Dam_R2-gatc.counts.txt.gz \
                --exp_antibody_files ${DATA_DIR}/${SMP}_LaminB1_R1-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_LaminB1_R2-gatc.counts.txt.gz \
                --ctrl_dam_files ${DATA_DIR}/${CTR}_Dam_R1-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_Dam_R2-gatc.counts.txt.gz \
                --ctrl_antibody_files ${DATA_DIR}/${CTR}_LaminB1_R1-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_LaminB1_R2-gatc.counts.txt.gz \
                --annotation /DATA/usr/b.v.steensel/Projects/LADrewiring/data/ucsc/hg38/ncbiRefSeq.txt \
                --output_dir VF230419_${SMP/RPE_CRISPRa_/}_vs_${CTR/RPE_CRISPRa_/} \
                --expression_tbl FALSE \
                --organism human \
                --smooth_window_size 301 \
                --my_genes "ABCB1" \
                --recombination_centered FALSE \
                --plot_title ${SMP/RPE_CRISPRa_/}_vs_${CTR/RPE_CRISPRa_/}
done

################################################################################################################################################

DATA_DIR="/DATA/scratch/usr/t.v.schaik/proj/sManzo_pADamID/ts210903_samples_Anna_pADamID_RPE/results/counts/bin-gatc"
SMP="RPE_dPC3Mi_TxR_s20"
CTR="RPE_dPC3_WT"
Rscript plot_domainograms_manjon_etal.R \
                --recombination FALSE \
                --exp_dam_files ${DATA_DIR}/${SMP}_Dam_R1-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_Dam_R2-gatc.counts.txt.gz \
                --exp_antibody_files ${DATA_DIR}/${SMP}_LaminB1_R1-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_LaminB1_R2-gatc.counts.txt.gz \
                --ctrl_dam_files ${DATA_DIR}/${CTR}_Dam_R3-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_Dam_R4-gatc.counts.txt.gz \
                --ctrl_antibody_files ${DATA_DIR}/${CTR}_LaminB1_R3-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_LaminB1_R4-gatc.counts.txt.gz \
                --annotation /DATA/usr/b.v.steensel/Projects/LADrewiring/data/ucsc/hg38/ncbiRefSeq.txt \
                --output_dir VF230419_${SMP/RPE_/}_vs_${CTR/RPE_/} \
                --expression_tbl FALSE \
                --organism human \
                --smooth_window_size 301 \
                --my_genes "ABCB1" \
                --recombination_centered FALSE \
                --plot_title ${SMP/RPE_/}_vs_${CTR/RPE_/}

################################################################################################################################################

DATA_DIR="/DATA/scratch/usr/t.v.schaik/proj/sManzo_pADamID/ts210903_samples_Anna_pADamID_RPE/results/counts/bin-gatc"
all_smp=( RPE_iCut_LBRKO RPE_iCut_LMNAKO RPE_iCut_LMNB1KO)
all_ctr=( RPE_iCut_WT RPE_iCut_WT RPE_iCut_WT)

for i in {0..2}; do 
        SMP="${all_smp[${i}]}"
        CTR="${all_ctr[${i}]}"

        Rscript plot_domainograms_manjon_etal.R \
                --recombination FALSE \
                --exp_dam_files ${DATA_DIR}/${SMP}_Dam_R1-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_Dam_R2-gatc.counts.txt.gz \
                --exp_antibody_files ${DATA_DIR}/${SMP}_LaminB2_R1-gatc.counts.txt.gz,${DATA_DIR}/${SMP}_LaminB2_R2-gatc.counts.txt.gz \
                --ctrl_dam_files ${DATA_DIR}/${CTR}_Dam_R3-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_Dam_R4-gatc.counts.txt.gz \
                --ctrl_antibody_files ${DATA_DIR}/${CTR}_LaminB2_R3-gatc.counts.txt.gz,${DATA_DIR}/${CTR}_LaminB2_R4-gatc.counts.txt.gz \
                --annotation /DATA/usr/b.v.steensel/Projects/LADrewiring/data/ucsc/hg38/ncbiRefSeq.txt \
                --output_dir VF230419_${SMP/RPE_iCut_/}_vs_${CTR/RPE_iCut_/} \
                --expression_tbl FALSE \
                --organism human \
                --smooth_window_size 301 \
                --my_genes "ABCB1" \
                --recombination_centered FALSE \
                --plot_title ${SMP/RPE_iCut_/}_vs_${CTR/RPE_iCut_/}
done

```