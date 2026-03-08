# Machine Learning Framework for Functional Classification of GOF and LOF Variants in KCNQ2

This repository contains the code and analysis developed for predicting the functional effects of missense variants in KCNQ2, a voltage-gated potassium channel associated with early-onset epilepsy.
The project was conducted within the BEATKCNQ consortium of the European Partnership for Personalised Medicine (EP PerMed), and is based on an adaptation of the machine-learning framework developed by Heyne et al. (2020) for predicting GOF and LOF variants in voltage-gated ion channels.

---

# Repository Structure

## evaluations

This folder contains scripts and outputs used to evaluate the performance of the adapted model.

The analyses include:
- model performance metrics
- comparison of evaluation strategies
- result summaries generated from the prediction outputs.

---

## heyne_files

This directory documents the files required from the original Heyne et al. (2020) repository.

To avoid duplicating large files and to preserve the original source, these files are not stored directly in this repository. Instead, a README file in this directory explains which files should be downloaded from the original repository and how they are used in the analysis.

Original repository: https://github.com/heyhen/funNCion

---

## kcnq2_pipeline

This directory contains the scripts used to generate the KCNQ2-specific feature set and prediction pipeline.

The scripts include:
- feature calculation
- variant density calculations
- preparation of feature tables for model input
- execution of the adapted prediction workflow.

These scripts reproduce the feature generation and model input used in the KCNQ2 analysis.

---

# Reference and Acknowledgments

Heyne H. O., et al. (2020): Predicting Functional Effects of Missense Variants in Voltage-Gated Sodium and Calcium Channels.

We thank Henrike O. Heyne for sharing the code and feature calculation pipelines, and for helpful discussions that contributed to this work.

We thank Sarah Weckhuysen and Maurizio Taglialatela for collecting, characterizing, and providing the KCNQ2 variant dataset used in this study, as part of the BEATKCNQ project.

---

# Author

Yuval Oren  
Supervised by Prof. Nir Ben-Tal

Tel-Aviv University, February 2026
