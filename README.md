# Brain Connectivity Analysis in MATLAB
***This code is for research purposes only.***

This repository contains the code associated with the publication: *paper to be published*.

## Overview

### `conn.m`
Main script that computes connectivity features from patients EEG following the classical steps: pre-processing, nodes and edges computations, and features extraction.

You need to specify:

(1) all parameters concerning your connectivity analyses:
- *f_bands*: list of frequency bands of interest
- *des_fs*: desired sampling frequency
- *channels*: list of EEG channels of interest
- *des_ref*: desired reference
- *reref*: list of re-referencing schemes of interest
- *rej_t*: artifact rejection threshold (µV)
- *window*: window  size (in seconds)
- *overlap*: overlap between windows (in seconds)
- *measrs*: list of association measures of interest
- *featrs*: number of connectivity features to extract

(2) the specifics regarding your database:
- *files*: path to your patients files
- *save*: path to the folder in which you want to save your analyses
- *patients*: number of patients in your cohort
- *start*: analysis start time (in seconds)
- *duration*: length of analysis (in seconds)
- *diff_pat*: list of the patients you wish to re-reference to the *desired ref* and/or to resample to the *desired fs* specified above

The data of each patient should be stored in a MAT structure with 2 fields:
- *channels*: structure of size 1 x *number of channels* with 2 fields: Name and Fs
- *data*: double of size *length of recordings* [samples] x *number of channels* containing the EEG voltages (µV)

### `data.m`
Class handling all processing on the raw data, including data preparation, pre-processing, and window splitting.

### `measures.m`
Class handling all association measures implementations (definition and parametrization) and the computations of the associations between nodes to generate the adjacency matrix.

### `network.m`
Class handling the extraction of network features from the adjacency matrix.

## Dependencies
MATLAB TOOLBOXES:
- Signal Processing Toolbox
- Deep Learning Toolbox
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox

OTHER TOOLBOXES:
- Brain Connectivity Toolbox. This is a freely available and open source Matlab toolbox for complex brain-network analysis. You can download it [here](https://sites.google.com/site/bctnet/). For more details, please refer to: *Complex network measures of brain connectivity: Uses and interpretations. Rubinov M, Sporns O (2010) NeuroImage 52:1059-69.*

## Usage
1. Ensure you have MATLAB installed (recommended version: R2024b or later).
2. Clone or download this repository and add the project folder to your MATLAB path.
3. Open MATLAB and navigate to the project directory.
4. Update parameters and database information in the `conn.m` file.
5. Run `conn.m` to compute brain connectivity features for every patient.

## Contributors
Lise Cottin

## License
If you use this algorithm for a publication (in a journal, in a conference, etc.), please cite the related publications (see below). The license attached to this toolbox is GPL v2, see https://www.gnu.org/licenses/gpl-2.0.txt. From https://www.gnu.org/licenses/gpl-2.0.html, it implies: This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

## Citation
If you use this code or data in your research, please cite the following paper: *paper to be published*.

and reference the archived version of the code:

> **Cottin, L., Damien, C., Anzalone, L., Aeby, A., Gaspard, N., & Nonclercq, A. (2026).**  
> *Code for: EEG Functional Connectivity Methods for Seizure-related Disorders* [Computer software].  
> Zenodo.  
> [![DOI](https://zenodo.org/badge/1105454056.svg)](https://doi.org/10.5281/zenodo.17779210) 

A `CITATION.cff` file is included in this repository.  
Click **“Cite this repository”** on the right sidebar of the GitHub page to automatically generate a BibTeX or APA citation.

### BibTeX

```bibtex

@software{Cottin2025CodeEEGConnectivity,
  author    = {Lise Cottin and Charlotte Damien and Luca Anzalone and Alec Aeby and Nicolas Gaspard and Antoine Nonclercq},
  title     = {Code for: EEG Functional Connectivity Methods for Seizure-related Disorders},
  year      = {2026},
  doi       = {10.5281/zenodo.17779211},
  url       = {https://github.com/username/repo},
  version   = {0}
}
