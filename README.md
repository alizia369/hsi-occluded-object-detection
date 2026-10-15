# Occluded Object Detection in Hyperspectral Images

**Ali Zia**
*PhD thesis, Griffith University, 2020 — arXiv version forthcoming*

[![Thesis](https://img.shields.io/badge/Thesis-Griffith%202020-blue)](https://research-repository.griffith.edu.au/items/6c4ccabf-f952-4780-8592-ac6f8ea62dae) ![arXiv](https://img.shields.io/badge/arXiv-forthcoming-lightgrey) ![Language](https://img.shields.io/badge/Language-MATLAB-orange) [![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

> **Status.** This is the original research code behind the publication,
> prepared for release. Duplicate, backup and autosave files, data and personal
> file paths were removed and third-party code was separated into
> `third_party/`; the algorithms themselves were not modified. It is research
> code rather than a polished toolbox, and it was developed on Windows MATLAB.

## Overview

Detecting an object that is partly hidden behind another is difficult from a
single RGB image, because occluder and occluded object often have similar
colour and the boundary between them carries little depth information. A
hyperspectral image provides two extra sources of evidence: the **materials**
in the scene, which can be separated by spectral unmixing, and **relative
depth**, which can be estimated from chromatic aberration and defocus blur
across bands.

This repository implements the occluded-object detection method from the
thesis *Exploring Extrinsic Spectral Properties of Hyperspectral Image for 3D
Computer Vision*. It fuses material abundance maps with depth cues on a graph
and localises the objects that are partially hidden. A standalone version of
this chapter has been submitted to arXiv; the identifier will be added here
once it is announced.

## Method

1. **Material distribution.** The cube is unmixed blindly with hierarchical
   clustering based on rank-two non-negative matrix factorisation
   (`hierclust2nmf`), and abundances are recovered with non-negative least
   squares (`nnlsm_blockpivot`). This gives one abundance map per material
   (`MainScript.m`, `MainScript2_snapScan.m`).
2. **Depth cues.** A sparse relative depth map is estimated at edges from
   band-wise defocus and chromatic aberration (`sparseDepthKLDivergence.m`,
   `getpatchesNew.m`), following the companion repository
   [`hsi-chromatic-defocus-depth`](https://github.com/alizia369/hsi-chromatic-defocus-depth).
3. **Graph fusion.** A graph Laplacian over the image encodes spectral
   similarity between neighbouring pixels (`graphLapL2.m`). Each material map is
   propagated under the depth constraints (`graphSolveLap.m`), so regions that
   share a material but lie at different depths are separated.
4. **Localisation.** The fused maps are thresholded, and the dominant connected
   regions are reported as bounding boxes of occluded objects
   (`MainScript3_CombiningCues.m`).

## Repository structure

```
setup_paths.m                    add all code folders to the MATLAB path
src/
  MainScript.m                   load a cube and unmix it (H2NMF)
  MainScript2_snapScan.m         the same for SnapScan cubes
  MainScript3_CombiningCues.m    fuse materials and depth, draw detections
  graphLapL2.m, graphSolveLap.m  graph construction and propagation
  sparseDepthKLDivergence.m      sparse depth cue
  MainScriptGraphLap*.m          fusion with alternative similarity measures
  TestScript*.m                  experiment scripts kept for reference
utils/                           cube loading and normalisation helpers
third_party/                     H2NMF, Canny, DBSCAN, Zhuo & Sim helpers, ENVI header reader
data/                            place your data here (not tracked)
```

## Requirements

- MATLAB (developed with R2019–R2020 on Windows)
- Image Processing Toolbox (`imbinarize`, `bwareafilt`, `regionprops`, …)
- Statistics and Machine Learning Toolbox (`pdist2`, `mvnrnd`)

## Getting started

The main files are MATLAB **scripts** that pass variables through the
workspace, so run them in order in one session:

```matlab
>> cd hsi-occluded-object-detection
>> setup_paths
>> MainScript2_snapScan        % or MainScript for a workspace.mat from the depth code
>> MainScript3_CombiningCues   % fusion and detection
```

Set the input folder at the top of the first script. `MainScript2_snapScan.m`
reads SnapScan cubes saved as `.mat` files (`snapScanSave2Matfiles.m`
converts the raw ENVI files). `MainScript.m` reads a `workspace.mat` written by
`hsi-chromatic-defocus-depth/src/sparse/MainScript.m`.

## Data

The lab scenes used in the thesis, with objects partly hidden behind others,
are not distributed with this repository; see `data/README.md` for the
expected inputs.

## Notes on reproducibility

- The number of materials `r` and the band range passed to the unmixing step
  were chosen per scene in the thesis.
- `MainScriptGraphLap{L2,SAM,SID,Proposed,Tang,rgb,SeperateEndmembers}.m`
  reproduce the comparison between similarity measures and inputs.

## Citation

If you use this code, please cite:

```bibtex
@misc{zia2026occluded,
  title         = {Occluded Object Detection in Hyperspectral Images},
  author        = {Zia, Ali},
  year          = {2026},
  archivePrefix = {arXiv},
  note          = {arXiv preprint, forthcoming}
}

@phdthesis{zia2020thesis,
  title  = {Exploring Extrinsic Spectral Properties of Hyperspectral Image for 3D Computer Vision},
  author = {Zia, Ali},
  school = {Griffith University},
  year   = {2020},
  url    = {https://research-repository.griffith.edu.au/items/6c4ccabf-f952-4780-8592-ac6f8ea62dae}
}
```

A `CITATION.cff` file is included, so GitHub's **Cite this repository** button gives the same reference.

## Licence and third-party code

The author's code is released under the **MIT licence** (see `LICENSE`).
Third-party code keeps its own terms, listed in `third_party/README.md`. In
particular, `third_party/zhuo_sim_defocus_2011/` is provided by its authors
**for academic use only**.

## Acknowledgements

This work was carried out as part of the author's PhD at the School of Information and Communication Technology, Griffith University. Thanks also to the authors of the third-party code listed above.

## Contact

Questions and bug reports are welcome through [GitHub issues](https://github.com/alizia369/hsi-occluded-object-detection/issues).
Ali Zia, La Trobe University — [www.ali-zia.me](https://www.ali-zia.me)

## Related repositories

This repository is one of a set released together from the same line of work on hyperspectral computer vision:

- [`hsi-3d-reconstruction`](https://github.com/alizia369/hsi-3d-reconstruction): 3D reconstruction from hyperspectral images (WACV 2015)
- [`hsi-relative-depth-estimation`](https://github.com/alizia369/hsi-relative-depth-estimation): Relative depth from a hyperspectral cube (DICTA 2015)
- [`hsi-chromatic-defocus-depth`](https://github.com/alizia369/hsi-chromatic-defocus-depth): Chromatic aberration + defocus for relative depth (IEEE TIP 2021)
- [`hsi-material-boundary-detection`](https://github.com/alizia369/hsi-material-boundary-detection): Material-based boundary detection (DICTA 2017)
- [`sugarcane-smut-hsi`](https://github.com/alizia369/sugarcane-smut-hsi): Early sugarcane smut detection (IVCNZ 2021)
- [`hsi-matlab-utils`](https://github.com/alizia369/hsi-matlab-utils): Shared hyperspectral I/O and normalisation helpers
