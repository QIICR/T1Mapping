# T1 Mapping

## Introduction
T1 Mapping is an extension of [3D Slicer software](http://www.slicer.org) that estimates effective tissue parameter maps (T1) from multi-spectral FLASH MRI scans with different flip angles. T1 mapping can be used to optimize parameters for a sequence, monitor diseased tissue, measure Ktrans in DCE-MRI and etc. 

See documentation at
[T1 Mapping Wiki](http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/T1Mapping)

## Functionality
* Take multi-spectral FLASH images with an arbitrary number of flip angles as input, and estimate the T1 values of the data for each voxel.
* Read repetition time(TR), echo time(TE) and flip angles from the Dicom header directly.
* Prostate, brain, head & neck, cervix, breast and etc. 

## Acknowledgments
Acknowledgments: This work is part of the National Alliance for Medical Image Computing (NA-MIC), funded by the National Institutes of Health through the NIH Roadmap for Medical Research, and by National Cancer Institute as part of the Quantitative Imaging Network initiative (U01CA154601) and QIICR (U24CA180918).

## References
1. [QIBA Phantom](https://sites.duke.edu/dblab/qibacontent/)
2. [Equations for Variable Flip Angle T1 Mapping](http://europepmc.org/articles/pmc3620726)
