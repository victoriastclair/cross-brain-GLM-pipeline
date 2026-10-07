# cross-brain-GLM-pipeline
Scripts for main analysis presented in St Clair et al. (under review)

Code repository for St Clair et al., (under review). Modelling children's brain activity in context: A cross-brain GLM approach to parent-child fNIRS hyperscanning. Oxford Open Neuroscience.

Please contact Dr Victoria St Clair with any questions (v.stclair@bbk.ac.uk).

**Requirements**

_Software_

  - MATLAB 2020B or later (installation: https://uk.mathworks.com/help/install/install-products.html)
  - Statistics and Machine Learning Toolbox (used by robustfit in the short-separation regression)

_Toolboxes/Functions_

  - Homer2 Toolbox (https://homer-fnirs.org/download/), used for preprocessing the raw fNIRS data into the _ppr.mat files read by this pipeline. Huppert, T. J., Diamond, S. G., Franceschini, M. A., & Boas, D. A. (2009). HomER: A review of time-series analysis methods for near-infrared spectroscopy of the brain. Applied Optics, 48(10), D280–D298. https://doi.org/10.1364/ao.48.00d280
  - SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/), used for the canonical haemodynamic response function (spm_hrf). Penny, W. D., Friston, K. J., Ashburner, J. T., Kiebel, S. J., & Nichols, T. E. (Eds.). (2011). Statistical parametric mapping: The analysis of functional brain images. Elsevier.
  - Short-separation regression functions: Abdalmalak, A. et al. (2022). Effects of systemic physiology on mapping resting-state networks using functional near-infrared spectroscopy. Frontiers in Neuroscience, 16, 803297. https://doi.org/10.3389/fnins.2022.803297. Questions about these functions ('AdjustTemporalShift_fnirs_course.m' and 'PhysiologyRegression_GLM_fnirs_course.m') should be directed to the original authors: Professor Rickson Mesquita (r.c.mesquita@bham.ac.uk) and Dr Sergio Luiz Novi Jr (novisl@ifi.unicamp.br).

**Example data**

After cloning the repository, set the paths in config.m and add SPM12 to the MATLAB path. The pipeline expects one folder per dyad (Dyad*/C*/C*_ppr.mat for the child and Dyad*/M*/M*_ppr.mat for the mother), plus three behavioural files per dyad in the analysis path: the eye-tracking world-view timestamps (*_wv.csv), the eye-tracking condition annotations (*_ann.csv), and the ELAN coding of parent scaffolding (*.txt). Channels excluded during manual quality checks are read from ChannelExclusions.xlsx. Templates for these file structures are provided in /example_data/.

**Structure**

***xGLM Pipeline***
1. config.m = study-specific parameters to be set before running pipeline (paths, ROI channels, short-separation channels, behavioural codes, sampling rate)
2. GLM_pipeline.m = main script
  - option to run with or without short-separation channel regression (line 27);
    - excludes noisy channels and applies short-separation regression (prepareData.m, excludeChannels.m);
    - averages channels within each ROI for HbO and HbR (averageROIs.m);
    - builds condition regressors (collaboration, individual) from the fNIRS triggers;
    - aligns eye-tracking and fNIRS time, and builds the scaffolding regressor from the ELAN coding (collaboration only);
    - fits three GLMs per child ROI and chromophore: behaviour only (conditions + scaffolding), brain only (conditions + mother's ROI signals per condition), and behaviour + brain;
    - stores betas, R2 and adjusted R2 for each model in the tmp struct, and scaffolding descriptives in fullSummary and dyadInfo
