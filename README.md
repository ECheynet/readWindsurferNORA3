# Met-ocean data extraction with NORA3

[![View Met-ocean data extraction from WINDSURFER/NORA3 on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/103055-met-ocean-data-extraction-from-windsurfer-nora3)
[![DOI](https://zenodo.org/badge/435617642.svg)](https://zenodo.org/badge/latestdoi/435617642)
[![Donation](https://camo.githubusercontent.com/a37ab2f2f19af23730565736fb8621eea275aad02f649c8f96959f78388edf45/68747470733a2f2f77617265686f7573652d63616d6f2e636d68312e707366686f737465642e6f72672f316339333962613132323739393662383762623033636630323963313438323165616239616439312f3638373437343730373333613266326636393664363732653733363836393635366336343733326536393666326636323631363436373635326634343666366536313734363532643432373537393235333233303664363532353332333036313235333233303633366636363636363536353264373936353663366336663737363737323635363536653265373337363637)](https://www.buymeacoffee.com/echeynet)

## Summary
Hourly 2D wave spectra and instantaneous wave field from WINDSURFER/NORA3 (3 km grid) hindcast are extracted using the OPeNDAP framework.

## Content
The repository contains:
  - The function windsurfer.m which gets the metocean conditions for a specific time and location
  - The function get2DSS.m, which gets the directional wave spectra 
  - The file world.mat, which contains the coordinates of coastlines. This file is used for visualization purpose only.
  -  The Matlab livescript Documentation_gridded_data.mlx, which illustrates the use of windsurfer.m.
  -  The Matlab livescript Documentation_wave_spectra.mlx, which illustrates the use of get2DSS.m

The repository is unlikely bug-free, so if you have any question or comment, please ask!
