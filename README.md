# Met-ocean data extraction from WINDSURFER/NORA3

[![View Met-ocean data extraction from WINDSURFER/NORA3 on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/103055-met-ocean-data-extraction-from-windsurfer-nora3)

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
