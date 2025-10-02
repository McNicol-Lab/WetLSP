# WetLSP v2

**Fork:** Updated algorithms (July 2025) to derived a land surface phenology product from PlanetScope imagery for FLUXNET-CH4 v2 sites. Data developed and initially shared for wetland methane research applications in the McNicol ecophilab, the ESIIL AI for Natural Methane Working Group, and the USGS Powell Synthesis for Wetlands.

To run updated code:  

- Clone the repo into your workspace
- Create directories listed at top of PLSP_Parameters_gm_v1.json and updated filepaths 
- Use `00_img_download.py` to retrieve data via PlanetScope API into created `rawImage` dir
- Run series of 5 processing scripts

Quirks:
* Note that original geojson files must be stored within `/AMFLX` subdirectory 
* `_gm_v1` workflow skips water mask steps in `03_LSP_script`, should not break process, but use if available


**Origin:** Algorithms to derive a land surface phenology product from PlanetScope imagery for AmeriFlux and NEON sites

- Publication: Moon, M., Richardson, A.D., Milliman, T. and Friedl, M.A., 2022. A high spatial resolution land surface phenology dataset for AmeriFlux and NEON sites. Scientific Data, 9(1), p.448. https://www.nature.com/articles/s41597-022-01570-5
- Data: Moon, M., A.D. Richardson, T. Milliman, and M.A. Friedl. 2023. Land Surface Phenology, Eddy Covariance Tower Sites, North America, 2017-2021. ORNL DAAC, Oak Ridge, Tennessee, USA. https://doi.org/10.3334/ORNLDAAC/2033

