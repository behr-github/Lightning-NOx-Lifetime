# Lightning-NOx-Lifetime
Relevant analysis code used in Nault et al. GRL 2017

## Overview of important functions
* The gridding and application of averaging kernels was done using `load_and_grid_OMNO2.m`, `load_and_grid_DOMINO.m`, and `gc_column_omi_ak.m`
  * The `load_and_grid_*.m` functions handle gridding the NASA SP or KNMI DOMINO data to the GEOS-Chem grid cells, and weighting each contribution appropriately
  * The `gc_column_omi_ak.m` file handles both binning the averaging kernels to GEOS-Chem grid cells and applying them; which mode is executed is controlled by the global variables `run_mode`
* `integrate_geoschem_profile.m` was used to obtain GEOS-Chem tropospheric NO<sub>2</sub> columns *without* averaging kernels
* The functions contained under `BPCH_Functions` are those created by Sebastian D. Eastham and published [on the GEOS-Chem wiki](http://wiki.seas.harvard.edu/geos-chem/index.php/Matlab_software_tools_for_use_with_GEOS-Chem#Download)
* We used `read_geos_output.m` to import the regular `ctm.bpch` files into Matlab, and a combination of `read_gc_nd51.py` and `convert_py_nd51_struct.m` to read the ND51 diagnostic outputs
* Users looking to read in the netCDF files provided with our paper can read them into Matlab using `nc2gcstruct.m`
  * If you use another language, use `nc2gcstruct.m` to understand how variables/attributes in the netCDF file map to fields in the Matlab structures used elsewhere in the code
* The code contained in `DC3 Utils` was used to handle matching GEOS-Chem results to DC3 aircraft observations

## Contact
If you have questions, either open an issue on this repo or contact the corresponding author.
