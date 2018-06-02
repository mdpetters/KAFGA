## KAppa Functional Group Analysis (KAFGA)
<img src="docs/Logo.png" width="100">

KAFGA is MATLAB/GNU Octave code to predict the hygroscopicity parameter kappa for organic compounds from functional group composition.

__Short summary__ Organic particles suspended in air serve as nucleation seeds for droplets in atmospheric clouds. Over time their chemical composition changes towards more functionalized compounds. This work presents a model that can predict an organic compounds' ability promote the nucleation of cloud drops from its functional group composition. Hydroxyl, carboxyl, aldehyde, hydroperoxide, carbonyl, and ether moieties promote droplet nucleation. Methylene and nitrate moieties inhibit droplet nucleation.

## Documentation
The algorithms are described in the [manuscript](docs/gmd-9-111-2016.pdf) and associated [supplement](docs/gmd-2015-172-si.pdf). Model calculations graphed in Figures 1-3 (main text) and Figure 1 (supplement) are in the example folders. The most straightforward example is fs01.<br>

_Model input files_
 - each compound is represented by a model input file in folder fs01comp
 - an input file consists of 4 columns and 9 rows
 - rows and columns correspond to main and subgroups, respectively
 - Specifically:

                           Sub1        Sub2        Sub3        Sub4
               Alkane      CH3         CH2         CH          C
               Alcohol     OH          ---         ---         ---
               Water       H2O         ---         ---         ---
               Carbonyl    CH3C(=O)    CH2(C=O)    ---         ---
               Aldehyde    H(C=O)      ---         ---         ---
               Ether       CH3(O)      CH2(O)      CH(O)       THF
               Acid        C(=O)OH     HC(=O)OH    ---         ---
               Nitrate     CH2(ONO2)   CH(ONO2)    C(ONO2)     ---
               Peroxide    CH2(OOH)    CH(OOH)     C(OOH)      ---

 - Example input file for ethanol CH3-CH2OH

                             1           1           0            0
                             1           0	        0            0
                             0           0           0            0
                             0           0           0            0
                             0           0           0            0
                             0           0           0            0
                             0           0           0            0
                             0           0           0            0
                             0           0           0            0
 - Input files for all compounds in gmd-2015-072 are provided in examples

 - A component input structure is needed. These are generated in the
   script, for example, fs01_compounds.m. The fields are: 'name':
   compound name, 'file': name of input text file (see above), 'sol'
   observed solubility', MS_obs: observed molecular weight, rho_obs:
   observed density, Cx, Hx, Ox, Nx: # of C, H, N, O atoms,  D:
   diameter for which CCN activity is computed. Note that observed
   quantities are not needed for the model. To see an example for
   initialization/generation of model compounds without observations
   see the structure returned for f03_compounds().

_Model execution and output_
 - For example model execution see file fs01.m
 - Multiple files can be present in the input directory
 - The output directory must exist. It is possible to combine them
   into a single directory
 - A file with interaction coefficients anm, located in src/txt can be
   specified at runtime. Current options are anm_standard and
   anm_revised. anm standard corresponds to standard UNIFAC,
   anm_revised to the coefficients used in gmd-2015-172
 - To execute now compounds, create a copy of for example f03 and
   rename it. If the copy is not in the example folder, the sourcepath
   needs to be changed in the execution file, as well as the src files
   main.m, load_coefficients.m and unifac.m
 - A summary output file is written. The file is appended each time
   the model runs. The file lists: filenamame (compound), # of C, # of H, # of O, # of N, molecular weight, predicted molar volume, Flory-Huggins kappa, Raoult kappa, UNIFAC kappa (reported in the paper), water
   activity of saturated solution, molefraction of xw of phase bounary
   1,molefraction of xw of phase bounary 2.
 - The complete output is also saved as a .mat file for quick loading
 - An ASCII file (.out) is written that contains xw, xs, gamma_w, gamma_s
 - Example f02 contains an initialization script for serial and
   parallel execution. Parallel execution is not supported in
   Octave. Serial execution for f02 produces the same result in octave
   and matlab. Postprocessing in fs01 is not supported in
   Octave.

## Citation
This work was supported by the United States Department of Energy, Office of Science, Biological and Environment Research, Grant number DE-SC0010470. If you use _KAFGA_ in your research, please cite

Petters, M. D., Kreidenweis, S. M., and Ziemann, P. J. (2016). <i> Prediction of cloud condensation nuclei activity for organic compounds using functional group contribution methods, </i> Geoscientific Model Development, 9, 111-124, https://doi.org/10.5194/gmd-9-111-2016.
