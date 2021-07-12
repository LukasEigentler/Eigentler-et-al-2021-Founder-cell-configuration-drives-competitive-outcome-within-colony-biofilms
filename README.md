# Eigentler-et-al-2021-Founder-cell-configuration-drives-competitive-outcome-within-colony-biofilms 
This repository contains code associated with the paper Eigentler et al. (2021) BiorXiv preprint DOI: 10.1101/2021.07.08.451560.
All code is provided through GNU General Public License v3.0. Please see the license file for more information

The folder "Model_simulations" contains Matlab code to numerically solve the mathematical models presented in the paper and calculate the access to free space score associated with the model simulations. The main file is "model_simulation.m" and requires Matlab's PDEToolbox.
Author: Lukas Eigentler (leigentler001@dundee.ac.uk)

The folder "Image_analysis" contains a macro that analyses biofilm images contained in .lif files using Fiji/ImageJ.
Author: Graeme Ball (g.ball@dundee.ac.uk)

The folder "Access_to_free_space_score_experiments" contains a Matlab script to calculate the access to free space score based on biofilm images from experimental assays. The main file is "image_analysis_AFS.m" and requires Matlab's Image Processing Toolbox. 
Author: Lukas Eigentler (leigentler001@dundee.ac.uk). Please note the separate license file for the Matlab function "DeltaE_mod.m".

More information can be found directly in the code scripts or in readme files contained within the subfolders.


