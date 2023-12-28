# btd_deconv
*A deconvolution approach for functional neuroimaging data to jointly estimate HRFs and source signals*

This package reproduces the results of the paper: 
Erol, A., Soloukey, C., Generowicz, B. et al. Deconvolution of the Functional Ultrasound Response in the Mouse Visual Pathway Using Block-Term Decomposition. Neuroinform (2022). https://doi.org/10.1007/s12021-022-09613-3

All codes are written in MATLAB. Tensorlab should be installed (https://www.tensorlab.net) for using this repository. 

Within this package, it is possible to test the proposed algorithm both over simulated data (under the folder *Simulations [+Demo]*) and real functional ultrasound time-series (under the folder *Real Data*). After adding Tensorlab to the search path of MATLAB, the test codes can be run without any extra adjustments, i.e., the structure of the repository should be kept as it is.

The software is so far tested on Windows and macOS.

#

### Table of Contents

* Simulations [+Demo]: Here, it is possible to *(i)* run a simple demo over example simulated data ('simple_demo.m') to quickly go over the steps of the algorithm and observe what are the expected inputs/outputs, and *(ii)* run all the Monte-Carlo simulations ('run_all_simulations.m').

* Real Data: Contains two sub-folders dedicated to analyzing two different functional ultrasound datasets, one with a single-stimulus condition and one with multiple stimulus conditions.

* utils: This folder has the implementation of the proposed method. The main algorithm is implemented within the function 'btd_deconv.m', which calls some of the other functions in the folder such as the tensor structurings. The remaining of the functions are used for visualizing the results and are called within the test codes.

* Results: Contains the results (i.e., plots) from both simulations and functional ultrasound experiments as presented in the paper.  

#

Contact a.kazaz@tudelft.nl for any questions.
