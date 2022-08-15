# btd_deconv
A deconvolution approach for functional neuroimaging data to jointly estimate HRFs and source signals

This package reproduces the results of the paper titled "Deconvolution of the Functional Ultrasound Response in the Mouse Visual Pathway Using Block-Term Decomposition" by Aybüke Erol, Chagajeg Soloukey, Bastian Generowicz, Nikki van Dorp, Sebastiaan Koekkoek, Pieter Kruizinga and Borbála Hunyadi (currently under review).

All codes are written in MATLAB. Tensorlab should be installed (https://www.tensorlab.net) for using this repository. 

The folder 'utils/' has the implementation of the proposed method. The main algorithm is 'btd_deconv.m', which calls some of the other functions in the folder, such as the tensor structurings. The remaining of the functions are used for visualizing the results and are called within the test codes.

The test codes are presented in two separate folders: one for simulations and one for real functional ultrasound data. After adding Tensorlab to the search path of MATLAB, the test codes can be run without any extra adjustments, i.e., the structure of the repository should be kept as it is.

In the folder 'Simulations [+Demo]', it is possible to (i) run a simple demo over example simulated data ('simple_demo.m') to quickly go over the steps of the algorithm and observe what are the expected inputs/outputs, and (ii) run all the Monte-Carlo simulations ('run_all_simulations.m').

In the folder 'Real Data', there are two sub-folders dedicated to analyzing two different functional ultrasound datasets, one with a single-stimulus condition and one with multiple stimulus conditions. 

Finally, the 'Results' folder has the results (i.e., plots) from both simulations and functional ultrasound experiments as presented in the paper.  

Contact a.erol@tudelft.nl for any questions.
