# Kernel-embedded-Functional-Subspace-Detector
Attacks Detector using Kernel Functional Subspaces Theory
** Read me **

** kernel-embedded functional subspace detector (KEFSD) **

** Implementation was done in MATLAB R2024a and requires matpower**

** Author: Andrew Omiloli
** Date: 6/16/2026



1.  Generate clean residuals (attack free) from estimator or physical plant modelled in Simulink and save as mat file e.g 'residual_nominal.mat'
      > comment the detector, attack input and run the model 'Model_with_Detector.slx' to generate the clean residuals 

2. Set number of training samples (n) and number of functional basis elements (n_basis) in the run file, 'Run_file.m'
       > n = 200; n_basis = 20;
    > Note:  n_basis <= n

3. Set window length (Window), thresholds (tau_KEFSD and tau_norm) in 'Run_file.m'
       > Window = 20;
       > tau_KEFSD = 0.03
       > tau_norm  = 0.01;
    > Note1: 'tau_norm' is the 2norm residual detector used before window samples is full while 'tau_KEFSD' is the threshold used after full window samples
    > Note2: These thresholds should be tuned for each system to obtain reasonable detector performance metrics


4. Optional: Other parameters used for learning the principal functional basis elements  can be changed in the run file, for example
             > gamma_min = 1e-6
             > gamma_max  = 10
             > n_gamma    = 30
             > lambda_pca = 1e-3

5. Run the file 'Run_file.m'


6. open the model 'Model_with_Detector.slx' and uncomment the detector blocks and attack blocks and click run, if successfully run, the following signals gets logged in the workspace  
           > detector_output
           > alarm 
           > ground_truth
	> Note1: the attacks should be enabled at x duty cycle. 

7. Run the file 'Detector_Performance.m' to see detector performance scores, if undesirable, change thresholds to improve performance 


%% How to tune thresholds 
1. Test with random attacks enabled at x duty cycle
2. Examine the detector output signal at the region the detector energy flattens 
3. Choose desired threshold before the flat occurs 
4. Update thresholds in 'Run_file.m' and click run
5. Run 'Detector_Performance.m'
6. Repeat as necessary. 



