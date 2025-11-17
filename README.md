# SLR_shortTR_SE
MATLAB code accompanying Balasubramanian et al. (2026) "Large gains in SNR through the application of Shinnar-Le Roux RF pulse design to short-TR spin-echo fMRI acquisitions at 7T"

1) This code was written and tested on MATLAB R2015a on Linux; if you are using a very different version of MATLAB or a different OS, you may have to make small changes to the code (e.g., changing the file separator character from "/" to "\", etc.)

2) John Pauly's "rf_tools" MATLAB library is required; it can be downloaded from https://rsl.stanford.edu/research/software.html (scroll down to "RF Pulse Design", below which is the "Download Matlab library" link). If you already have rf_tools installed, or wish to install it somewhere other than in the current directory, please edit the two lines starting with "addpath" in create_RF_waveforms.m accordingly. (Non-Linux users may also need to edit the second "addpath" line, depending on their OS.)

3) In the MATLAB command window, cd to the "private" directory and type: "mex innerloop_blochsim_1D_2PULSE_MEX.c". Note that if this compilation fails, the parent code will still run (by instead calling the unmex'ed version of the function innerloop_blochsim_1D_2PULSE.m), but will be painfully slow!

4) Try running example01, example02 and example03 (noting the expected run times for the mex'ed version of the code).

5) If you find this code useful in your work, please don't forget to cite the accompanying manuscript (Balasubramanian et al., "Large gains in SNR through the application of Shinnar-Le Roux RF pulse design to short-TR spin-echo fMRI acquisitions at 7T". Magn Reson Imaging 2026).
