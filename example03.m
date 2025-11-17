function example03

RF_SHAPE = 'nsinc';
FA_EXC_IDEAL = 145; % degrees
numTRs = 27; % steady state

fprintf('If mex''ed, this should take <2m on a modern computer\n\n');
tic, blochsim_1D_2PULSE(RF_SHAPE, FA_EXC_IDEAL, numTRs); toc

fprintf('\n');
fprintf('Fig 1: Excitation (145 deg) and Refocusing SINC pulses\n');
fprintf('Fig 2: Transverse magnetization vs z at steady state\n');
fprintf('Fig 3: Longitudinal magnetization vs z at steady state\n');
