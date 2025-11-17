function example01

RF_SHAPE = 'slr';
FA_EXC_IDEAL = 145; % degrees
numTRs = 1; % one repetition

fprintf('If mex''ed, this should take <5s on a modern computer\n\n');
tic, blochsim_1D_2PULSE(RF_SHAPE, FA_EXC_IDEAL, numTRs); toc

fprintf('\n');
fprintf('Fig 1: Excitation (145 deg) and Refocusing SLR pulses\n');
fprintf('Fig 2: Transverse magnetization vs z after 1 repetition\n');
fprintf('Fig 3: Longitudinal magnetization vs z after 1 repetition\n');
