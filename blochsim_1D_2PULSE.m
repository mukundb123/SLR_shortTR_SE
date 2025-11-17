function blochsim_1D_2PULSE(RF_SHAPE, FA_EXC_IDEAL, numTRs)
% blochsim_1D_2PULSE(RF_SHAPE, FA_EXC_IDEAL, numTRs)
%
% INPUTS:
% -- RF_SHAPE: 'slr' or 'nsinc' (hanning-windowed sinc)
% -- FA_EXC_IDEAL: flip angle (in degrees) for the excitation pulse
% -- numTRs: number of repetitions to simulate
%
% calls (in current directory) create_RF_waveforms and (in private/)
% innerloop_blochsim_1D_2PULSE(_MEX)
%
% (Mukund Balasubramanian, 2025)

switch RF_SHAPE
  case 'slr',
    USE_SINC = false;
  case 'nsinc', 
    USE_SINC = true;
  otherwise,
    error('invalid RF_SHAPE entered');
end

DOWNSAMPLE_RF = 1; % increase this for a speedup (with lower accuracy)
PTS_PER_VOX = 250; % decrease this for a speedup (with lower accuracy)
NUM_CYCLES_VOX = 16; % used to determine the strength of crushing
B1_factor = 1; % non-unity values simulate imperfect B1+ values

% constants
gamma_bar = 4258;
gamma = 2*pi*gamma_bar;

% params
TH = 8; % slice thickness in mm (for now, same for EXC and REF)
VOX_SIZE = 0.250; % size of voxels (in mm) along the slice direction
dz = VOX_SIZE/PTS_PER_VOX; % size of fine-grid subvoxels (in mm)

T1 = 1500;  % ms
T2 = 50;    % ms
TE = 30;    % ms
TR = 300;   % ms
tau = TE/2; % ms

E1 = exp(-tau/T1);
E1p = exp(-(TR-TE)/T1);

% ------------------------------------------------------------ 
% compute ideal case
% ------------------------------------------------------------ 
mm_INIT = 1;
for ii = 1 : numTRs
  gg_IDEAL = mm_INIT * sin(FA_EXC_IDEAL/180*pi);
  mm_IDEAL = mm_INIT * cos(FA_EXC_IDEAL/180*pi);
  mm_IDEAL = (1-E1) + mm_IDEAL * E1; % relax: EXC --> TE/2
  gg_IDEAL = -gg_IDEAL; % ideal REF
  mm_IDEAL = -mm_IDEAL; % ideal REF
  mm_IDEAL = (1-E1) + mm_IDEAL * E1; % relax: TE/2 --> TE
  mm_INIT = (1-E1p) + mm_IDEAL * E1p; % relax: TE --> TR end
end

% ------------------------------------------------------------ 
% create both RF_EXC and RF_REF
% ------------------------------------------------------------ 
if USE_SINC
  [RF_EXC, RF_REF, dt, tbG] = create_RF_waveforms('nsinc', FA_EXC_IDEAL);
else
  [RF_EXC, RF_REF, dt, tbG] = create_RF_waveforms('slr', FA_EXC_IDEAL);
end
% NOTE: create_RF_waveforms currently assumes that EXC and REF have
% the same length, dwelltime and bandwidth-time product

% right now, we can only handle real-valued RFs
if max(abs(imag(RF_EXC))) > sqrt(eps)
  error('RF_EXC should be real-valued');
else
  RF_EXC = real(RF_EXC);
end

% right now, we can only handle real-valued RFs
if max(abs(imag(RF_REF))) > sqrt(eps)
  error('RF_REF should be real-valued');
else
  RF_REF = real(RF_REF);
end

% waveform params
N_EXC = length(RF_EXC);
N_REF = length(RF_REF);
T_EXC = N_EXC * dt; % ms
T_REF = N_REF * dt; % ms
BWT_EXC = tbG;
BWT_REF = tbG;

% CRUSH_FACTOR in terms of NUM_CYCLES
NUM_CYCLES_SLICE = NUM_CYCLES_VOX * (TH / VOX_SIZE);
CRUSH_FACTOR = NUM_CYCLES_SLICE / BWT_EXC;

% ms --> sec
dt = dt / 1000;       % sec
T_EXC = T_EXC / 1000; % sec
T_REF = T_REF / 1000; % sec

% NB: the SLR slice-center FA doesn't exactly match the requested FA
FA_EXC = abs(sum(RF_EXC)*dt)/pi*180;
FA_REF = abs(sum(RF_REF)*dt)/pi*180;

% let's downsample the RF to speedup the sims
dt = dt * DOWNSAMPLE_RF;
tt = (0 : dt : T_EXC-dt)';
RF_EXC = RF_EXC(1:DOWNSAMPLE_RF:N_EXC);
RF_REF = RF_REF(1:DOWNSAMPLE_RF:N_REF);

% normalize amplitude to achieve desired FA (at slice center)
RF_EXC = RF_EXC / abs(sum(RF_EXC)*dt) * (FA_EXC/180*pi);
RF_REF = RF_REF / abs(sum(RF_REF)*dt) * (FA_REF/180*pi);
% so sum(RF_EXC or RF_REF) * dt = desired_FA (in radians)

% apply B1 inhomogeneity factor
RF_EXC = RF_EXC * B1_factor;
RF_REF = RF_REF * B1_factor;

% TH in cm, DUR in s --> grad in gauss/cm
grad = BWT_EXC / (gamma_bar * (TH/10) * T_EXC);
% for now, we assume that EXC and REF have the same gradient AMPLITUDE
% (which they will if they have the same BWT, TH and DURATION)

% define the subvoxel z-grid
zz_list = (-TH : dz : TH)'/10; % divide by 10 for mm --> cm

try
% try calling the mex'ed innerloop ...
[ff_list, gg_list, mm_list, mm_init_list] = ...
 innerloop_blochsim_1D_2PULSE_MEX(zz_list, grad, dt, RF_EXC, ...
  numTRs, E1, CRUSH_FACTOR, RF_REF, E1p);
catch
% ... if that fails, call the UNmex'ed innerloop
fprintf('WARNING: could not find innerloop_blochsim_1D_2PULSE_MEX\n');
fprintf('instead running the unmexed innerloop_blochsim_1D_2PULSE\n');
[ff_list, gg_list, mm_list, mm_init_list] = ...
 innerloop_blochsim_1D_2PULSE(zz_list, grad, dt, RF_EXC, ...
  numTRs, E1, CRUSH_FACTOR, RF_REF, E1p);
end

% downsample: integrate signal from points within each voxel
zz_down = [];
ff_down = [];
gg_down = [];
mm_down = [];
counter = 0;
for zz_mm = -TH : VOX_SIZE : TH-VOX_SIZE
  counter = counter + 1;
  zz_cm_start = zz_mm/10;
  zz_cm_end = (zz_mm+VOX_SIZE)/10;
  indx = find( (zz_list >= zz_cm_start) & (zz_list < zz_cm_end) );
  zz_down(counter) = mean(zz_list(indx));
  ff_down(counter) = mean(ff_list(indx));
  gg_down(counter) = mean(gg_list(indx));
  mm_down(counter) = mean(mm_list(indx));
end
zz_list = zz_list(:)*10; % convert from cm to mm
zz_down = zz_down(:)*10; % convert from cm to mm

% plot slice profile for Mx and My
figure(2); clf
ph2 = plot(zz_down, -ff_down, 'r'); hold on
set(ph2, 'linewidth', 2);
ph1 = plot(zz_down, -gg_down, 'b');
set(ph1, 'linewidth', 2);
ph = plot([-TH -TH/2 -TH/2     +TH/2    +TH/2 +TH], ...
            [0  0    -gg_IDEAL -gg_IDEAL 0     0], 'k--');
set(ph, 'linewidth', 1.5);
plot([-TH +TH], [0 0], 'k:');
xlim([-TH +TH])
legend('Mx', 'My', 'ideal', 'location', 'northwest');
xlabel('z (mm)')
ylabel('transverse magnetization (a.u.)');
title('slice profile');

% plot slice profile for Mz
figure(3); clf
ph = plot(zz_down, mm_down, 'm'); hold on
set(ph, 'linewidth', 2);
ph = plot([-TH -TH/2 -TH/2     +TH/2    +TH/2 +TH], ...
            [1  1    +mm_IDEAL +mm_IDEAL 1     1], 'k--');
set(ph, 'linewidth', 1.5);
plot([-TH +TH], [0 0], 'k:');
xlim([-TH +TH])
legend('Mz', 'ideal', 'location', 'northwest')
xlabel('z (mm)')
ylabel('longitudinal magnetization (a.u.)');
title('slice profile');

return
