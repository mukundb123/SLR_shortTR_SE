function [ff_list, gg_list, mm_list, mm_init_list] = ...
         innerloop_blochsim_1D_2PULSE(zz_list, grad, dt, RF_EXC, ...
         numTRs, E1, CRUSH_FACTOR, RF_REF, E1p);
% [ff_list, gg_list, mm_list, mm_init_list] = ...
% innerloop_blochsim_1D_2PULSE(zz_list, grad, dt, RF_EXC, ...
% numTRs, E1, CRUSH_FACTOR, RF_REF, E1p);
%
% called by blochsim_1D_2PULSE (in the directory above)
%
% units:
% -- zz_list and TH: cm
% -- dt: sec
% -- grad: gauss/cm
%
% (Mukund Balasubramanian, 2025)

% constants
gamma_bar = 4258;
gamma = 2*pi*gamma_bar;

% params
SHIM_FACTOR = 0;

% allocs and initial conditions
ff_list = zeros(size(zz_list));
gg_list = zeros(size(zz_list));
mm_list = zeros(size(zz_list));
mm_init_list = ones(size(zz_list));

% loop through TRs
for ii = 1 : numTRs

% loop through z locations
for kk = 1 : length(zz_list)

  % get zz and fgm for this voxel
  zz = zz_list(kk);
  fgm = [0; 0; mm_init_list(kk)];

  % get the Delta_w from the slice-select gradient
  Delta_w = gamma * grad * (1+SHIM_FACTOR) * zz;

  % excitation RF pulse
  for jj = 1 : length(RF_EXC)
    M = RF_matrix(Delta_w, RF_EXC(jj), dt);
    fgm = M * fgm;
  end 
  
  % rewind: T_EXC = dt * length(RF_EXC)
  fg = fgm(1) + i*fgm(2);
  fg = fg * exp(i * gamma * grad * (1-SHIM_FACTOR) * zz * dt * length(RF_EXC) / 2);
  fgm(1) = real(fg);
  fgm(2) = imag(fg);

  % transverse relaxation between EXC and REF
  % (ignore for now)

  % longitudinal relaxation between EXC and REF
  fgm(3) = (1-E1) + fgm(3) * E1;

  % crush left before REF
  fg = fgm(1) + i*fgm(2);
  fg = fg * exp(i * gamma * grad * zz * dt * length(RF_EXC) * CRUSH_FACTOR);
  fgm(1) = real(fg);
  fgm(2) = imag(fg);
  
  % refocusing RF pulse
  for jj = 1 : length(RF_REF)
    M = RF_matrix(Delta_w, RF_REF(jj), dt);
    fgm = M * fgm;
  end 

  % crush right after REF
  fg = fgm(1) + i*fgm(2);
  fg = fg * exp(i * gamma * grad * zz * dt * length(RF_EXC) * CRUSH_FACTOR);
  fgm(1) = real(fg);
  fgm(2) = imag(fg);

  % transverse relaxation between REF and TE
  % (ignore for now)

  % longitudinal relaxation between REF and TE
  fgm(3) = (1-E1) + fgm(3) * E1;

  % store magnetization at TE for the kkth zz location
  ff_list(kk) = fgm(1);
  gg_list(kk) = fgm(2);
  mm_list(kk) = fgm(3);

  % Mz at end of TR will be the new Mz_init
  mm_init_list(kk) = (1-E1p) + fgm(3) * E1p;

end % matches for kk = 1 : length(zz_list)

end % matches for ii = 1 : numTRs

return
