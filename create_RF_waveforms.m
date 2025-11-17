function [rfEncOut, rfOtherOut, dt, tbG] = create_RF_waveforms(RF_SHAPE, FA_EXC)
% [rfEncOut, rfOtherOut, dt, tbG] = create_RF_waveforms(RF_SHAPE, FA_EXC)
%
% INPUTS:
% -- RF_SHAPE: 'slr' or 'nsinc' (hanning-windowed sinc)
% -- FA_EXC: flip angle (in degrees) for the excitation pulse
%
% OUTPUTS:
% -- rfEncOut: EXCITATION pulse waveform
% -- rfOtherOut: REFOCUSING pulse waveform
% -- dt: dwell time (in ms)
% -- tbG: bandwidth-time product of the excitation (and refocus) pulse
%
% (Will Grissom and Mukund Balasubramanian, 2025)

% number of points used in SLR algoritm
N = 256;

% time-bandwidth product ...
tbG = 8;     % ... for excitation pulse
tbOther = 8; % ... for refocusing pulse

% duration of RF pulses (in ms)
T = 10;

% dwell time for RF pulses (in ms)
dt = 4e-3;

% desired slice thickness (in mm)
TH = 8;

% const
GAMMA_BAR = 42.58*1e6; % Hz/T

% calc Nout, the number of points in the final EXC waveform 
Nout = round(T/dt);
T = Nout * dt; % adjust T slightly for discretization

% for EXC, we now only need to solve for the gradient amplitude:
GS_EXC = (tbG/(T*1e-3)) / (GAMMA_BAR * TH*1e-3);
GS_EXC = GS_EXC * 1e3; % T/m --> mT/m

% for now, we choose to set GS_REF equal to GS_EXC
GS_REF = GS_EXC;
  
% for REF, we now only need to solve for its duration:
Tother = tbOther / (GAMMA_BAR * GS_REF*1e-3 * TH*1e-3) * 1e3;

% calc NoutOther, the number of points in the final REF waveform 
NoutOther = round(Tother/dt);
Tother = NoutOther * dt; % adjust Tother slightly for discretization

% slightly adjust GS_REF for the discretization above
GS_REF = (tbOther/(Tother*1e-3)) / (GAMMA_BAR * TH*1e-3);
GS_REF = GS_REF * 1e3; % T/m --> mT/m

% print statements:
switch RF_SHAPE
  case 'slr',
    fprintf('SLR: FA_EXC = %.1f degs N = %d\n', FA_EXC, N);
  case 'nsinc', 
    fprintf('hanning-windowed sinc: FA_EXC = %.1f degs N = %d\n', FA_EXC, N);
  otherwise,
    error('invalid RF_SHAPE entered');
end
fprintf('TBW_EXC = %.1f TBW_REF = %.1f\n', tbG, tbOther);
fprintf('T_EXC = %.3f ms T_REF = %.3f ms\n', T, Tother);
fprintf('BW_EXC = %.1f Hz BW_REF = %.1f Hz\n', ...
        tbG/(T*1e-3), tbOther/(Tother*1e-3));
fprintf('TH = %.1f mm\n', TH);
fprintf('GS_EXC = %.1f mT/m GS_REF = %.1f mT/m (ratio: %.3f)\n', ...
        GS_EXC, GS_REF, GS_EXC/GS_REF);
fprintf('Nout = %.3f NoutOther = %.3f\n', T/dt, Tother/dt);

% add path to RF tools
addpath rf_tools/ % JP's tools: gets dinf, b2a, cabc2rf, abr...
addpath rf_tools/mex5/

% ------------------------------------------------------------ 
% EXCITATION PULSE
% ------------------------------------------------------------ 

switch RF_SHAPE

  case 'slr', 
    
    exFlip = FA_EXC*pi/180; % excitation flip angle
    bsf = sin(exFlip/2); % for excitation pulse

    % in-slice and out-of-slice magnetization ripples
    d1 = 0.01; d2 = 0.01; 

    % now calculate the corresponding polynomial ripples
    d1Candidates = -0.1:0.001:0.1;
    bd1 = (1 + d1Candidates)*sin(exFlip/2);
    ad1 = sqrt(1-bd1.^2);
    MxyErr = abs(2*ad1.*bd1 - 2*cos(exFlip/2)*sin(exFlip/2)); 
    [~,NegInd] = min(abs(MxyErr(1:100) - d1));
    Negd1 = -d1Candidates(NegInd);
    [~,PosInd] = min(abs(MxyErr(101:end) - d1));
    Posd1 = d1Candidates(100+PosInd);
    fprintf('desired magnetization ripple d1 = %.4f\n', d1);
    d1 = max(Negd1,Posd1);
    fprintf('polynomial ripple d1: %.4f for %.1f deg FA\n', d1, FA_EXC);
    d2 = d2/sqrt(2); % stopband ripple: independent of flip angle

    % design the filter
    ftw = dinf(d1,d2)/tbG; % fractional transition width
    f_firls = [0 (1-ftw)*(tbG/2) (1+ftw)*(tbG/2) (N/2)]/(N/2);
    m_firls = [1 1 0 0];
    w_firls = [1 d1/d2];
    b_firls = firls(N-1,f_firls,m_firls,w_firls); % the filter
    b = bsf*b_firls;

    a = b2a(b);
    target = fft(b(:),8*N).*exp(1i*angle(fft(flipud(a(:)),8*N)));
    b = lsqr(@(x,tflag)fftoversamp(x,N,8,tflag),target);

    rfEnc = cabc2rf(a,fliplr(b(:).')); % iSLR for min-power pulse
    rfEnc = rfEnc(:);

  case 'nsinc', 

    % hanning-windowed sinc
    tt = linspace(-T/2, +T/2, N+1);
    tt(end) = [];
    tt = tt(:);
    rfEnc = sinc(tbG*tt/T);
    rfEnc = rfEnc .* (1+cos(2*pi*tt/T))/2;
    rfEnc = rfEnc / sum(rfEnc) * (FA_EXC/180*pi);

end

% ------------------------------------------------------------ 
% REFOCUS PULSE
% ------------------------------------------------------------ 

switch RF_SHAPE

  case 'slr', 

    % in-slice and out-of-slice magnetization ripples
    d1_Other = 0.01; d2_Other = 0.01;

    % design the refocus pulse:
    rfOther = dzrf(N, tbOther, 'se', 'ls', d1_Other, d2_Other);
    rfOther = rfOther(:);

  case 'nsinc', 

    % hanning-windowed sinc
    tt = linspace(-T/2, +T/2, N+1);
    tt(end) = [];
    tt = tt(:);
    rfOther = sinc(tbOther*tt/T);
    rfOther = rfOther .* (1+cos(2*pi*tt/T))/2;
    rfOther = rfOther / sum(rfOther) * pi;
    
end

% ------------------------------------------------------------ 
% RESAMPLE RFs (NOTE: this does nothing if the input N = T/dt)
% ------------------------------------------------------------ 

rfEncOut = interp1((0:N)'/N*T, [rfEnc; rfEnc(end)], ...
                   (0:Nout-1)'*dt, 'spline', 0);
rfEncOut = rfEncOut / abs(sum(rfEncOut)) * abs(sum(rfEnc)); 
rfEncOut = rfEncOut / (dt * 1e-3);

rfOtherOut = interp1((0:N)'/N*Tother, [rfOther; rfOther(end)], ...
                     (0:NoutOther-1)'*dt,'spline',0);
rfOtherOut = rfOtherOut / abs(sum(rfOtherOut)) * abs(sum(rfOther)); 
rfOtherOut = rfOtherOut / (dt * 1e-3);

% ------------------------------------------------------------ 
% EXPRESS RF AMPLITUDES IN UNITS OF MICROTESLA
% ------------------------------------------------------------ 

uT_factor = 2*pi*GAMMA_BAR*1e-6; % conversion factor to uT
rfEncOut_uT = rfEncOut / uT_factor;
rfOtherOut_uT = rfOtherOut / uT_factor;

% ------------------------------------------------------------ 
% PLOT RFs
% ------------------------------------------------------------ 

figure(1); clf
subplot(2,1,1);
ph = plot(linspace(0,T,Nout), real(rfEncOut_uT)); hold on
set(ph, 'linewidth', 2);
switch RF_SHAPE
  case 'slr', set(ph, 'color', 'r');
  otherwise,  set(ph, 'color', 'b');
end
plot([T/2 T/2], ylim, 'k--');
plot([0 T], [0 0], 'k-');
xlim([0 max([T Tother])]);
xlabel('time (ms)')
ylabel('B_1^+ (\muT)')
title(sprintf('excitation RF pulse (FA = %d^{\\circ})', FA_EXC));

subplot(2,1,2);
ph = plot(linspace(0,Tother,NoutOther), real(rfOtherOut_uT)); hold on
set(ph, 'linewidth', 2);
switch RF_SHAPE
  case 'slr', set(ph, 'color', 'r');
  otherwise,  set(ph, 'color', 'b');
end
plot([Tother/2 Tother/2], ylim, 'k--');
plot([0 Tother], [0 0], 'k-');
xlim([0 max([T Tother])]);
xlabel('time (ms)')
ylabel('B_1^+ (\muT)')
title(sprintf('refocusing RF pulse'));

return
