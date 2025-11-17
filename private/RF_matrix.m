function M = RF_matrix(Delta_w, RF, dt)
% implements Appendix B of Mulkern and Williams (1993)
% (w --> Delta_w, w1 --> RF, b --> beta)
%
% NOTE: RF is assumed to be real-valued, corresponding to
% a pulse along x in the rotating frame
%
% (BoB Mulkern and Mukund Balasubramanian, 2025)

if ~isreal(RF)
  error('a real-valued RF was expected!')
end

M = zeros(3,3);

beta = sqrt(Delta_w^2 + RF^2);

% should we be using abs and EPS below?
if (Delta_w == 0) & (beta == 0)
  Delta_w_over_beta = 1;
else
  Delta_w_over_beta = Delta_w / beta;
end

% should we be using abs and EPS below?
if (RF == 0) & (beta == 0)
  RF_over_beta = 1;
else
  RF_over_beta = RF / beta;
end

M(1,1) = 1 - (Delta_w_over_beta)^2 * (1-cos(beta*dt));
M(2,2) = cos(beta*dt);
M(3,3) = 1 - RF_over_beta^2 * (1-cos(beta*dt));

M(1,2) = Delta_w_over_beta * sin(beta*dt);
M(2,1) = -M(1,2);

M(1,3) = RF_over_beta * Delta_w_over_beta * (1-cos(beta*dt));
M(3,1) = +M(1,3);

M(2,3) = RF_over_beta * sin(beta*dt);
M(3,2) = -M(2,3);

return
