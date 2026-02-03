var C1 C2 F1 F2 Y1 Y2 Y pF R Pi;
varexo eps_y1 eps_y2 eps_r eps_pf;

parameters pi1 pi2 beta1 beta2 sigma1 sigma2 h1 h2 omega1 omega2 barF1 barF2;
parameters rho_y1 rho_y2 Rbar pFbar rho_r rho_pf sigma_r sigma_pf sigma_y1 sigma_y2 Y1bar Y2bar;

@#include "micro_macro_params.m"

model;
  C1 + pF * F1 = Y1;
  C2 + pF * F2 = Y2;

  pF = (1-omega1)/omega1 * (C1 - h1*C1(-1)) / (F1 - barF1);
  pF = (1-omega2)/omega2 * (C2 - h2*C2(-1)) / (F2 - barF2);

  (C1 - h1*C1(-1))^(-sigma1) = beta1 * R / Pi(+1) * (C1(+1) - h1*C1)^(-sigma1);
  (C2 - h2*C2(-1))^(-sigma2) = beta2 * R / Pi(+1) * (C2(+1) - h2*C2)^(-sigma2);

  log(Y1) = (1-rho_y1)*log(Y1bar) + rho_y1*log(Y1(-1)) + sigma_y1*eps_y1;
  log(Y2) = (1-rho_y2)*log(Y2bar) + rho_y2*log(Y2(-1)) + sigma_y2*eps_y2;
  log(R) = (1-rho_r)*log(Rbar) + rho_r*log(R(-1)) + sigma_r*eps_r;
  log(pF) = (1-rho_pf)*log(pFbar) + rho_pf*log(pF(-1)) + sigma_pf*eps_pf;

  Pi = 1;
  Y = pi1*Y1 + pi2*Y2;
end;

initval;
  C1 = 0.7;
  C2 = 0.55;
  F1 = 0.3;
  F2 = 0.25;
  Y1 = Y1bar;
  Y2 = Y2bar;
  pF = pFbar;
  R = Rbar;
  Pi = 1;
  Y = pi1*Y1bar + pi2*Y2bar;
end;

steady;
check;
stoch_simul(order=1, periods=240);
