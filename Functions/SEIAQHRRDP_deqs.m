function deqs = SEIAQHRRDP_deqs(t, y, Coef, N, rho, gamma, nu, phi, delta, theta, tau, pi, SD_delay, SD_remove, zeta_factor, window, chi)

% ### ODEs
% dSdt   = -((chi_of_t*rho*beta*S)/N)*(A+I) - (alpha*alpha_ON*S) + (zeta*zeta_of_t*P) + ((pi/(1-tau)))*(Ra + Q + Rh);
% dEdt   = ((chi_of_t*rho*beta*S)/N)*(A+I) - (gamma*E);
% dIdt   = (gamma*E) - (nu*I) - (theta*I) - (delta*I);
% dAdt   = (nu*I) - (phi*A);
% dQdt   = (theta*I) - (((pi/(1-tau)))*Q);
% dHdt   = (delta*I) - ((lambda + (((sigma*lambda) - lambda)*(1 - exp(-epsilon*t))))*H)...
%           - (((kappa/sigma) + ((PF*(kappa-(kappa/sigma)))*(exp(-psi*t))) + (((1-PF)*(kappa-(kappa/sigma)))*(exp(-mu*t))))*H);
% dRadt = (phi*A) - (((pi/(1-tau)))*Ra);
% dRhdt  = ((lambda + (((sigma*lambda) - lambda)*(1 -exp(-epsilon*t))))*H)...
%           - (((pi/(1-tau)))*Rh);
% dDdt   = (((kappa/mu) + (((kappa-(kappa/mu)))*(exp(-psi*t))))*H);
% dPdt   = (alpha*alpha_ON*S) - (zeta*zeta_of_t*P);

% ### States to solve
% y(1)  = S:   number of succeptibles
% y(2)  = E:   number of exposed
% y(3)  = I:   number of infectious (not quarantined)
% y(4)  = A:   number of infectious asymptomatic
% y(5)  = Q:   number of quarantined, active cases not requiring hospitalization
% y(6)  = H:   number of hospitalized, active cases requiring hospitalization
% y(7)  = Ra:  number of recovered asymptomatic cases
% y(8)  = Rh:  number of recovered cases requiring hospitalization
% y(9)  = D:   number of dead
% y(10) = P:   number of protected (by public health measures, ie social distancing)

% ### Coefficients to solve
% Coef(1) = alpha: effective protection rate
% Coef(2) = zeta: effective protection leak rate
% Coef(3) = beta: effective contact rate
% Coef(4) = lambda: basic intervention recovery rate
% Coef(5) = epsilon: recovery intervention development rate constant;
% Coef(6) = sigma: intervention recovery efficacy;
% Coef(7) = kappa: basic intervention death rate
% Coef(8) = psi: death intervention development rate constant
% Coef(9) = mu: intervention death efficacy;

% ### Fixed Coefficients
% chi_of_t: seasonal variability in infectivity modeled after influenza from 1997-present
% rho: population density
% gamma: effective exposed to infected rate
% nu: effective asymptomatic rate
% phi: (infectious asymptomatic period)^-1
% delta: effective hospitalization rate
% theta: effective infected to quarantined rate
% PF: Fraction of death intervention span accounted by fast rate constant
% tau: fraction of recovered who develop immunity
% pi: (duration of immunity)^-1
% SD_delay: time to social distancing measures
% SD_remove: time to remove social distancing measures from start of measures

%% Import Flu data to estimate seasonal variance in infectivitiy
% From CDC FluView Website
% https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html
% Averaged "WeightedILI" column over 1997-2020

chi_of_t = interp1(1:length(chi),chi,t);
%% Define alpha(t): Enable "P" once social distancing measures are enacted
alpha_t = window(1):window(2);
alpha_fxn(1,window(1):window(1)+SD_delay) = 0;
alpha_fxn(1,window(1)+SD_delay:window(2)) = 1;
alpha_ON = interp1(alpha_t, alpha_fxn, t);

%% Define zeta(t): Determine when social distancing measures are lessened
zeta_t = window(1):window(2);
zeta_fxn=[];

if size(SD_remove,2) > 1
    zeta_fxn(1,window(1):window(1)+SD_delay+SD_remove(1)) = 1;
    for i = 2:size(zeta_factor,2)
        zeta_fxn = [zeta_fxn,ones(1,SD_remove(i))*zeta_factor(i-1)];
    end
    zeta_fxn = [zeta_fxn,ones(1,window(2)-(window(1)+SD_delay+sum(SD_remove)))*zeta_factor(end)];
else
    zeta_fxn(1,window(1):window(1)+SD_delay+SD_remove) = 1;
    zeta_fxn(1,window(1)+SD_delay+SD_remove:window(2)) = zeta_factor;
end

zeta_of_t = interp1(zeta_t, zeta_fxn, t);

%% Define ODEs
dSdt   = -((chi_of_t*rho*Coef(3)*y(1))/N)*(y(4)+y(3)) - (Coef(1)*alpha_ON*y(1)) + (Coef(2)*zeta_of_t*y(10)) + ((pi/(1-tau)))*(y(7) + y(5) + y(8));
dEdt   = ((chi_of_t*rho*Coef(3)*y(1))/N)*(y(4)+y(3)) - (gamma*y(2));
dIdt   = (gamma*y(2)) - (nu*y(3)) - (theta*y(3)) - (delta*y(3));
dAdt   = (nu*y(3)) - (phi*y(4));
dQdt   = (theta*y(3)) - ((pi/(1-tau)))*y(5);
dHdt   = (delta*y(3)) - ((Coef(4) + (((Coef(6)*Coef(4)) - Coef(4))*(1 -exp(-Coef(5)*t))))*y(6))...
          - (((Coef(7)/Coef(9)) + (((Coef(7)-(Coef(7)/Coef(9))))*(exp(-Coef(8)*t))))*y(6));
dRadt  = (phi*y(4)) - (((pi/(1-tau)))*y(7));
dRhdt  = ((Coef(4) + (((Coef(6)*Coef(4)) - Coef(4))*(1 -exp(-Coef(5)*t))))*y(6))...
          - (((pi/(1-tau)))*y(8));
dDdt   = (((Coef(7)/Coef(9)) + (((Coef(7)-(Coef(7)/Coef(9))))*(exp(-Coef(8)*t))))*y(6));
dPdt   = (Coef(1)*alpha_ON*y(1)) - (Coef(2)*zeta_of_t*y(10));
      
deqs = [dSdt; dEdt; dIdt; dAdt; dQdt; dHdt; dRadt; dRhdt; dDdt; dPdt];
end

