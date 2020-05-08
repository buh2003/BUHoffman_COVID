function output = SEIAQHRRDP_ConvODE2(Coef_r, tspan, Npop, y0, rho, gamma, nu, phi, delta, theta, tau, pi, SD_delay, SD_remove, zeta_factor, chi)
    solpts = SEIAQHRRDP_ConvODE(Coef_r, tspan, Npop, y0, rho, gamma, nu, phi, delta, theta, tau, pi, SD_delay, SD_remove, zeta_factor, chi);
    output = solpts([5,6,8,9],:); 
     
%     [lambda_vector,kappa_vector] = lambdakappa(tspan,Coef_r);
%     output = [(cumtrapz(solpts(4,:)).*theta);(cumtrapz(solpts(4,:)).*delta);(cumtrapz(solpts(4,:)).*lambda_vector); solpts(9,:)];
end