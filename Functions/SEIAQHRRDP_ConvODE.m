function output = SEIAQHRRDP_ConvODE(Coef_r, tspan, Npop, y0, rho, gamma, nu, phi, delta, theta, tau, pi, SD_delay, SD_remove, zeta_factor, chi)
    sol = ode45(@(t,y)SEIAQHRRDP_deqs(t, y, Coef_r, Npop, rho, gamma, nu, phi, delta, theta, tau, pi, SD_delay, SD_remove, zeta_factor, [tspan(1) tspan(end)], chi), tspan, y0);
    output = deval(sol,tspan);
end

