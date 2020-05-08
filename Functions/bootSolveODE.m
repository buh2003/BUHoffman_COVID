function [out] = bootSolveODE(Coef2_UB, Coef, tspan, Npop, ic, rho, gamma, nu, phi, delta, theta, tau, pi, SD_delay, SD_remove, zeta_factor, chi, UB, LB, yvals)
UB(2) = Coef2_UB(1);
UB(2)

Coef_r = optimvar('Coef_r',9,"LowerBound",LB,"UpperBound",UB);

Coef1_randi = randi(round([LB(1), UB(1)]*1E5),1,1)*1E-5;
Coef2bounds = round([LB(2), UB(2)]*1E5);
Coef2_randi = randi(Coef2bounds,1,1)*1E-5;
Coef3_randi = randi(round([LB(3), UB(3)]*1E5),1,1)*1E-5;
Coef4_randi = randi(round([LB(4), UB(4)]*1E5),1,1)*1E-5;
Coef5_randi = randi(round([LB(5), UB(5)]*1E5),1,1)*1E-5;
Coef6_randi = randi(round([LB(6), UB(6)]*1E5),1,1)*1E-5;
Coef7_randi = randi(round([LB(7), UB(7)]*1E5),1,1)*1E-5;
Coef8_randi = randi(round([LB(8), UB(8)]*1E5),1,1)*1E-5;
Coef9_randi = randi(round([LB(9), UB(9)]*1E5),1,1)*1E-5;

Coefs_randi = [Coef1_randi,Coef2_randi,Coef3_randi,Coef4_randi,Coef5_randi,Coef6_randi,Coef7_randi,Coef8_randi,Coef9_randi];

myfcn = fcn2optimexpr(@SEIAQHRRDP_ConvODE2, Coef_r, tspan, Npop, ic, rho, gamma, nu, phi, delta, theta, tau, pi, SD_delay, SD_remove, zeta_factor, chi);
obj = sum(sum((myfcn - yvals).^2));
prob = optimproblem("Objective",obj);
constraint1 = Coef_r(7) <= Coef_r(4);
prob.Constraints.constraint1 = constraint1;

opts = optimoptions(@fmincon,'MaxIterations', 6e2,'MaxFunctionEvaluations',2e3,...
    'StepTolerance',1e-6,'OptimalityTolerance',1e-6,...
    'Algorithm','sqp');

Coef_r0.Coef_r = Coefs_randi;
[rsol,~,~,~,~] = solve(prob, Coef_r0,'Options',opts);
out = rsol.Coef_r;

end

