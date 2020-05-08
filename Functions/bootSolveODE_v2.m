function [out] = bootSolveODE_v2(Coef2_UB, tspan, Npop, ic, rho, gamma, nu, phi, delta, theta, tau, pi, SD_delay, SD_remove, zeta_factor, chi, UB, LB, yvals)

%% Resample NYSDOH from lognormal distribution to estimate CI
yvals_with_Noise = zeros(size(yvals));

%yvals = [Not_Hospitalized'; Hospitalized_active'; Discharged'; Deaths'];
% Not Hospitalized (estimated standard error = 0.2)
% Discharged (estimated standard error = 0.2)
% Deaths (estimated standard error = 0.2)
% These are cumulative, add noise where yval(n+1)>yval(n)
% for i = [1,3,4]
%     for j = 1:size(yvals,2)
%         if yvals(i,j) == 0
%             yvals_with_Noise(i,j) = yvals(i,j);
%         else
%             mean = yvals(i,j);
%             variance = (.1*yvals(i,j))^2;
%             mu_temp = log((mean^2)/sqrt(variance+mean^2));
%             sigma_temp = sqrt(log(variance/(mean^2)+1));
%             noise_temp = lognrnd(mu_temp,sigma_temp);
%             while noise_temp < yvals_with_Noise(i,j-1)
%                 mean = yvals(i,j);
%                 variance = (.1*yvals(i,j))^2;
%                 mu_temp = log((mean^2)/sqrt(variance+mean^2));
%                 sigma_temp = sqrt(log(variance/(mean^2)+1));
%                 noise_temp = lognrnd(mu_temp,sigma_temp);
%             end
%             j
%             yvals_with_Noise(i,j) = noise_temp;
%         end
%     end
% end
% 
% Hospitalized active (estimated standard error = 0.2)
% Is not cumulative, just generate log gaussian noise
% for j = 1:size(yvals,2)
%     if yvals(2,j) == 0
%         yvals_with_Noise(2,j) = yvals(2,j);
%     else
%         mean = yvals(2,j);
%         variance = (.2*yvals(2,j))^2;
%         mu_temp = log((mean^2)/sqrt(variance+mean^2));
%         sigma_temp = sqrt(log(variance/(mean^2)+1));
%         yvals_with_Noise(2,j) = lognrnd(mu_temp,sigma_temp);
%     end
% end
% figure;plot(yvals','r');hold on;plot(yvals_with_Noise','b')
% %%
% %yvals = [Not_Hospitalized'; Hospitalized_active'; Discharged'; Deaths'];
% % Not Hospitalized (estimated standard error = 0.2)
% % Discharged (estimated standard error = 0.2)
% % Deaths (estimated standard error = 0.2)
% % These are cumulative, add noise where yval(n+1)>yval(n)
% 
% for i = [1,3,4]
%     for j = 1:size(yvals,2)
%         if yvals(i,j) == 0
%             yvals_with_Noise(i,j) = yvals(i,j);
%         else
%             mean_temp = yvals(i,j);
%             sigma_temp = .1*yvals(i,j);
%             noise_temp = normrnd(mean_temp,sigma_temp);
%             while noise_temp < yvals_with_Noise(i,j-1)
%                 mean_temp = yvals(i,j);
%                 sigma_temp = .1*yvals(i,j);
%                 noise_temp = normrnd(mean_temp,sigma_temp);
%             end
%             yvals_with_Noise(i,j) = noise_temp;
%         end
%     end
% end
% 
% % Hospitalized active (estimated standard error = 0.2)
% % Is not cumulative, just generate log gaussian noise
% for j = 1:size(yvals,2)
%     if yvals(2,j) == 0
%         yvals_with_Noise(2,j) = yvals(2,j);
%     else
%         mean_temp = yvals(2,j);
%         sigma_temp = .1*yvals(2,j);
%         noise_temp = normrnd(mean_temp,sigma_temp);
%         yvals_with_Noise(2,j) = noise_temp;
%     end
% end
% % figure;plot(yvals','r');hold on;plot(yvals_with_Noise','b')

%%
%yvals = [Not_Hospitalized'; Hospitalized_active'; Discharged'; Deaths'];
% Not Hospitalized (estimated standard error = 0.2)
% Discharged (estimated standard error = 0.2)
% Deaths (estimated standard error = 0.2)
% These are cumulative, add noise to the differential, then calculate
% cumulative sum
yvals_diff = diff(yvals,1,2);
yvals_diff = [zeros(4,1) yvals_diff];

for i = [1,3,4]
    for j = 1:size(yvals,2)
        if yvals(i,j) == 0
            yvals_with_Noise(i,j) = yvals(i,j);
        else            
            mean_temp = yvals_diff(i,j);
            if mean_temp == 0
                yvals_with_Noise(i,j) = noise_temp + yvals_with_Noise(i,j-1);
            else
                variance = (.7*yvals_diff(i,j))^2;
                mu_temp = log((mean_temp^2)/sqrt(variance+mean_temp^2));
                sigma_temp = sqrt(log(variance/(mean_temp^2)+1));
                noise_temp = lognrnd(mu_temp,sigma_temp);
                yvals_with_Noise(i,j) = noise_temp + yvals_with_Noise(i,j-1);
            end
        end
    end
end

% Hospitalized active (estimated standard error = 0.2)
% Is not cumulative, just generate log gaussian noise
for j = 1:size(yvals,2)
    if yvals(2,j) == 0
        yvals_with_Noise(2,j) = yvals(2,j);
    else
        mean_temp = yvals(2,j);
        variance = (.1*yvals(2,j))^2;
        mu_temp = log((mean_temp^2)/sqrt(variance+mean_temp^2));
        sigma_temp = sqrt(log(variance/(mean_temp^2)+1));
        noise_temp = lognrnd(mu_temp,sigma_temp);
        yvals_with_Noise(2,j) = noise_temp;
    end
end
figure;plot(yvals','r');hold on;plot(yvals_with_Noise','b')

%%


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
obj = sum(sum((myfcn -  yvals_with_Noise).^2));
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

