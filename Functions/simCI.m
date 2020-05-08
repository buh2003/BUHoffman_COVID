function bootyCI = simCI(bootCoef,Npop, rho, gamma, nu, phi, delta, theta, omega,tau, pi, SD_delay, SD_remove, zeta_factor, window_sim,chi, ic)
%% Estimate CI of simulation
booty = [] ;
for i = 1:size(bootCoef,1)
    [~,y] = ode45(@(t,y) SEIAQHRRRDP_deqs(t, y, bootCoef(i,:)', Npop, rho, gamma, nu, phi, delta, theta, omega,tau, pi, SD_delay, SD_remove, zeta_factor, window_sim,chi), window_sim(1):window_sim(2), ic);
    booty(:,:,i)=y;
end
bootyCI=[];
for i = 1: size(booty,2)
    for j = 1: size(booty,1)
        dist_data = squeeze(booty(j,i,:));
        SEM = std(dist_data)/sqrt(length(dist_data));               % Standard Error
        ts = tinv([0.025  0.975],length(dist_data)-1);      % T-Score
        CI = mean(dist_data) + ts*SEM*2; 
        bootyCI{j,i} =CI;
   end
end
end

