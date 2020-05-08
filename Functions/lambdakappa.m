function [lambda_vector,kappa_vector] = lambdakappa(time,Coef)


lambda_of_t = @(t) ((Coef(4) + (((Coef(6)*Coef(4)) - Coef(4))*(1 -exp(-Coef(5)*t)))));
lambda_vector = feval(lambda_of_t,time);

kappa_of_t = @(t) ((Coef(7)/Coef(9)) + (((Coef(7)-(Coef(7)/Coef(9))))*(exp(-Coef(8)*t))));
kappa_vector = feval(kappa_of_t,time);

end

