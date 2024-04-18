function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively

Ct = [eye(6) zeros(6,9)];
Rt = eye(6)*0.000000095;    % Trust Factor

K = (covarEst * Ct')*pinv((((Ct * covarEst * Ct') + Rt))); %Kalman gain 
%% Update 
uCurr = uEst + (K * (z_t - (Ct * uEst))); %mean current
covar_curr = covarEst - (K * Ct * covarEst); %covariance current 



end