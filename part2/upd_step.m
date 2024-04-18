function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively


C = [zeros(3,6) eye(3) zeros(3,6)]; %C matrix

R = eye(3) * 0.000001; %Trust Factor
K = (covarEst * transpose(C))*pinv((((C * covarEst * transpose(C)) + R))); %Kalman gain 

%% Update Step
uCurr = uEst + (K * ( z_t - (C * uEst))); %mean current
covar_curr = covarEst - (K * C * covarEst); %covariance current 

end