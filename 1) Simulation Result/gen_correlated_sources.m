function [S] = gen_correlated_sources(corr_coeff,T,Q)
Cov=ones(Q)*corr_coeff+diag(ones(Q,1)*(1-corr_coeff)); % required covariance matrix
freq=randi(20,Q,1)+10; % random frequencies between 10Hz to 30Hz
phases = 2*pi*rand(Q,1); % random phases
t=10*pi/T:10*pi/T:10*pi;
Signals=sqrt(2)*cos(2*pi*freq*t+phases); % the basic signals
if corr_coeff < 1
    A=chol(Cov)'; % Cholesky Decomposition
    S=A*Signals;
else % Coherent Sources
    S=repmat(Signals(1,:),Q,1);
end

