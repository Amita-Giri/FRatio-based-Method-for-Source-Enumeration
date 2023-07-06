function [Y,S,Noise] = CreateMatrix_mstd(corr, T, True_nsources, gainMatrix, si, idxMC, N_sensors, snr)
if True_nsources== 0
    S = zeros(True_nsources,T);
else
    S = gen_correlated_sources(corr,T,True_nsources);
end
M = gainMatrix(:,si(:,idxMC)') * S;

noise_var = 1; meann = 0;
Noise = meann+randn(N_sensors,T).*sqrt(noise_var);
MEG_energy = (10^(snr/10))*noise_var;

tracemsmstran = MEG_energy*(N_sensors*(T));
if True_nsources== 0
    scale = 0;
else
    scale = sqrt(tracemsmstran/trace(M*M'));
end

Ms = M*scale;
S = S*scale;

Y = Ms + Noise;
end

