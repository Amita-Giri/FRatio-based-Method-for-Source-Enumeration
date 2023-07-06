%% lpc implemetataion in this
clear
close all
clc

%% Parameters
Signal_Range = 300:330;
Noise_Range = 10:200; % Noise covariance matrix calculated in Brainstorm is from -190 to 1ms. Therefore Noise_Range is defined 10:200
True_Num_Sources = 2;
AP_max_iters = True_Num_Sources+6;
mode = 'AP';
nmax = True_Num_Sources+2;

%% Amplitude
load('data_1_average_230508_1003.mat');
MEG_Range = 1:306;
N_sensors = length(MEG_Range);
Data = F;
Data = Data(MEG_Range,:);

%% Load Gain Matrix
load('headmodel_vol_os_meg.mat');
Gain = Gain(MEG_Range,:);
N_Dipoles = length(Gain)/3;

%% noise covariance;
load('noisecov_full.mat');
C_noise = NoiseCov(MEG_Range,MEG_Range);
a1 = C_noise;
C_noise = (C_noise + C_noise')/2;
[Un,Sn2] = svd(C_noise,'econ');
Sn = sqrt(diag(Sn2)); % singular values
tol = length(Sn) * eps(single(Sn(1))); % single precision tolerance
Rank_Noise = sum(Sn > tol);
Un = Un(:,1:Rank_Noise);
Sn = Sn(1:Rank_Noise);
% now rebuild the noise covariance matrix with just the non-zero components
C_noise = Un*diag(Sn.^2)*Un'; % possibly deficient matrix now
%regularization method: none (this is inverse noise covariance (thus, whitener matrix))
iW_noise = Un*diag(1./(Sn))*Un'; % inverse whitener
% % so Ww * Ww' = Cw, and iWw*iWw' = pinv(Cw)

%% QQ Plot and whitened Data
Multfactor = iW_noise*sqrt(166); % 166 files are there in data 1
PWMEG_MAT = Multfactor*Data;
sspw = sum(PWMEG_MAT(MEG_Range,Noise_Range).^2);
stdnoise = (std((PWMEG_MAT(MEG_Range,Noise_Range))'))';
temp = PWMEG_MAT(MEG_Range,Noise_Range);
pd = makedist('Normal');
figure; qqplot(temp(:),pd); grid on; axis on;

%% Estimate SNR
Signal_Data = PWMEG_MAT(MEG_Range,Signal_Range);
Noise_Data = PWMEG_MAT(MEG_Range,50:50+length(Signal_Range)-1);
Est_SNR = snr(Signal_Data,Noise_Data)

%% Load threshold value from subj04 simulation result
load('E:\AGiri_FratioMethod_MatlabCode\1) Simulation Result\OptimalThreshold_Subj01.mat');
V = OptimalThreshold';
SNR = -10:2:22; FinerSNR = -10:0.1:22;
MaxNSources = 7;

%% Curve Fit
X = meshgrid(SNR); Xq = meshgrid(FinerSNR);
X = X(1:MaxNSources,:); Xq = Xq(1:MaxNSources,:);
Y = meshgrid(1:length(SNR)); Yq = meshgrid(1:length(FinerSNR));
Y = Y(:,1:MaxNSources);  Yq = Yq(:,1:MaxNSources);
Y = Y'; Yq = Yq';

Vq = interp2(X,Y,V,Xq,Yq,'cubic'); % cubic interpolation
Fvalues_at_EstSNR = interp2(X,Y,V,Est_SNR,1:5,'cubic');

%% Source Localization
T = length(Signal_Range);
Y_White = Multfactor*Data(MEG_Range,Signal_Range);
Gain1 = Multfactor*Gain;

%% AIC and MDL method
[p,n] = size(Y_White);
[Msnr_EigVec_k,Msnr_EigVal_k] = eigs((1/n)*(Y_White*Y_White'),p);
Msnr_EigVal = diag(Msnr_EigVal_k); %abs is added because of very small negative values
Msnr_EigVec = Msnr_EigVec_k;
%Log Predictions
[AIC_P, MDL_P,aic_metric,mdl_metric] = AIC_MDL(Msnr_EigVal, p, n, 306);

%% F-Ratio based Method
for nsources = 0:nmax
    if nsources == 0
        [S_Reduced{nsources+1,1},~,q_Red,~] = AP_Free_Orient_Amita(Y_White, Gain1, nsources, AP_max_iters,mode,GridLoc);
    else
        S_Reduced{nsources+1,1} = S_Full{nsources,1};  q_Red = q_Fu;
    end
    [S_Full{nsources+1,1},~,q_Fu,~] = AP_Free_Orient_Amita(Y_White, Gain1, nsources+1, AP_max_iters,mode,GridLoc);

    gain_idxR = []; gain_idxF = [];
    q_Reduced = []; q_Full =[];
    Y_Reduced = zeros(N_sensors,T); Y_Full = zeros(N_sensors,T);
    for i = 1:length(S_Reduced{nsources+1,1})
        gain_idxR = [gain_idxR (S_Reduced{nsources+1,1}(1,i)-1)*3+1:S_Reduced{nsources+1,1}(1,i)*3];
        q_Reduced = blkdiag(q_Reduced, (q_Red{1,i}/norm(q_Red{1,i})));
    end
    qq_Reduced{nsources+1,1} = q_Reduced;
    A_Reduced = Gain1(:,gain_idxR)*q_Reduced; S_Est_Reduced{nsources+1,1} = (A_Reduced'*A_Reduced)\(A_Reduced'*Y_White);
    Y_Reduced = A_Reduced*S_Est_Reduced{nsources+1,1};
    for i = 1:length(S_Full{nsources+1,1})
        gain_idxF= [gain_idxF (S_Full{nsources+1,1}(1,i)-1)*3+1:S_Full{nsources+1,1}(1,i)*3];
        q_Full = blkdiag(q_Full, (q_Fu{1,i}/norm(q_Fu{1,i})));
    end
    qq_Full{nsources+1,1} = q_Full;
    A_Full = Gain1(:,gain_idxF)*q_Full; S_Est_Full{nsources+1,1} = (A_Full'*A_Full)\(A_Full'*Y_White);
    Y_Full = A_Full*S_Est_Full{nsources+1,1};
    Grand_Reduced = sum((Y_White-Y_Reduced).^2);
    Grand_Full = sum((Y_White-Y_Full).^2);
    SSE_Reduced{nsources+1,1} = Grand_Reduced;
    SSE_Full{nsources+1,1} = Grand_Full;

    %% F-Ratio ANlysis
    req = 4+T;
    DOF_Reduced = N_sensors*T-req*(nsources);
    DOF_Full = N_sensors*T-req*(nsources+1);
    Fratiosupek_Combined(nsources+1,1) = ((sum(SSE_Reduced{nsources+1,1}))/(DOF_Reduced))/((sum(SSE_Full{nsources+1,1}))/(DOF_Full));
    nsources
end

%% Decision Criteria
for nsources = 0:nmax
    if Fratiosupek_Combined(nsources+1,1) < Fvalues_at_EstSNR(nsources+1,1)
        Est_NumSources = nsources;
        break;
    end
end

[AIC_P MDL_P Est_NumSources]