clear
close all
clc

%% Parameters
amp = 200;
SourceNum = 5;
SNR = -10:2:22; FinerSNR = -10:0.1:22;
MC_repetitions = 100; % Monte-Carlo repetitions
MaxNSources = 7;

%% Load threshold value from subj04 simulation result
load('E:\AGiri_FratioMethod_MatlabCode\1) Simulation Result\OptimalThreshold_Subj01.mat');
V = OptimalThreshold(1:length(SNR),1:MaxNSources)';
V(1,:) = V(2,:);


%% Curve Fit
X = meshgrid(SNR); Xq = meshgrid(FinerSNR);
X = X(1:MaxNSources,:); Xq = Xq(1:MaxNSources,:);
Y = meshgrid(1:length(SNR)); Yq = meshgrid(1:length(FinerSNR));
Y = Y(:,1:MaxNSources)-1;  Yq = Yq(:,1:MaxNSources)-1;
Y = Y'; Yq = Yq';
Vq = interp2(X,Y,V(1:MaxNSources,:),Xq,Yq,'cubic'); % cubic interpolation

%% Accuracy estimation
Fratiomean = []; Fratiostd = [];
AICmean = []; AICstd = [];
MDLmean = []; MDLstd = [];
for SourceNum = 0:5
    load(['LPC_VOLsour' num2str(SourceNum) '.mat'])
    tolerance = 0.001;
    for i=1:MC_repetitions
        for j=0:nmax-1
            [~,m1] = min(abs(FinerSNR-Est_SNR(i,1)));
            Thresh = Vq(j+1,m1);
            if Fratiosupek_Combined{i,1}(j+1,1) < Thresh
                Estsources(i,1) = j;
                break
            else
                Estsources(i,1) = NaN;
            end
        end
    end

    %% AIC and MDL
    %Calculate eigens:
    parfor i = 1:length(Signal_Data)
        Y_White=Signal_Data{i,1};
        [p,n] = size(Y_White);
        [Msnr_EigVec_k,Msnr_EigVal_k] = eigs((1/n)*(Y_White*Y_White'),p);
        Msnr_EigVal = diag(Msnr_EigVal_k); %abs is added because of very small negative values
        Msnr_EigVec = Msnr_EigVec_k;
        %Log Predictions
        [AIC_P(i,1), MDL_P(i,1),aic_metric,mdl_metric] = AIC_MDL(Msnr_EigVal, p, n, 300);
    end

    Fratiomean = [Fratiomean mean(Estsources)]; Fratiostd = [Fratiostd std(Estsources)];
    AICmean = [AICmean mean(AIC_P)]; AICstd = [AICstd std(AIC_P)];
    MDLmean = [MDLmean mean(MDL_P)]; MDLstd = [MDLstd std(MDL_P)];
end


x = linspace(0,5,6);
y =   Fratiomean;
err = Fratiostd;
errorbar(x,y,err,'-s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
grid on
xlabel('True Number of Sources')
ylabel('Estimated Number of Sources')
