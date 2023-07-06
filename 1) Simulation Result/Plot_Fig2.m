clear
close all
clc

%% Parameters
SNR = -10:2:10; FinerSNR = -10:0.1:10;
MaxNSources = 6;

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

figure
surf(X,Y,V(1:MaxNSources,:))
title('Original Sampling');
xlabel('SNR (dB)')
ylabel('Number of Sources')
zlabel('Optimal Threshold')
view([0 90])

Vq = interp2(X,Y,V(1:MaxNSources,:),Xq,Yq,'cubic'); % cubic interpolation

figure
surf(Xq,Yq,Vq);
title('Cubic Interpolation Over Finer Grid');
xlabel('SNR (dB)')
ylabel('Number of Sources')
zlabel('Optimal Threshold')
view([0 90])