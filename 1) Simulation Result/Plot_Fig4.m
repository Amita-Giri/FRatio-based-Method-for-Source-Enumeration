clear all
close all
clc

%% Parameters to set
SNR = -4;
Corr = 0.5;

%%
ModelingError = 'X';
Correlation = Corr*10;
MaxnumSources = 5;
dirr1 = 'E:\AGiri_FratioMethod_MatlabCode\1) Simulation Result\AllDataFiles';

%% Assessing the performance
load 'OptimalThreshold_Subj01.mat'
SNRs = -10:2:22;

idxSNR = find(SNRs == SNR);
idxAantomy = 0;
y1=[]; y2=[]; y3=[]; y4 = []; z1=[]; z2=[]; z3=[]; z4 = [];
for Anatomy = [1 2 3 4]
    idxAantomy = idxAantomy+1;
    if SNR<0
        load(fullfile(dirr1, ['DOFMError' num2str(ModelingError) 'SNR_' num2str(abs(SNR)) 'CORR0' num2str(Correlation) '_Subj0' num2str(Anatomy) '.mat']));
    else
        load(fullfile(dirr1, ['DOFMError' num2str(ModelingError) 'SNR' num2str(abs(SNR)) '_CORR0' num2str(Correlation) '_Subj0' num2str(Anatomy) '.mat']));
    end

    %% AP
    for True_nsources=0:MaxnumSources
        for idxMC = 1:MC_repetitions
            for k= 1:min(True_nsources+2,6)
                %% F-Ratio
                if Fratio_AP{True_nsources+1,idxMC}(k,1) < OptimalThreshold(idxSNR,k)
                    Est_no_sources_AP{idxSNR,idxAantomy}(True_nsources+1,idxMC)=k-1;
                    break
                else
                    Est_no_sources_AP{idxSNR,idxAantomy}(True_nsources+1,idxMC)=NaN;
                end
            end
        end
        if Anatomy==4
            y1 = [y1 nanmean(Est_no_sources_AP{idxSNR,idxAantomy}(True_nsources+1,:))];
            z1 = [z1 nanstd(Est_no_sources_AP{idxSNR,idxAantomy}(True_nsources+1,:))];
        elseif Anatomy==5
            y2 = [y2 nanmean(Est_no_sources_AP{idxSNR,idxAantomy}(True_nsources+1,:))];
            z2 = [z2 nanstd(Est_no_sources_AP{idxSNR,idxAantomy}(True_nsources+1,:))];
        elseif Anatomy==6
            y3 = [y3 nanmean(Est_no_sources_AP{idxSNR,idxAantomy}(True_nsources+1,:))];
            z3 = [z3 nanstd(Est_no_sources_AP{idxSNR,idxAantomy}(True_nsources+1,:))];
        elseif Anatomy==7
            y4 = [y4 nanmean(Est_no_sources_AP{idxSNR,idxAantomy}(True_nsources+1,:))];
            z4 = [z4 nanstd(Est_no_sources_AP{idxSNR,idxAantomy}(True_nsources+1,:))];
        end
    end
end
idxSNR

%% AIC and MDL method
for True_nsources=0:MaxnumSources
    for idxMC = 1:MC_repetitions
        [True_nsources idxMC]
        %Calculate eigens:
        Y_White=ABY{True_nsources+1,idxMC};
        [p,n] = size(Y_White);
        [Msnr_EigVec_k,Msnr_EigVal_k] = eigs((1/n)*(Y_White*Y_White'),p);
        Msnr_EigVal = diag(Msnr_EigVal_k); %abs is added because of very small negative values
        Msnr_EigVec = Msnr_EigVec_k;
        %Log Predictions
        [AIC_P(True_nsources+1,idxMC), MDL_P(True_nsources+1,idxMC),aic_metric,mdl_metric] = AIC_MDL(Msnr_EigVal, p, n, nmax+1);
    end
end

ax = 0:MaxnumSources;
ax = ax';
%F-Ratio
NSDR1=(y4-z4)';
PSDR1=(y4+z4)';
figure(1)
plot(ax,(NSDR1),'g'); 
hold on;
plot(ax,(PSDR1),'g');
patch([ax; flipud(ax)],[NSDR1; flipud(PSDR1)], 'g', 'FaceAlpha',0.2, 'EdgeColor','none');
plot(ax, y4','color','g','LineWidth',1.5)
hold on
% AIC
NSDR1=(mean(AIC_P')-std(AIC_P'))';
PSDR1=(mean(AIC_P')+std(AIC_P'))';
figure(1)
plot(ax,(NSDR1),'b'); 
hold on;
plot(ax,(PSDR1),'b');
patch([ax; flipud(ax)],[NSDR1; flipud(PSDR1)], 'b', 'FaceAlpha',0.2, 'EdgeColor','none');
plot(ax, mean(AIC_P')','color','b','LineWidth',1.5)
hold on
% MDL
NSDR1=(mean(MDL_P')-std(MDL_P'))';
PSDR1=(mean(MDL_P')+std(MDL_P'))';
figure(1)
plot(ax,(NSDR1),'r'); 
hold on;
plot(ax,(PSDR1),'r');
patch([ax; flipud(ax)],[NSDR1; flipud(PSDR1)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
plot(ax, mean(MDL_P')','color','r','LineWidth',1.5)
hold off
grid on
xlabel('True Number of Sources','fontweight','bold','FontSize',14)
ylabel('Estimated Number of Sources','fontweight','bold','FontSize',14)
title(['(a) SNR = ', num2str(SNR), ' dB, \rho = ', num2str(corr)]);
xticks([0 1 2 3 4 5])
legend({'','','','F-Ratio','','','','AIC','','','','MDL'},'Location','northwest')

%%
figure
plot(0:5,y4,'LineWidth', 2)
hold on;
plot(0:5,mean(AIC_P'),'LineWidth', 2)
hold on;
plot(0:5,mean(MDL_P'),'LineWidth', 2)
% hold on;
% plot(0:5,mean(Diff_P'),'LineWidth', 2)
hold off
legend({'F-Ratio','AIC','MDL'},'Location','northwest')
grid on
xlabel('True Number of Sources','fontweight','bold','FontSize',14)
ylabel('Estimated Number of Sources','fontweight','bold','FontSize',14)



