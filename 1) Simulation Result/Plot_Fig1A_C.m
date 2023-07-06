%% Thrteshold values are dependent upon the number of sources
clear
close all
clc

%% Parameters
dirr1 = 'E:\AGiri_FratioMethod_MatlabCode\1) Simulation Result\AllDataFiles';
THRESH = 1:0.01:1.5;
MaxnumSources = 5;

%%
Anat = 1;
ModelingError = 'X';
SNRs =-8:4:8;
idxCorr = 0;
Corr = 0.5;
for Correlation = Corr*10
    idxCorr = idxCorr+1;
    idxSNR =0;
    for SNR = SNRs
        idxSNR = idxSNR+1;
        if SNR<0
            load(fullfile(dirr1, ['DOFMError' num2str(ModelingError) 'SNR_' num2str(abs(SNR)) 'CORR0' num2str(Correlation) '_Subj0' num2str(Anat) '.mat']));
        else
            load(fullfile(dirr1, ['DOFMError' num2str(ModelingError) 'SNR' num2str(abs(SNR)) '_CORR0' num2str(Correlation) '_Subj0' num2str(Anat) '.mat']));
        end
        [OptimumThresh,BestAccuracy,Accuracy] = AccuracyandThreshold(MaxnumSources,THRESH,MC_repetitions,Fratio_RAP,Fratio_AP);
        ThresholdRAP{1,idxCorr}(idxSNR,:) = OptimumThresh(:,1)';
        ThresholdRAP{1,idxCorr}(:,1) = ThresholdRAP{1,idxCorr}(:,2);
        ThresholdAP{1,idxCorr}(idxSNR,:) = OptimumThresh(:,2)';
        ThresholdAP{1,idxCorr}(:,1) = ThresholdAP{1,idxCorr}(:,2);
        Accuracyfull{idxSNR,idxCorr} = Accuracy;
        %% Accuracy
        AccuracyRAP{1,idxCorr}(idxSNR,:) = BestAccuracy(:,1)';
        AccuracyAP{1,idxCorr}(idxSNR,:) = BestAccuracy(:,2)';
        idxSNR
        %% Thrteshold values are dependent upon the number of sources
        figure
        for i=1:MaxnumSources+1
            yy=[];
            xx = THRESH;
            for ax=1:length(THRESH)
                yy = [yy Accuracyfull{idxSNR,idxCorr}{i,ax}(1,3)];
            end
            zz = (i-1)*ones(1,length(THRESH));
            plot(xx,yy,'LineWidth',2)
            hold on
        end
        grid on
        xlabel('F-Ratio Threshold Value','FontSize', 15)
        ylabel('Number of Sources Estimation Accuracy (%)','FontSize', 15)
        set(gca,'fontweight','bold','FontSize', 15)
        legend('Q=0','Q=1','Q=2','Q=3','Q=4','Q=5','FontSize', 15)
       title(['Subj ' num2str(Anat) ', SNR=', num2str(SNR) ', Corr=' num2str(Corr) ', and MError=', num2str(ModelingError) ],'FontSize', 10)
        hold off
        view([0 90])
        x0=100;
        y0=10;
        width=500;
        height=350;
        set(gcf,'position',[x0,y0,width,height])
    end
end

