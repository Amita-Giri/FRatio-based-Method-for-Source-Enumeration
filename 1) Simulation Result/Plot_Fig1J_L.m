clear
close all
clc

%% Parameters
dirr1 = 'E:\AGiri_FratioMethod_MatlabCode\1) Simulation Result\AllDataFiles';
THRESH = 1:0.01:1.5;
MaxnumSources = 5;

%% Thrteshold values are independent upon the Modeling error
idxMError = 0;
Corr = 0.5;
Anat = 1;
Correlation = Corr*10;
for ModelingError = ['X' 'Y' 'Z']
    idxMError = idxMError+1;
    idxSNR =0;
    for SNR = 0 %0
        idxSNR = idxSNR+1;
        if SNR<0
            load(fullfile(dirr1, ['DOFMError' num2str(ModelingError) 'SNR_' num2str(abs(SNR)) 'CORR0' num2str(Correlation) '_Subj0' num2str(Anat) '.mat']));
        else
            load(fullfile(dirr1, ['DOFMError' num2str(ModelingError) 'SNR' num2str(abs(SNR)) '_CORR0' num2str(Correlation) '_Subj0' num2str(Anat) '.mat']));
        end
        [OptimumThresh,BestAccuracy,Accuracy] = AccuracyandThreshold(MaxnumSources,THRESH,MC_repetitions,Fratio_RAP,Fratio_AP);
        ThresholdRAP{1,idxMError}(idxSNR,:) = OptimumThresh(:,1)';
        ThresholdRAP{1,idxMError}(:,1) = ThresholdRAP{1,idxMError}(:,2);
        ThresholdAP{1,idxMError}(idxSNR,:) = OptimumThresh(:,2)';
        ThresholdAP{1,idxMError}(:,1) = ThresholdAP{1,idxMError}(:,2);
        Accuracyfull{idxSNR,idxMError} = Accuracy;
        %% Accuracy
        AccuracyRAP{1,idxMError}(idxSNR,:) = BestAccuracy(:,1)';
        AccuracyAP{1,idxMError}(idxSNR,:) = BestAccuracy(:,2)';
        idxSNR
    end
end

for j=1:length(Accuracyfull{1,1}(:,1))
    figure
    for i=1:3 %Modeling error
        yy=[];
        xx = THRESH;
        for ax=1:length(THRESH)
            yy = [yy Accuracyfull{idxSNR,i}{j,ax}(1,3)];
        end
        plot(xx,yy,'LineWidth',2)
        hold on
    end
    grid on
    xlabel('F-Ratio Threshold Value','FontSize', 15)
    ylabel('Number of Sources Estimation Accuracy (%)','FontSize', 15)
    set(gca,'fontweight','bold','FontSize', 15)
    legend('Model Error in X','Model Error in Y','Model Error in Z','FontSize', 15)
    title(['Subj ' num2str(Anat) ', SNR=', num2str(SNR) ', Corr=' num2str(Corr) ', and NOS=', num2str(j-1)],'FontSize', 10)
    hold off
    view([0 90])
    x0=100;
    y0=10;
    width=500;
    height=350;
    set(gcf,'position',[x0,y0,width,height])
end