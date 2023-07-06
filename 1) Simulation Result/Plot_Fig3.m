clear all
close all
clc

%%
ModelingError = 'X';
SNR = 0;
Corr = 0.1;
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

%% Bar plot
figure
v = 0:5;
v = v(ones(1,4),:);
BAR_DATA = [y1' y2' y3' y4'];
ERR_DATA = [z1' z2' z3' z4'];
errhigh = ERR_DATA;
errlow = ERR_DATA;
bcv = bar(v',BAR_DATA)
hold on
xtips1 = bcv(1).XEndPoints;
xtips2 = bcv(2).XEndPoints;
xtips3 = bcv(3).XEndPoints;
xtips4 = bcv(4).XEndPoints;
av = [xtips1; xtips2; xtips3; xtips4];
er = errorbar(av',BAR_DATA,errlow,errhigh,'linestyle','none');  
set(gca,'xticklabel',{'0','1','2','3','4','5'})
legend({'Trained Anatomy #Subj01','Tested Anatomy #Subj02','Tested Anatomy #Subj03','Tested Anatomy #Subj04'},'Location','northwest')
xlabel('True Number of Sources','FontSize',14)
ylabel('Estimated Number of Sources','FontSize',14)
graphTitle = " # SNR = " + SNR + " dB, Corr is " + corr;
title(graphTitle);
grid on