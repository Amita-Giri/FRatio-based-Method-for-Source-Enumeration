clear 
close all
clc

%% Compute complete Table for different SNR and NOS
idxSNR = 0;
for SNR = -10:2:22
    SNRTEM = SNR;
    idxSNR = idxSNR+1;
    %%
    Anat = 1;
    MaxnumSources = 6;
    dirr1 = 'E:\AGiri_FratioMethod_MatlabCode\1) Simulation Result\AllDataFiles';
    THRESH = 1:0.01:10;
    idxNOS = 0;
    for NOS = 0:MaxnumSources
        idxNOS = idxNOS+1;
        Txx = []; Tyy = []; Tzz = [];
        idxMError = 0;
        for ModelingError = ['X' 'Y' 'Z']
            idxMError = idxMError+1;
            idxCorr = 0;
            for Correlation = [1 5 9]
                idxCorr = idxCorr+1;
                if SNR<0
                    load(fullfile(dirr1, ['DOFMError' num2str(ModelingError) 'SNR_' num2str(abs(SNR)) 'CORR0' num2str(Correlation) '_Subj0' num2str(Anat) '.mat']));
                else
                    load(fullfile(dirr1, ['DOFMError' num2str(ModelingError) 'SNR' num2str(abs(SNR)) '_CORR0' num2str(Correlation) '_Subj0' num2str(Anat) '.mat']));
                end
                [Accuracy] = Thresholdfunction(MaxnumSources,THRESH,MC_repetitions,Fratio_AP);
                Accuracyfull{idxSNR,idxMError,idxCorr} = Accuracy;
                if idxMError == 1
                    Txx = [Txx Accuracyfull{idxSNR,idxMError,idxCorr}(idxNOS,:)'];
                elseif idxMError == 2
                    Tyy = [Tyy Accuracyfull{idxSNR,idxMError,idxCorr}(idxNOS,:)'];
                else
                    Tzz = [Tzz Accuracyfull{idxSNR,idxMError,idxCorr}(idxNOS,:)'];
                end
            end
        end
        %% Compute average value
        AVG = (sum(Txx,2)+sum(Tyy,2)+sum(Tzz,2))./9;
        [Maxval,ind] = max(AVG);
        OptimalThreshold(idxSNR,idxNOS) = THRESH(1,ind);
        Acc(idxSNR,NOS+1) = Maxval;
        SATVE(idxSNR,idxNOS) = THRESH(1,ind);
    end
    idxSNR
    clear Accuracy_RAP
    clear Accuracy_AP
end

%% save workspace variables
% save OptimalThreshold_Subj01.mat OptimalThreshold


