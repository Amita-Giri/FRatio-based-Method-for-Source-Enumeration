clear
close all
clc

tic
%% Parameters
ModelingError = 'Z';
Corr = 0.5;
Anatomy = 4; % [1:4]
MC_repetitions = 100; % Monte-Carlo repetitions
T = 50; % time samples
True_numS = 0:5;
Correlation = Corr*10;

%% Main Simulation part
for SNR = -10:2:10
    SNR
    modeldir = ['E:\AGiri_FratioMethod_MatlabCode\1) Simulation Result\subj0' num2str(Anatomy) 'NN'];
    for True_nsources = True_numS
        idxNSources = True_nsources+1;
        AP_max_iters = True_nsources + 6; % AP maximal number of iterations
        nmax = True_nsources+2;

        %% Load Ground Truth (GT) sampling points (15002) on the cortex and the Gain matrix
        load(fullfile(modeldir, 'headmodel_surf_openmeeg_MEG_groundtruth_constrained.mat'));
        GT_Gain = Gain; clear Gain;
        GT_Gain_std=std(GT_Gain);
        GT_Gain = GT_Gain./repmat(GT_Gain_std,size(GT_Gain,1),1);
        N_sensors = size(GT_Gain,1); % number of sensors

        %% Load Gain matrix of 1 mm posterior modeling error condition
        if ModelingError == 'X'
            load(fullfile(modeldir, 'headmodel_surf_openmeeg_MEG_regrid_x_1mm_posterior_constrained'));
        elseif ModelingError == 'Y'
            load(fullfile(modeldir, 'headmodel_surf_openmeeg_MEG_regrid_y_1mm_right_constrained'));
        else
            load(fullfile(modeldir, 'headmodel_surf_openmeeg_MEG_regrid_z_1mm_upward_constrained'));
        end
        gain_std=std(Gain);
        GAIN = Gain./repmat(gain_std,size(Gain,1),1);

        %% Active source index at each repititions
        load SI4new % random source indexes are saved,////////si = randi(ndipoles,nsources,MC_repetitions);
        si = si(1:True_nsources,1:MC_repetitions); % active source index at each repitions

        %% F-Statistics AP Implementation on AP Algorithm
        gainMatrix = GAIN;
        ngrids = size(GAIN,2);
        parfor idxMC = 1:MC_repetitions
            [Y{idxNSources,idxMC},S_True{idxNSources,idxMC},Noise{idxNSources,idxMC}] = CreateMatrix_mstd(Corr, T, True_nsources, GT_Gain, si, idxMC, N_sensors, SNR); % MEG Data
            disp(" SNR #" + SNR + " Correlation #" + Corr + " True_nsources #" + True_nsources + " MC Repititions #" + idxMC);

            %% F-Statistics
            [Fratio_RAP{idxNSources,idxMC},S_RAP_Reduced{idxNSources,idxMC}] = Fstatistics(Y{idxNSources,idxMC},ngrids,gainMatrix,AP_max_iters,T,N_sensors,nmax,'RAP');
            [Fratio_AP{idxNSources,idxMC},S_AP_Reduced{idxNSources,idxMC}] = Fstatistics(Y{idxNSources,idxMC},ngrids,gainMatrix,AP_max_iters,T,N_sensors,nmax,'AP');
        end
    end
    %     if SNR<0
    %         save(['DOFMError' num2str(ModelingError) 'SNR_' num2str(abs(SNR)) 'CORR0' num2str(Correlation) '_Subj0' num2str(Anatomy) '.mat'])
    %     else
    %         save(['DOFMError' num2str(ModelingError) 'SNR' num2str(abs(SNR)) '_CORR0' num2str(Correlation) '_Subj0' num2str(Anatomy) '.mat'])
    %     end
end
toc



