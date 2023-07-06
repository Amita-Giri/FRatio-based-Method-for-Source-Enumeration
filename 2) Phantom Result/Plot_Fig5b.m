clear
close all
clc

%% Amplitude
temp = 1; % =0 represents no source and =1,2,3,4 represents that much number of dipoles
n_sources = 2; %Number of sources can be varied from 1 to 10 or more
range = 155:195;
nmax = n_sources+1;
Repetitions = 100;
T_Samples = 350;
AP_max_iters = n_sources+6;
mode = 'AP';
amp = '200'; % change pahntom directory
Phantom_Dir='E:\AGiri_FratioMethod_MatlabCode\2) Phantom Result\LPC_Filtered_PhantData';

%% Load Gain Matrix
load(strcat(Phantom_Dir,'\headmodel_vol_meg_sphere'));
N_sensors = 306;
MEG_range = 1:N_sensors;
Gain = Gain(MEG_range,:);

%% filenames
truedipoles = load(strcat(Phantom_Dir,'\TrueDipoles\dipoles_230321_1550.mat'));
NoiseCovMat = load(strcat(Phantom_Dir,'\noisecov_full.mat'));

%% Load dipoles
Dipole = truedipoles.Dipole;
N_Dipoles = size(Dipole,2);

%% Load LPC filtered MEG Average Phantom Data along 20 trials
MEG_MAT=cell(N_Dipoles,1);
listing = dir(Phantom_Dir);
dip_num=0;
for f=1:length(listing)
    fname=listing(f).name;
    TF = contains(fname,'average');
    if TF
        dip_num=dip_num+1;
        load(strcat(Phantom_Dir,'/',fname));
        MEG_MAT{dip_num}=F;
        MEG_MATD{dip_num,1} = MEG_MAT{dip_num,1}(MEG_range,:);
    end
end
MEG_MAT = MEG_MATD; clear MEG_MATD;

%% Active source index at each repititions
load Source_Ind.mat
for i = 1:Repetitions
    Source_Ind1(1:n_sources,i) = Source_Ind(1:n_sources,i);
end
Source_Ind = Source_Ind1;

%% Random Delay
if n_sources>1
    load(['Delay' num2str(n_sources) '.mat'])
end

%% Measurement Data
for i = 1:Repetitions
    MEAS_MAT{i,1} = zeros(length(MEG_range),T_Samples);
    Noise{i,1} = zeros(length(MEG_range),length(range));
    for j = 1:n_sources
        if j==1
            MEAS_MAT{i,1} = MEAS_MAT{i,1} + MEG_MAT{Source_Ind(j,i),1}(:,1:T_Samples);
        else
            MEAS_MAT{i,1} = MEAS_MAT{i,1} +  circshift(MEG_MAT{Source_Ind(j,i),1}(:,1:T_Samples),Delay(j-1,i),2);
        end
        Noise{i,1} = Noise{i,1} + MEG_MAT{Source_Ind(j,i),1}(:,10+1:10+length(range));
        True_Loc{i,j}=Dipole(Source_Ind(j,i)).Loc;
    end
    Y{i,1} = (MEAS_MAT{i,1});
end

%% noise covariance;
C_noise = NoiseCovMat.NoiseCov(MEG_range,MEG_range); % all of the channels requested
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

%% Apply Pre-Whitenig, Check PWMEG_MAT has a normal distribution
Multfactor = iW_noise*sqrt(20)./sqrt(n_sources);

Data = Multfactor*Y{1,1}(:,51:250);
x = linspace(-50,150,200);
figure; plot(x,Data')
xlabel('Time (ms)')
ylabel('fT')

