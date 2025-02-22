clc; clear; close all;


load('Param_R8_sigma8_K64_T600_v0.1_var.mat');

%% Scenario configurations
[Config.length,Config.width,Config.freq,Config.time] = size(X4DT);
Config.SamplingRatio = 0.1;
Config.pattern = 2;
Config.CPrank = 5;
Config.AlgInit = 0;


outputConfig = SamplingPattern(Config);



%% Down-sample
Ycompact = cell(1,Config.time);
Wtens = outputConfig.tensor;
for tt = 1:Config.time 
    Wmatt = squeeze(Wtens(:,:,tt));
    idxt = find(Wmatt);
    Xground_truth_t = squeeze(X4DT(:,:,:,tt));
    Xground_truth_t_mode_3 = tens2mat(Xground_truth_t,[],3);
    Ycompact{tt} = Xground_truth_t_mode_3(idxt,:);

end


%% IncreDSC
% algorithm configurations
IncreDSCInput.length = Config.length;
IncreDSCInput.width = Config.width;
IncreDSCInput.data = Ycompact;
IncreDSCInput.SamplingTensor = outputConfig.tensor;

IncreDSCConfig.ForgettingFactor = 0.9;
IncreDSCConfig.CPrank = Config.CPrank;
IncreDSCConfig.NMF_init = 1;



X4DEst.IncreDSC = IncreDSC_LT(IncreDSCInput,IncreDSCConfig);


%% Calculating Reconstruction Error
NMSE.IncreDSC = zeros(1,Config.time);


for tt = 1:Config.time

    
    NMSE.IncreDSC(tt) = frob(X4DEst.IncreDSC(:,:,:,tt) - X4DT(:,:,:,tt))^2/(frob(X4DT(:,:,:,tt))^2);
    
    
end

