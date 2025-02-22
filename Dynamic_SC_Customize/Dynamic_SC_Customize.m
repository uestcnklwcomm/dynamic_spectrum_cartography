clc; clear; close all;
addpath(genpath('function'));




%% Scenario configurations
Config.EmitterNumber = 8;
Config.sigma = 8; %shadowing co-variance
Config.decorr = 50; %shadowing de-correlation distance
Config.speed = 0.01;
Config.width = 51;
Config.length = 51;
Config.freq = 64;
Config.time = 600;

X4DT = RadioMapGenerator(Config);

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


%% Static Approaches： LL1-TPS; NMF-TPS; Plain TPS

%% LL1-TPS
LL1Input.data = Ycompact;
LL1Input.SamplingTensor = outputConfig.tensor;
LL1Config.MLR = 5;
X4DEst.LL1TPS = LL1TPS(LL1Input,LL1Config);

%% NMF-TPS
NMFInput.data = Ycompact;
NMFInput.SamplingTensor = outputConfig.tensor;
NMFConfig.HALS_FineTune = 1;



X4DEst.NMFTPS = NMFTPS(NMFInput,NMFConfig);

%% Plain-TPS
X4DEst.PlainTPS = PlainTPS(Ycompact,outputConfig.tensor);

%% CPD-OL
CPDOLInput.data = Ycompact;
CPDOLInput.SamplingTensor = outputConfig.tensor;
CPDOLInput.Batchsize = outputConfig.Batchsize;

CPDOLConfig.mode = 0;
CPDOLConfig.CPrank = 15;
CPDOLConfig.Iteration = 10;


X4DEst.CPDOL = LowCPRankCompletion(CPDOLInput,CPDOLConfig);



% %% TOUCAN
% save('SamplingMask.mat','Wtens');
% system('python TOUCAN_main.py');
% X4DEst.TOUCAN = load('X4DHatTOUCAN.mat','X4DHatTOUCAN');
%% BatchDSC
BatchDSCInput.length = Config.length;
BatchDSCInput.width = Config.width;
BatchDSCInput.data = Ycompact;
BatchDSCInput.SamplingTensor = outputConfig.tensor;
BatchDSCInput.Batchsize = outputConfig.Batchsize;




BatchDSCConfig.CPrank = Config.CPrank;
BatchDSCConfig.Iteration = 10;
BatchDSCConfig.NMF_init = 1;



if isfield(BatchDSCConfig,'InitAB') && BatchDSCConfig.InitAB
    BatchDSCConfig.alg_row = outputConfig.alg_row;
    BatchDSCConfig.alg_col = outputConfig.alg_col;
    BatchDSCConfig.alg_mask = outputConfig.AlgebraicMask;
end


[X4DEst.BatchDSC,~] = BatchDSC(BatchDSCInput,BatchDSCConfig);


%% IncreDSC
% algorithm configurations
IncreDSCInput.length = Config.length;
IncreDSCInput.width = Config.width;
IncreDSCInput.data = Ycompact;
IncreDSCInput.SamplingTensor = outputConfig.tensor;

IncreDSCConfig.ForgettingFactor = 0.9;
IncreDSCConfig.CPrank = Config.CPrank;
IncreDSCConfig.NMF_init = 1;



[X4DEst.IncreDSC,~] = IncreDSC(IncreDSCInput,IncreDSCConfig);




%% Calculating Reconstruction Error
NMSE.LL1TPS = zeros(1,Config.time);
NMSE.NMFTPS = zeros(1,Config.time);
NMSE.PlainTPS = zeros(1,Config.time);
NMSE.CPDOL = zeros(1,Config.time);
NMSE.IncreDSC = zeros(1,Config.time);
NMSE.BatchDSC = zeros(1,Config.time);

for tt = 1:Config.time

    NMSE.LL1TPS(tt) = frob(X4DEst.LL1TPS(:,:,:,tt) - X4DT(:,:,:,tt))^2/(frob(X4DT(:,:,:,tt))^2);
    NMSE.NMFTPS(tt) = frob(X4DEst.NMFTPS(:,:,:,tt) - X4DT(:,:,:,tt))^2/(frob(X4DT(:,:,:,tt))^2);
    NMSE.PlainTPS(tt) = frob(X4DEst.PlainTPS(:,:,:,tt) - X4DT(:,:,:,tt))^2/(frob(X4DT(:,:,:,tt))^2);
    NMSE.CPDOL(tt) = frob(X4DEst.CPDOL(:,:,:,tt) - X4DT(:,:,:,tt))^2/(frob(X4DT(:,:,:,tt))^2);

    % NMSE.TOUCAN(tt) = frob(X4DEst.TOUCAN.X4DHatTOUCAN(:,:,:,tt) - X4DT(:,:,:,tt))^2/(frob(X4DT(:,:,:,tt))^2);

    NMSE.BatchDSC(tt) = frob(X4DEst.BatchDSC(:,:,:,tt) - X4DT(:,:,:,tt))^2/(frob(X4DT(:,:,:,tt))^2);
    
    NMSE.IncreDSC(tt) = frob(X4DEst.IncreDSC(:,:,:,tt) - X4DT(:,:,:,tt))^2/(frob(X4DT(:,:,:,tt))^2);
    

    
end


%% Visualization
Visualization = 1;

if Visualization
    Vis.time = 550;
    Vis.bin = 56;
    
    GT = squeeze(X4DT(:,:,Vis.bin,Vis.time)); % Ground-truth
    LL1Vis = squeeze(X4DEst.LL1TPS(:,:,Vis.bin,Vis.time));
    NMFVis = squeeze(X4DEst.NMFTPS(:,:,Vis.bin,Vis.time));
    TPSVis = squeeze(X4DEst.PlainTPS(:,:,Vis.bin,Vis.time));
    CPDVis = squeeze(X4DEst.CPDOL(:,:,Vis.bin,Vis.time));
    % TOUCANVis = squeeze(X4DEst.TOUCAN.X4DHatTOUCAN(:,:,Vis.bin,Vis.time));
    IncreVis = squeeze(X4DEst.IncreDSC(:,:,Vis.bin,Vis.time));
    BatchVis = squeeze(X4DEst.BatchDSC(:,:,Vis.bin,Vis.time));
    
    val_min = min(min(GT));
    val_max = max(max(GT));
    
    LL1Vis(LL1Vis<val_min) = val_min;
    LL1Vis(LL1Vis>val_max) = val_max;
    
    NMFVis(NMFVis<val_min) = val_min;
    NMFVis(NMFVis>val_max) = val_max;

    TPSVis(TPSVis<val_min) = val_min;
    TPS(TPSVis>val_max) = val_max;

    CPDVis(CPDVis<val_min) = val_min;
    CPDVis(CPDVis>val_max) = val_max;
    
    

    % TOUCANVis(TOUCANVis<val_min) = val_min;
    % TOUCANVis(TOUCANVis>val_max) = val_max;
    
    BatchVis(BatchVis<val_min) = val_min;
    BatchVis(BatchVis>val_max) = val_max;
    
    IncreVis(IncreVis<val_min) = val_min;
    IncreVis(IncreVis>val_max) = val_max;
    
    
    figure(1)
    
    subplot(241)
    contourf(10*log10(GT),100,'linecolor','None');
    colormap jet;
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    title('Ground truth')
    set(gca,'FontName','Times New Roman','FontSize',10,'LineWid',1);
    axes('position',[0.2,0.02,.6,.3])
    axis off
    my_handle = colorbar('east');
    my_handle.Title.String='dB';
    
    subplot(242)
    contourf(10*log10(LL1Vis),100,'linecolor','None');
    colormap jet;
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    title(['LL1-TPS,NMSE(dB) =', num2str(10*log10(NMSE.LL1TPS(Vis.time)),'%3.3f')])
    set(gca,'FontName','Times New Roman','FontSize',10,'LineWid',1);
    axes('position',[0.2,0.02,.6,.3])
    axis off
    
    
    subplot(243)
    contourf(10*log10(NMFVis),100,'linecolor','None');
    colormap jet;
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    title(['NMF-TPS,NMSE(dB) =', num2str(10*log10(NMSE.NMFTPS(Vis.time)),'%3.3f')])
    set(gca,'FontName','Times New Roman','FontSize',10,'LineWid',1);
    axes('position',[0.2,0.02,.6,.3])
    axis off
    
    subplot(244)
    contourf(10*log10(TPSVis),100,'linecolor','None');
    colormap jet;
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    title(['Plain-TPS,NMSE(dB) =', num2str(10*log10(NMSE.PlainTPS(Vis.time)),'%3.3f')])
    set(gca,'FontName','Times New Roman','FontSize',10,'LineWid',1);
    axes('position',[0.2,0.02,.6,.3])
    axis off

    subplot(245)
    contourf(10*log10(CPDVis),100,'linecolor','None');
    colormap jet;
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    title(['CPDOL,NMSE(dB) =', num2str(10*log10(NMSE.CPDOL(Vis.time)),'%3.3f')])
    set(gca,'FontName','Times New Roman','FontSize',10,'LineWid',1);
    axes('position',[0.2,0.02,.6,.3])
    axis off


    % subplot(246)
    % contourf(10*log10(TOUCANVis),100,'linecolor','None');
    % colormap jet;
    % set(gca,'xtick',[],'xticklabel',[])
    % set(gca,'ytick',[],'yticklabel',[])
    % title(['TOUCAN,NMSE(dB) =', num2str(10*log10(NMSE.TOUCAN(Vis.time)),'%3.3f')])
    % set(gca,'FontName','Times New Roman','FontSize',10,'LineWid',1);
    % axes('position',[0.2,0.02,.6,.3])
    % axis off
    % 
    
    subplot(247)
    contourf(10*log10(BatchVis),100,'linecolor','None');
    colormap jet;
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    title(['BatchDSC,NMSE(dB) =', num2str(10*log10(NMSE.BatchDSC(Vis.time)),'%3.3f')])
    set(gca,'FontName','Times New Roman','FontSize',10,'LineWid',1);
    axes('position',[0.2,0.02,.6,.3])
    axis off
    
    
    subplot(248)
    contourf(10*log10(IncreVis),100,'linecolor','None');
    colormap jet;
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    title(['IncreDSC,NMSE(dB) =', num2str(10*log10(NMSE.IncreDSC(Vis.time)),'%3.3f')])
    set(gca,'FontName','Times New Roman','FontSize',10,'LineWid',1);
    axes('position',[0.2,0.02,.6,.3])
    axis off
end

