function [X4DHat,outputVar] = IncreDSC(Input,AlgConfig)
%% Proposed Incremental dynamic spectrum cartography algorithm
    %% loading scenario configurations
    if ~isfield(Input, 'data')
        error('Input size/data no detected');
    end

    Yall = Input.data; % with size |Omega_t| * |K| * |T|;
    K = size(Yall{1},2);
    T = length(Yall);

    
    if ~isfield(Input,'SamplingTensor')
        error('Sampling tensor not detected.');
    end
    Wtens = Input.SamplingTensor;
    [I,J,~,~] = size(Wtens);

    % SampleIndexAll = Input.SampleIndex;

    %% loading algorithm configurations
    
    if ~isfield(AlgConfig,'CPrank')
        warning('CPrank undefined, setting as 5 by default.');
        Fr = 5;
    else

        Fr = AlgConfig.CPrank;
    end

    if ~isfield(Input,'ForgettingFactor')
        
        lambda = 0.9;

    else
        lambda = AlgConfig.ForgettingFactor;
    end

    if lambda<0 || lambda > 1
        warning('Forgetting factor must be non-negative and no larger than 1, reset as 0.9 automatically.')
        lambda = 0.9;
    end

    

    if ~isfield(AlgConfig, 'EmitterNumber')
        %% model-order selection via SVD
        Ymat = Yall{1};
        singvalue = svd(Ymat);
        %% checking number of emitters via SVD
        for tau = 1:length(singvalue)
            toptauvals = sum(singvalue(1:tau));
            occ_ratio = toptauvals/sum(singvalue);
            if occ_ratio > 0.9999
                Rest = tau;
                break;
            end
        end
        
    else
        Rest = AlgConfig.EmitterNumber;
    end
     
    %% initializing CP and cache factors
    for rr = 1:Rest
        A{rr} = randn(I,Fr);
        B{rr} = randn(J,Fr);
        c{rr} = randn(1,Fr);
        for ii = 1:I
            Phi{rr}{ii} = 0;
            phi{rr}{ii} = 0;
        end
        for jj = 1:J
            Psi{rr}{jj} = 0;
            psi{rr}{jj} = 0;
        end

    end

    SpectraAll = cell(1,T);

    if ~isfield(AlgConfig,'NMF_init')
        NMF_init = 1;
    else
        NMF_init = AlgConfig.NMF_init; % 1. Unraveling emitters' parameters via NMF; 0. SMF-LL1
    end
    
    if ~isfield(AlgConfig,'HALS_FineTune')
        HALS_FineTune = 1;
    else
        HALS_FineTune = AlgConfig.HALS_FineTune;
    end
    
    if ~isfield(AlgConfig,'HALS_FineTuneIter')
        HALS_FineTuneIter = 50;
    else
        HALS_FineTuneIter = AlgConfig.HALS_FineTuneIter;
    end

    X4DHat = zeros(I,J,K,T);
    
    for tt = 1:T
        Yt = Yall{tt};
        Wmatt = squeeze(Wtens(:,:,tt));
        SampleIndext = find(Wmatt(:));

        if NMF_init 

            SelectInd = SPA(Yt,Rest);
            Gt = Yt(:,SelectInd);
            %column-wise normalization 
            Gt = Gt / diag(max(Gt));
            Pt = (Gt\Yt)';    

           
    
            if HALS_FineTune % fine-tunning via HALS
    
                [Gt, Pt] = NMF_HALS(Yt,Rest,HALS_FineTuneIter,Gt,Pt);
            
            end

        else

            [Gt,Pt] = SMF_LL1(I,J,Yt,Fr,Rest,SampleIndext);

        end

        SpectraAll{tt} = Pt;


        %% Unifying ambiguities
        if tt > 1
            Pt_historical = SpectraAll{tt-1};
            [permute,scale,~,Pt] = Unifying_Ambiguities(Pt_historical,Pt);
            SpectraAll{tt} =  Pt;
            Gt = Gt*permute/scale;

            SpectraAll{tt} =  Pt;
  
        end
        
        Xhatt = zeros(I,J,K);
        for rr = 1:Rest

            Gmat{rr} = zeros(I,J);
            Gmat{rr}(SampleIndext) = Gt(:,rr);

            [A{rr},B{rr},c{rr},phi{rr},Phi{rr},psi{rr},Psi{rr}] = IncreDSC_Update(A{rr},B{rr},c{rr},phi{rr},Phi{rr},psi{rr},Psi{rr},Gmat{rr},lambda);
    
            Srt = A{rr} * diag(c{rr}) * B{rr}';
            Xhatt = Xhatt + outprod(Srt, Pt(:,rr));
        
        end

        X4DHat(:,:,:,tt) = Xhatt;
       
    end
    outputVar.FactorA = A;
    outputVar.FactorB = B;
    outputVar.FactorC = c;
    outputVar.Spectra = SpectraAll;


end