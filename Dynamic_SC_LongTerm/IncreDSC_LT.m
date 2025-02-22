function [X4DHat,outputVar] = IncreDSC_LT(Input,AlgConfig)
%% Proposed Long-term Incremental dynamic spectrum cartography algorithm
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
        SampleIndexAll{tt} = SampleIndext;
        %% recheck via SVD
        Ymat = Yall{tt};
        singvalue = svd(Ymat);
        for tau = 1:length(singvalue)
            toptauvals = sum(singvalue(1:tau));
            occ_ratio = toptauvals/sum(singvalue);
            if occ_ratio > 0.9999
                Rest = tau;
                break;
            end
        end

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
        SLFAll{tt} = Gt;


        %% Unifying ambiguities
        if tt > 1

            % calibrating power spectra first
            Pt_historical = SpectraAll{tt-1};
            [Pt,AmbiguitiesInfo] = Parameter_rearrange(Pt_historical,Pt);
            SpectraAll{tt} =  Pt;
            
            % calibrating observed spatial loss field
            Gt_historical = SLFAll{tt-1};
            AmbiguitiesInfo.scaling = AmbiguitiesInfo.scaling^(-1);
            Gt = Parameter_rearrange(Gt_historical,Gt,AmbiguitiesInfo);
            SLFAll{tt} = Gt;



            SampleIndex_historical = SampleIndexAll{tt-1};
            coupled_idx = intersect(SampleIndex_historical, SampleIndext);

            [~,coupled_idx_t] = ismember(coupled_idx,SampleIndext);
            [~,coupled_idx_historical] = ismember(coupled_idx,SampleIndex_historical);

            coupled_idx_t(coupled_idx_t>0);
            coupled_idx_historical(coupled_idx_historical>0);

            
            
            Gt_historical = SLFAll{tt-1};

            Gcoupled_historical = Gt_historical(coupled_idx_historical,:);
            Gcoupled_t = Gt(coupled_idx_t,:);

      
            [~,~,congruence,~] = Unifying_Ambiguities(Gcoupled_historical,Gcoupled_t);
            

            %reseting CP and auxiliary variables using a heuristic threshold

            th = 0.4;
            A = Factor_reset(A,congruence,1,th);
            B = Factor_reset(B,congruence,1,th);
            c = Factor_reset(c,congruence,1,th);
            Phi = Factor_reset(Phi,congruence,0,th);
            phi = Factor_reset(phi,congruence,0,th);
            Psi = Factor_reset(Psi,congruence,0,th);
            psi = Factor_reset(psi,congruence,0,th);

           
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