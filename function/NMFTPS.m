function X4DHat = NMFTPS(Input,AlgConfig)
    if  ~isfield(Input, 'data')
        error('Input data not detected.');
    end

    Yall = Input.data; % size |Omega_t| * |K| * |T|;
    K = size(Yall{1},2);
    T = length(Yall);


    if ~isfield(Input,'SamplingTensor')
        error('Sampling mask not detected.');
    end

    Wtens = Input.SamplingTensor;

    [I,J,~,~] = size(Wtens);
    gridLen = I-1;
    gridResolution = 1;%
    x_grid = 0:gridResolution:gridLen;
    y_grid = 0:gridResolution:gridLen;
    [Xmesh_grid, Ymesh_grid] = meshgrid(x_grid, y_grid);
    Xgrid = Xmesh_grid + 1i*Ymesh_grid;
    
    xgrid = real(Xgrid(:));
    ygrid = imag(Xgrid(:));
    lambda_tps = 0e-4;
    dB=0;



    if ~isfield(AlgConfig,'EmitterNumber')
        %% model-order selection via SVD
        Ymat = Yall{1};
        singvalue = svd(Ymat);
 
        %% checking number of emitters via SVD
        for tau = 1:length(singvalue)
            toptauvals = sum(singvalue(1:tau));
            occ_ratio = toptauvals/sum(singvalue);
            if occ_ratio > 0.99
                Rest = tau;
                break;
            end
        end
        
    else
        Rest = AlgConfig.EmitterNumber;
    end

    %% configuring LL1 scenario
    
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
    
    for tt = 1:T
        Yt = Yall{tt};
        Wmatt = squeeze(Wtens(:,:,tt));
        SampleIndext = find(Wmatt(:));
        SelectInd = SPA(Yt,Rest);
        Somega = Yt(:,SelectInd);
        C = (Somega\Yt)';

        if HALS_FineTune
            [Somega,C] = NMF_HALS(Yt,Rest,HALS_FineTuneIter,Somega,C);
        end

        
        
        
        %% TPS completes each SLF

        Xhatt = zeros(I,J,K);
        x_IND = real(Xgrid(SampleIndext));
        y_IND = imag(Xgrid(SampleIndext));


        for rr = 1:Rest
            sromega = Somega(:,rr);
            stps = TPS(x_IND,y_IND,sromega,xgrid,ygrid,lambda_tps,dB);
            Xhatt = Xhatt + outprod(reshape(stps,[I,J]),C(:,rr));
        end

        X4DHat(:,:,:,tt) = Xhatt;
 
        
    end
end