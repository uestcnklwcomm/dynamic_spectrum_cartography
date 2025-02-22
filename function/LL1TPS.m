function X4DHat = LL1TPS(Input,AlgConfig)
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



    if ~isfield(AlgConfig,'MLR')
        Lr = 6;
    else
        Lr = AlgConfig.MLR;
    end
    
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
    
    if ~isfield(AlgConfig,'NMF_init')
        NMF_init = 1;
    else
        NMF_init = AlgConfig.NMF_init;
    end

    if ~isfield(AlgConfig,'iter')
        iter = 50;
    end
    
    for tt = 1:T
        
        Yt = Yall{tt};
        Wmatt = squeeze(Wtens(:,:,tt));
        SampleIndext = find(Wmatt(:));

        M = size(Yt,1);
        P = zeros(M,I*J);
        for mm = 1:M
            row_selection = SampleIndext(mm);
            P(mm,row_selection) = 1;
        end

        if NMF_init
            SelectInd = SPA(Yt,Rest);
            S_init = Yt(:,SelectInd);
            C_init = (S_init\Yt)';
        else
            [S_init,C_init] = NMF_HALS(Yt,Rest,50);
        end
        
        S = zeros(I*J,Rest);
        S(SampleIndext,:) = S_init;
        C = C_init;

        lambda = 1e-8;
        svalue = 1e-16;

        for ii = 1:iter
            PS = P*S;
            LS = 1/norm(PS'*PS);
            C = C - rand*LS * (lambda*C + (C*PS' - Yt') * PS);
            C(C<svalue) = svalue;
            LC = 1/norm(C'*C);
            S = S - rand*LC * (P'*(PS*C'- Yt)*C); 
            for rr = 1:Rest
                sr = reshape(S(:,rr),[],1);
                SrMat = reshape(sr,I,J);
                [Us,Ss,Vs] = svds(SrMat,Lr);
                SrMat = Us*Ss*Vs';
                SrMat(SrMat<svalue) =svalue;
                sr = reshape(SrMat,[],1);
                S(:,rr) = sr;
            end
        end        
        Somega = S(SampleIndext,:);
        
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