function X4DHat = LowCPRankCompletion(Input,AlgConfig)

    if  ~isfield(Input, 'data')
        error('Input data not detected.');
    end

    Yall = Input.data; % size |Omega_t| * |K| * |T|;
    K = size(Yall{1},2);
    T = length(Yall);

    

    if ~isfield(AlgConfig,'mode')
        Batchmode = 1;
    else
        Batchmode = AlgConfig.mode;
    end
    
    %% excluding all-zero frequency bins
    
    if ~isfield(Input,'Batchsize') && Batchmode
        error('Batchsize Undefined');
    end

    Batchsize = Input.Batchsize;
    BatchActivate = Batchsize:1:T;
   
    
    Kallzero = [];

    for kk = 1:K
        yk = Yall{1}(:,kk);
        if isempty(yk)
            Kallzero = [Kallzero,k];
        end
    end
    
   
    
    for tt = 1:T
        Yt = Yall{tt};
        Yt(:,Kallzero) = [];
        Yall{tt} = Yt;
    end
    
    K = size(Yall{1},2);
    

    if ~isfield(Input,'SamplingTensor')
        error('Sampling mask not detected.');
    end

    Wtens = Input.SamplingTensor;

    [I,J,~,~] = size(Wtens);

    if ~isfield(AlgConfig,'CPrank')
        F = 15;
    else
        F = AlgConfig.CPrank;
    end

    if ~isfield(AlgConfig,'Iteration')
        iter = 10;
    else
        iter = AlgConfig.Iteration;
    end
    
   
    
    if Batchmode %% DW-CPD；
        X4DHat = zeros(I,J,K,T);
        Yfreq = cell(1,K);
       
        for tt = 1:T
            Ymatt = Yall{tt};
            Wmatt = squeeze(Wtens(:,:,tt));
            SampleIndext = find(Wmatt(:));
    
            for kk = 1:K

                Yfreqslab = zeros(I,J);
                Yfreqslab(SampleIndext) = Ymatt(:,kk);
                Yfreq{kk}(:,:,tt) = Yfreqslab;
                
                if ismember(tt,BatchActivate)
                    if tt == Batchsize
                        U{kk} = CPDOL_UpdateBatch(Yfreq{kk},F,iter,Batchsize);
                    else
                        U{kk} = CPDOL_UpdateBatch(Yfreq{kk},F,iter,Batchsize,U{kk});
                    end
                    X4DHat(:,:,kk,tt-Batchsize+1:tt) = cpdgen(U{kk});
                end
            end
        end
    
    else %% I-CPD

        X4DHat = zeros(I,J,K,T);

        if ~isfield(AlgConfig,'ForgettingFactor')
            lambda = 0.9;
        else
            lambda = AlgConfig.ForgettingFactor;
        end
        
        A = cell(1,K);
        B = cell(1,K);
        for kk = 1:K 
            A{kk} = randn(I,F);
            B{kk} = randn(J,F);
            for ii = 1:I
                P{kk}{ii} = 0;
                q{kk}{ii} = 0;
            end

            for jj = 1:J
                R{kk}{jj} = 0;
                s{kk}{jj} = 0;
            end
        end

        for tt = 1:T
            Ymatt = Yall{tt};
            Wmatt = squeeze(Wtens(:,:,tt));
            SampleIndext = find(Wmatt(:));
            for kk = 1:K
                Yk = zeros(I,J);
                Yk(SampleIndext) = Ymatt(:,kk);
                [A{kk},B{kk},P{kk},q{kk},R{kk},s{kk},c{kk}] = CPDOL_UpdateIncre(A{kk},B{kk},P{kk},q{kk},R{kk},s{kk},Yk,lambda);
                X4DHat(:,:,kk,tt) = A{kk}*diag(c{kk})*(B{kk})';
            end
        end
    
    end
end