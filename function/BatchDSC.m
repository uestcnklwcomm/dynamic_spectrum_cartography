function [X4DHat,outputVar] = BatchDSC(Input,AlgConfig)
    %% loading scenario configurations
    % warning off;
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
    

    if ~isfield(Input,'Batchsize')
        error('Batchsize undefined.')
    end

    Batchsize = Input.Batchsize;
    BatchDSCActivated = Batchsize:1:T;

    
    
    % SampleIndexAll = Input.SampleIndex;

    if ~isfield(AlgConfig,'CPrank')
        error('Hyperparameter CP-rank undefined.');
    end

    
    Fr = AlgConfig.CPrank;

    if ~isfield(AlgConfig,'Iteration')
        ALS_iter = 10;
    else
        ALS_iter = AlgConfig.Iteration;
    end

    if ~isfield(AlgConfig, 'EmitterNumber')
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
    
    A = cell(1,Rest);
    B = cell(1,Rest);
    C = cell(1,Rest);


  
    for rr = 1:Rest
        A{rr} = randn(I,Fr);
        B{rr} = randn(J,Fr);
        C{rr} = randn(Batchsize,Fr);
    end
        

    if ~isfield(AlgConfig,'InitAB')
        AlgInitAB = 0;
    else
        AlgInitAB = AlgConfig.InitAB;
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
     

    SpectraAll = cell(1,T);

    for tt = 1:T
        Yt = Yall{tt};
        Wmatt = squeeze(Wtens(:,:,tt));
        SampleIndext = find(Wmatt(:));

        if AlgConfig.NMF_init % 1. Unraveling emitters' parameters via NMF; 0. SMF-LL1

            SelectInd = SPA(Yt,Rest);
            Gt = Yt(:,SelectInd);
            %column-wise normalization 
            Gt = Gt / diag(max(Gt));
            Pt = (Gt\Yt)';    
    
            if HALS_FineTune % fine-tunning via HALS
    
                [Gt, Pt] = NMF_HALS(Yt,Rest,HALS_FineTuneIter,Gt,Pt);
            
            end

        else

            [Gt,Pt] = SMF_LL1(I,J,Yt,Fr+1,Rest,SampleIndext);

        end

        SpectraAll{tt} = Pt;


        %% Unifying ambiguities
        if tt > 1
            Pt_historical = SpectraAll{tt-1};
            % [~,permute,scale,Pt] = cpderr(Pt_historical,Pt);
            [permute,scale,~,Pt] = Unifying_Ambiguities(Pt_historical,Pt);
            
        
            Gt = Gt*permute/scale;


            SpectraAll{tt} =  Pt;
  
        end
        
        Xhatt = zeros(I,J,K);
        for rr = 1:Rest

            Gmatrt = zeros(I,J);
            Gmatrt(SampleIndext) = Gt(:,rr);
            Gtens{rr}(:,:,tt) = Gmatrt;

           
            if ismember(tt,BatchDSCActivated) 

                Wbatch = Wtens(:,:,tt-Batchsize+1:tt);
                Gbatch = Gtens{rr}(:,:,tt-Batchsize+1:tt);
                
                %% checking whether algebraic initialization is available

                if tt == Batchsize && AlgInitAB
                    alg_row = AlgConfig.alg_row;
                    alg_col = AlgConfig.alg_col;
                    alg_G = Gbatch(alg_row,alg_col,:);

                    num_coupled_tensor = nchoosek(numel(alg_row),2) * nchoosek(numel(alg_col),2);

                    if numel(alg_row) < 1 || numel(alg_col) < 1 || num_coupled_tensor < nchoosek(Fr,2)
                        error('Insufficient amount of subtensors');
                    end


                    U_alg = Alg_CPD(alg_G,Fr);

                    C{rr} = U_alg{end};
                    C{rr} = C{rr}/diag(max(C{rr}));
         

                    alg_mask = AlgConfig.alg_mask;
                    alg_SampleIndex = find(alg_mask(:));

                    Gbatch3 = tens2mat(Gbatch,[],3);
                    Gomega = Gbatch3(alg_SampleIndex,:);
                    
                    krAB = Gomega/(C{rr}');
                  

                    for ff = 1:Fr
                        sparseABmat = zeros(I,J);
                        sparseABmat(alg_SampleIndex) = krAB(:,ff);
                        
                        [A{rr}(:,ff),B{rr}(:,ff)] = complete_rank1_matrix(sparseABmat);
        
                    end

                end


                [A{rr},B{rr},C{rr}] = BatchDSC_Update(Gbatch,Wbatch,Fr,ALS_iter,A{rr},B{rr},C{rr});
                % [A{rr},B{rr},C{rr}] = BatchCPD_BCD_Update_R(Gbatch,Wbatch,Fr,ALS_iter,A{rr},B{rr});

                Sestr = A{rr}*diag(C{rr}(end,:))*(B{rr})';
                Xhatt =  Xhatt + outprod(Sestr,Pt(:,rr));

            end
   
        end



        X4DHat(:,:,:,tt) = Xhatt;
       
    end
    outputVar.FactorA = A;
    outputVar.FactorB = B;
    outputVar.FactorC = C;
    outputVar.Spectra = SpectraAll;


end