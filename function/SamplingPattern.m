function outputConfig = SamplingPattern(Config)

    I = Config.length;
    J = Config.width;
    T = Config.time;
    F = Config.CPrank;
    % R = Config.EmitterNumber;
    rho = Config.SamplingRatio;
    flag = Config.pattern; % 1: each column uniformly at random; 2. coupled sampling pattern

    SampleIndexAll = cell(1,T);
    Wtens = zeros(I,J,T);

    switch flag
        case 1
            for tt = 1:T
                Wmatt = zeros(I,J);

                for jj = 1:J
                    idx_j = randperm(I,round(I*rho));
                    Wmatt(idx_j,jj) = 1;
                end
                SampleIndexAll{tt} = find(Wmatt(:));
                Wtens(:,:,tt) = Wmatt;
            end

            outputConfig.idx = SampleIndexAll;
            outputConfig.tensor = Wtens;
            outputConfig.Batchsize = round(1/rho);

        case 2
            %% generating patterns for CPD algebraic intialization; number of patterns D = 2
            
            
            alg_mask = zeros(I,J);
            
            

            row_comb = nchoosek(1:I,2);
            col_comb = nchoosek(1:J,2);
            num_subtensor = nchoosek(F,2);

            row_selection = row_comb([randperm(size(row_comb,1),num_subtensor)],:);
            col_selection = col_comb(randperm(size(col_comb,1),num_subtensor)',:);

            for kk = 1:num_subtensor
                alg_mask([row_selection(kk,:)],[col_selection(kk,:)]) = 1;
            end

            
            
            %% sampling pattern for algebraically initializing factor A and B

            if isfield(Config,'AlgInit') && Config.AlgInit

                coupled_mask = zeros(I,J);
                Avail_col_idx = 1:J;
                sel_col = randperm(numel(Avail_col_idx),2);
                coupled_col = Avail_col_idx(sel_col);
                coupled_col_idx = cell(1,I);
                coupled_col_idx{1} = coupled_col;
                coupled_mask(1,sel_col) = 1;
    
                for ii = 2:I
                    Avail_col_idx = 1:J;
                    sel_coupled_col = randperm(numel(coupled_col),1);
                    inherit_coupled_col = coupled_col(sel_coupled_col);
                    Avail_col_idx(coupled_col) = [];
    
                    if numel(Avail_col_idx) < 1
                        new_col = randperm(J,1);
                    else
                        sel_new_col = randperm(numel(Avail_col_idx),1);
                        new_col = Avail_col_idx(sel_new_col);
                    end
    
                    coupled_col = [coupled_col,new_col];
                    coupled_col = unique(coupled_col);
                    coupled_col_idx{ii} = [inherit_coupled_col,new_col];
    
                    coupled_mask(ii,coupled_col_idx{ii}) = 1;
           
                end

                coupled_idx = find(coupled_mask(:));
                alg_mask(coupled_idx) = 1;

                outputConfig.alg_row = row_selection;
                outputConfig.alg_col = col_selection;

                outputConfig.AlgebraicMask = alg_mask;

            end

            

            alg_idx = find(alg_mask(:));            

            AvaiIndex = 1:I*J;

            AvaiIndex(alg_idx) = [];

            SampPerIt = round(I*J*rho) - numel(alg_idx);     
            SampleIndexAcc = [];
            
            bb = 1;
            
            while numel(AvaiIndex) > 0

                if SampPerIt > numel(AvaiIndex)
                    IndexNew = 1:numel(AvaiIndex);
                    idx_Retain = randperm(numel(SampleIndexAcc),SampPerIt - numel(AvaiIndex));
                    SampleNew = [SampleIndexAcc(idx_Retain),AvaiIndex];
                                
                else

                    IndexNew = randperm(numel(AvaiIndex),SampPerIt);
                    SampleNew = AvaiIndex(IndexNew);

                end
                
                AvaiIndex(IndexNew) = [];
            
                SampleIndex{bb} = SampleNew;
                SampleIndexAcc = [SampleIndexAcc,SampleNew];

                bb = bb + 1;
            end

            cyc_len = length(SampleIndex);

            Batchsize = max(cyc_len,F);
            outputConfig.Batchsize = Batchsize;

      


            for tt = 1:T
                cyc_idx = mod(tt,cyc_len) + cyc_len*(mod(tt,cyc_len) == 0);
                Wmatt = zeros(I,J);
                Wmatt(alg_idx) = 1;
                idx_t = SampleIndex{cyc_idx};
                Wmatt(idx_t) = 1;

                SampleIndexAll{tt} = find(Wmatt(:));
                SampleIndexAll{tt} = unique(SampleIndexAll{tt});
                Wtens(:,:,tt) = Wmatt; 
                
            end

            outputConfig.tensor = Wtens;


           
            
            
    end

end
