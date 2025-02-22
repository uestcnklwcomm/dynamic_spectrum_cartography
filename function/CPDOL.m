function [X4DT,outputVar] = CPDOL(Input,AlgConfig)

%% loading scenario configurations
    if ~isfield(Input, 'length') || ~isfield(Input, 'width') || ~isfield(Input, 'data')
        error('Input size/data not dectected');
    end

    I = Input.length;
    J = Input.width;

    Yall = Input.data; % with size |Omega_t| * |K| * |T|;
    K = size(Yall{1},2);
    T = length(Yall);

    if ~isfield(Input,'Batchsize')
        error('Batchsize not detected.');
    end

    Batchsize = Input.Batchsize;
    
    if ~isfield(AlgConfig,'CPrank')
        error('Hyperparameter CP-rank undefined.');
    end

    
    F = AlgConfig.CPrank;

    if ~isfield(AlgConfig,'Iteration')
        ALS_iter = 10;
    else
        ALS_iter = AlgConfig.Iteration;
    end

    
    for kk = 1:K
        A{kk} = randn(I,F);
        B{kk} = randn(J,F);
        C{kk} = randn(Batchsize,F);
    end

    Wtens = Input.SamplingTensor;

    for tt = 1:T
        Wmatt = squeeze(Wtens(:,:,tt));
        for kk = 1:K
            
        end
    end
    % SampleIndexAll = Input.SampleIndex;

    
end