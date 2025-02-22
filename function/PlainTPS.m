function X4DHat = PlainTPS(SampleData,SamplingTensor)
    if isempty(SampleData)
        error('Input shape/data not detected.');
    end

    Yall = SampleData;
    
    K = size(Yall{1},2);
    T = length(Yall);

    if isempty(SamplingTensor)
        error('Sampling mask not detected');
    end

    Wtens = SamplingTensor;

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
    
    Xhat_tps_tt = zeros(I,J,K);
    X4DHat = zeros(I,J,K,T);
    

    for tt = 1:T
        Yt = Yall{tt};
        Wmatt = squeeze(Wtens(:,:,tt));
        SampleIndext = find(Wmatt(:));
        x_IND = real(Xgrid(SampleIndext));
        y_IND = imag(Xgrid(SampleIndext));
        
        for kk = 1:K
            ykk = Yt(:,kk);
            stps = TPS(x_IND,y_IND,ykk,xgrid,ygrid,lambda_tps,dB);
            Xhat_tps_tt(:,:,kk) = reshape(stps,[I,J]);
        end

        X4DHat(:,:,:,tt) = Xhat_tps_tt;
    end

end