function X4DT = RadioMapGenerator(Config)

    if ~isfield(Config,'EmitterNumber')
        error('Number of emitters undefined!')
    end

    R = Config.EmitterNumber; 
    
    if ~isfield(Config,'width') || ~isfield(Config,'length') || ~isfield(Config,'freq') || ~isfield(Config,'time')
        warning('Incomplete scenario configurations; set by default as 51*51*64*600')
        
        gridLen = 50;
        gridWid = 50;
        K = 2^6;
        T = 600;
    
    else 
        gridLen = Config.length - 1;
        gridWid = Config.width - 1;
        K = Config.freq;
        T = Config.time;
    end

    if ~isfield(Config, 'Directional')
        directional = 0;
    else
        directional = Config.directional;  %0 for omin-directional; 1 for directional;
    end
    
    if ~isfield(Config, 'speed')
        v = 0.01;
    else
        v = Config.speed;
    end
    
    if ~isfield(Config, 'sigma')
        sigma_s = 8;
    else
        sigma_s = Config.sigma;
    end

    if ~isfield(Config, 'decorr')
        dc = 50;
    else
        dc = Config.decorr;
    end

    
    loss_f = @(x,d,alpha) min(1,(x/d).^(-alpha)); %synthetic pathloss component
    alpha = 2;%path loss coefficient
    d0 = 2;
    %% Generating spatial meshgrids
    gridResolution = 1;%
    x_grid = 0:gridResolution:gridLen;
    y_grid = 0:gridResolution:gridWid;
    [Xmesh_grid, Ymesh_grid] = meshgrid(x_grid, y_grid);
    Xgrid = Xmesh_grid + 1i*Ymesh_grid;
    [I,J] = size(Xgrid); 
    
    
    %% Generating PSDs
    indK = (1:K)';
    Sx =@(f0,a) sinc((indK-f0)/a).^2.*( abs((indK-f0)/a)<=1);  % basis sinc wave function
    Ctrue = zeros(K,R);
    ind_psd = 2:2:K-2;
    
    for rr = 1:R
        am = 1.5 + 0.5*rand(3,1); % random amplitude
        centerf0 = ind_psd(randperm(length(ind_psd),3));
        cr = am(1)*Sx(centerf0(1),3*(1+rand)) + (rand>0.5)*am(2)*Sx(centerf0(2),3*(1+rand)) + (rand>0.5)*am(3)*Sx(centerf0(3),3*(1+rand));
        Ctrue(:,rr) = cr;
    end
    
    
    
    
        %% Generating emitters' traces
    for rr = 1:R
        start_location_rr = rand*50 + 1j*rand*50;
        clock_wise_rr = (rand < .5);
        rect_length_rr = rand*(50 - real(start_location_rr));
        rect_width_rr = rand*(50 - imag(start_location_rr));
        for tt = 1:T
            [xt,yt] = RectCircle(v*tt,rect_length_rr,rect_width_rr,clock_wise_rr);
            xi{rr}(tt) = real(start_location_rr) + xt;
            yi{rr}(tt) = imag(start_location_rr) + yt;
        end
        if directional %generating radius angle
            Ang1{rr} = 0.5*rand*pi;
            Ang2{rr} = Ang1{rr} + 0.5*pi + 0.5*rand*pi;
            Ang2{rr} = max(Ang2{rr},pi);
        end       
        location_set{rr} = xi{rr} + 1i*yi{rr};
    end
    
    %% Generating dynamic SLFs
    
    X4DT = zeros(I,J,K,T);

    for rr = 1:R
        shadowdB = Shadowing(Xgrid,sigma_s,exp(-1/dc)); %larger sigma_s for heavier shadowing
        shadow_linear{rr} = 10.^(shadowdB/10);
    end
    
    for tt = 1:T
        Xt = zeros(I,J,K);
        for rr = 1:R
            location = location_set{rr}(tt);
            loss_mat = abs(Xgrid - location);
            Sr = loss_f(loss_mat,d0,alpha).*shadow_linear{rr};
            Sr = Sr/max(max(Sr));
            Sr = Sr/norm(Sr,'fro');
            if directional
                Angrr = Algemat(Xgrid,location,Ang1{rr},Ang2{rr});
                Sr = Sr.*Angrr;
            end
            Xt = Xt + outprod(Sr, Ctrue(:,rr));
        end
        X4DT(:,:,:,tt) = Xt;
    end

end