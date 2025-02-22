function U = CPDOL_UpdateBatch(Y,F,iter,Batchsize,U)
    
    [I,J,T] = size(Y);
    mu = 1e-8;

    if T < Batchsize
        Batchsize = T;
    end
    
    Y = Y(:,:,end-Batchsize+1:end);

    Wtens = zeros(I,J,Batchsize);

    for bb = 1:Batchsize
        Wmatb = zeros(I,J);
        Yb = squeeze(Y(:,:,bb));
        SampleIndexb = find(Yb(:));
        Wmatb(SampleIndexb) = 1;
        Wtens(:,:,bb) = Wmatb;
    end
     
    W1 = tens2mat(Wtens,[],1);
    W2 = tens2mat(Wtens,[],2);
    W3 = tens2mat(Wtens,[],3);

    Y1 = tens2mat(Y,[],1);
    Y2 = tens2mat(Y,[],2);
    Y3 = tens2mat(Y,[],3);
    

    if nargin < 5 || isempty(U)
        A = randn(I,F);
        B = randn(J,F);
    else
        A = U{1};
        B = U{2};
    end
    
    C = zeros(Batchsize,F);

    for it = 1:iter

        cbBA = kr(B,A)';
        for bb = 1:Batchsize

            logg = W3(:,bb)==1;
            cbsparse = zeros(size(cbBA));
            cbsparse(:,logg) = cbBA(:,logg);
            cb = (cbsparse*cbsparse' + mu*eye(F) )\( cbsparse*Y3(:,bb) );
            C(bb,:) = cb';
        end

        aiCB = kr(C,B)';

        for ii = 1:I

            logg = W1(:,ii)==1;
            aisparse = zeros(size(aiCB));
            aisparse(:,logg) = aiCB(:,logg);
            ai = ( aisparse*aisparse' + mu*eye(F) )\( aisparse*Y1(:,ii) );
            A(ii,:) = ai';

        end

        bjCA = kr(C,A)';

        for jj = 1:J

            logg =W2(:,jj)==1;
            bjsparse =zeros( size(bjCA) );
            bjsparse(:,logg) = bjCA(:,logg);
            bj = ( bjsparse*bjsparse' + mu*eye(F))\( bjsparse*Y2(:,jj) );
            B(jj,:) = bj';


        end

        
    end

    U{1} = A;
    U{2} = B;
    U{3} = C;
    
end