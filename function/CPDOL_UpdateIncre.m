function [A,B,P,q,R,s,c] = CPDOL_UpdateIncre(A,B,P,q,R,s,Yt,lambda)

    mu = 1e-8;
        
    Frank = size(A,2); %CP-rank
    [I,J] = size(Yt);
    yvec = Yt(:);
    
    SampleIndex = find(yvec);

    Wmat = zeros(I,J);
    Wmat(SampleIndex) = 1;

    Wvec = Wmat(:);

    if nargin < 8 || isempty(lambda)
        lambda = 0.9;
    end

    
    DBA = Wvec.*kr(B,A);
    c = (DBA'*DBA + mu*eye(Frank))^(-1) * DBA'*yvec;
    
    BDC = B.* c';

    for ii = 1:I
        q{ii} = lambda*q{ii} + Yt(ii,:).*Wmat(ii,:)*BDC;
        P{ii} = lambda*P{ii} + BDC'.*Wmat(ii,:)*BDC;
        A(ii,:) = q{ii}*(P{ii} + mu*eye(Frank))^-1;
    end

    ADC = A.* c';

    for jj = 1:J
        s{jj} = lambda*s{jj} + ADC'.*Wmat(:,jj)'*Yt(:,jj);
        R{jj} = lambda*R{jj} + ADC'.*Wmat(:,jj)'*ADC;
        B(jj,:) = s{jj}'*(R{jj} + mu*eye(Frank))^-1;
    end


end