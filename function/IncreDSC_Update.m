function [A,B,c,s,R,q,P] =  IncreDSC_Update(A,B,c,s,R,q,P,Yt,lambda)
    
    mu = 1e-8;
    % svalue = 1e-16;
    [I,J] = size(Yt);
    Frank = size(A,2); 
    yvec = Yt(:);
    Wt = zeros(I,J);
    SampleIndex = find(yvec);
    yomega = yvec(SampleIndex);
    Wt(SampleIndex) = 1;


    %% update C
    
    krAB = kr(B,A);
    PartkrAB = krAB(SampleIndex,:);
    LC = 1/norm(PartkrAB'*PartkrAB);
    gC = (mu*c + (c*PartkrAB' - yomega') * PartkrAB);
    c = c - rand * LC * gC;
    % c(c<svalue) = svalue;
        
        %% update A
    BDC = B.* c;
    for ii = 1:I
        s{ii} = lambda * s{ii} + Yt(ii,:) .* Wt(ii,:) * BDC;
        R{ii} = lambda * R{ii} + BDC' .* Wt(ii,:) * BDC;
        
        aii = s{ii}*(R{ii} + mu*eye(Frank))^-1;


        A(ii,:) = aii;

    end
    % A(A<svalue) = svalue;


        %% update B
    ADC = A.* c;
    for jj = 1:J
        
        q{jj} = lambda * q{jj} + ADC' .* Wt(:,jj)' * Yt(:,jj);
        P{jj} = lambda * P{jj} + ADC' .* Wt(:,jj)' * ADC;
    
        bjj = q{jj}'*(P{jj} + mu*eye(Frank))^-1;


        B(jj,:) = bjj;

        
    end
    % B(B<svalue) = svalue;
    
end