function [A,B,C] = BatchDSC_Update(Y,W,Rrank,iter,A0,B0,C0)
mu = 1e-8;
% svalue = 1e-16; % non-negative projection
warning off;

[I1,J1,K1] = size(Y);
[I2,J2,K2] = size(W);

if I1 ~= I2 || J1 ~= J2 || K1~= K2
    error('Sampling size does not match with data size');
end


if nargin < 5
    A0  = randn(I1,Rrank);
    B0 = randn(J1,Rrank);
    C0 = randn(K1,Rrank);
end




Y1 = tens2mat(Y,[],1);
Y2 = tens2mat(Y,[],2);
Y3 = tens2mat(Y,[],3);
W1 = tens2mat(W,[],1);
W2 = tens2mat(W,[],2);
W3 = tens2mat(W,[],3);


A = A0;
B = B0;
% % % A = randn(I,Rrank);
% % % B = randn(J,Rrank);
% % U = cpd(Y,Rrank);
% % A = U{1};
% % B = U{2};
% % C = U{3};
C = zeros(K1,Rrank);
C(1,:) = C0(2,:);
% C = C0;


for iitt = 1:iter
    %% updateC
    krAB = kr(B,A);
    

    for ki = 1:K1
        SampleIndexk = find(W3(:,ki));
        PartkrABk = krAB(SampleIndexk,:);
        yomega = Y3(SampleIndexk,ki);

        if iitt == 1 
            if ki >= 2
                ck = C(ki-1,:);
            else
                ck = C(ki,:);
            end
        else

            ck = C(ki,:);
        end

        
        LCk = 1/norm(PartkrABk);
       
        gCk = (mu * ck + (ck * PartkrABk' - yomega') * PartkrABk);
        ck = ck - rand*(LCk^2) * gCk;
%         end
        C(ki,:) = ck;
    end
   % C(C<svalue) = svalue;
    %% update A;
    aiCB = kr(C,B)';
    for ii = 1:I1
        
        logg = W1(:,ii)==1;
        aisparse = zeros(size(aiCB));
        aisparse(:,logg) = aiCB(:,logg);
        ai = (aisparse*aisparse' + mu*eye(Rrank) )\( aisparse*Y1(:,ii) );
        A(ii,:) = ai';

    end
   % A(A<svalue) = svalue;
    %% update B
    bjCA = kr(C,A)';

    for jj = 1:J1

        logg = W2(:,jj)==1;
        bjsparse = zeros( size(bjCA) );
        bjsparse(:,logg) = bjCA(:,logg);
        bj = ( bjsparse*bjsparse' + mu*eye(Rrank) )\( bjsparse*Y2(:,jj) );
        B(jj,:) = bj';
        
    end
   % B(B<svalue) = svalue;
end


end

