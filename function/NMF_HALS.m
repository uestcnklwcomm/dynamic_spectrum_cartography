function [A,B ] = NMF_HALS( X,R,iter,Ai,Bi )
%NMF_HALS 此处显示有关此函数的摘要
%   此处显示详细说明
[m,n] = size(X);
if nargin<4
    A = randn(m,R).^2;
    B = randn(n,R).^2;
else
    A = Ai;
    B = Bi;
end
lambda = 1e-10;
for ii = 1:iter
    for rr=1:R
        ar = A(:,rr);
        br = B(:,rr);
        Yr = X-A*B'+ar*br';
        ar =Yr*br/( norm(br)^2+ lambda) ;
        ar(ar<0)=0;
        br =Yr'*ar/( norm(ar)^2 + lambda);
        br(br<0)=0;
        
        A(:,rr) = ar;
        B(:,rr) = br;
        
    end
end

end

