function Y = AgleMat( X,location,ang1,ang2 )
%AGLEMAT  
%   X input location matrix
%   location: x + 1i*y
%  Y is 0-1 matrix
[m,n] = size(X);
Y = zeros(m,n);
for ii = 1:m
    for jj = 1:n
        ang = angle(X(ii,jj)-location);
        if ang>=ang1&ang<=ang2
            Y(ii,jj)=1;
        end
    end
end


end

