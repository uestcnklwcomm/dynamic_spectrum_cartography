function [ x,y ] = RectCircle(l,a,b,p )
%RECTCIRCLE total length of path
%   rectangle length and width
%  p = 1: anti-clockwise
%  p = 0: clockwise
if nargin<4
    p = 1;
end
m = floor(l/(a+b));
if p==1
    if mod(m,2)==0
        x = min(l-m*(a+b),a);
        y = max(l-m*(a+b)-a,0);
    else
        x = a - min(l-m*(a+b),a);
        y = b - max(l-m*(a+b)-a,0);
    end
    
    
    
else
    if mod(m,2)==0
        x = max(l-m*(a+b)-b,0);
        y = min(l-m*(a+b),b);
    else
        x = a - max(l-m*(a+b)-b,0);
        y = b - min(l-m*(a+b),b);
    end
    
end

end

