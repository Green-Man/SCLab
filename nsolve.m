function [ x ] = nsolve( f, fp, e)
%NSOLVE Newton's method to find roots of the function f
%   Detailed explanation goes here
    
    x = 20; %TODO: estimate it better
    while abs(f(x)) > e
        x = (fp(x)*x-f(x))./fp(x);
    end;

end

