function [ r ] = nsolve( f, fp, e)
%NSOLVE Newton's method to find roots of the function f
%   Detailed explanation goes here
    
    xroot = 20; %TODO: estimate it better
    while abs(f(xroot)) > e
        xroot = -f(xroot)./fp(xroot);
    end;

    r = xroot;

end

