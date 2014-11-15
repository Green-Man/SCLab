function [ x r ] = nsolve( f, fp, x0, e )
%NSOLVE Newton's method to find roots of the function f
%   Detailed explanation goes here
    r = true;
    x = x0;
    c = 0;
    while abs(f(x)) > e
        x = x - (f(x))./fp(x);
        c = c + 1;
        if c>1000
            r = false;
            break;
        end;
    end;

end


