% ----- Bisection -----

function r = bisection(f, a, b, tol)
%   Returns NaN if f(a) and f(b) have the same sign (no root bracketed).
if nargin < 4
    tol = 1e-6;
end
if sign(f(a)) == sign(f(b))
    r = NaN;
    return
end
while (b - a)/2 > tol
    c = (a + b)/2;
    if sign(f(c)) == sign(f(a))
        a = c;       % root is in the right half
    else
        b = c;       % root is in the left half
    end
end
r = (a + b)/2;
end
