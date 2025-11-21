function [sys_trim] = trim_tf(sys, tol)

if nargin<2
    tol = 1e-5;
end
[n,d] = tfdata(sys);
n = cellfun(@(x) {x.*(abs(x)>tol)}, n);
d = cellfun(@(x) {x.*(abs(x)>tol)}, d);
sys_trim = tf(n, d, sys.Ts);
end

