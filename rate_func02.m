%rate function for the second test IVP (eqn 7), undamped second order system
%X' = [0,-1;1,0]*X, X(0) = [1;0]
%INPUTS:
%t: the current time (not used, the system doesn't depend on t)
%X: the current state as a column vector [x1;x2]
%OUTPUTS:
%dXdt: the rate of change of X at (t,X)
function dXdt = rate_func02(t,X)
    dXdt = [0,-1;1,0]*X;
end
