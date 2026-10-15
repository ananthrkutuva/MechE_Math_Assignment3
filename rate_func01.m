%rate function for the first test IVP (eqn 5)
%x' = -5x + 5cos(t) - sin(t), x(0) = 1
%INPUTS:
%t: the current time
%X: the current value of x
%OUTPUTS:
%dXdt: the rate of change of x at (t,X)
function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end
