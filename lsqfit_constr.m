function p3 = lsqfit_constr(x,y)
%--------------------------------------------------------------------------
input = [ones(length(x),1) x x.^2 x.^3];
alpha0 = pinv(input)*y;

obj= @(alpha) sum((y- input*alpha).^2);

A = [0 0 0 -1];
b = +1e-4;
function [c,ceq] = constr(alpha)
c = ((alpha(3)^2) -3*alpha(2)*alpha(4) -1e-5);
ceq = [];
end
optOptions = struct('MaxFunctionEvaluations', 10000000);
p3 = fmincon(obj,alpha0,A,b,[],[],[],[],@constr,optOptions);

end
%--------------------------------------------------------------------------

