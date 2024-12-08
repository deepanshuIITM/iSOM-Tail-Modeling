function p3 = weig_lsqfit_constr(x,y,w)
%--------------------------------------------------------------------------
input = [ones(length(x),1) x x.^2 x.^3];
alpha0 = pinv(input)*y;
w = (w-min(w))./range(w);
obj= @(alpha) sum((w.^1)*(y- input*alpha).^2);

function [c,ceq] = constr(alpha)
c = [-alpha(4)-10^-5;((alpha(3))^2 - 3*alpha(2)*alpha(4)-10^-5)];
ceq = [];
end
p3 = fmincon(obj,alpha0,[],[],[],[],[],[],@constr);
%--------------------------------------------------------------------------

end