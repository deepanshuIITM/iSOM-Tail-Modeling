% TAIL MODELLING FOR DAMPER 
function [yd, Data]=mass_damp(n,xx)
R=0.01; J=0.01;y0=27;
if (nargin == 1)
mu = 1; sigma = 0.025; 
% pd = makedist('Normal','mu',mu,'sigma',sigma);
% tpd = truncate(pd,0.95,1.05);
% r1 = random(tpd,n,1);
% r2 = random(tpd,n,1);
r1 = normrnd(mu,sigma,n,1);
r2 = normrnd(mu,sigma,n,1);
elseif (nargin == 2)
r1 = xx(:,1);
r2 = xx(:,2);
end 

for i = 1:n
%     yd(i,1) = abs(1-(1/r2(i,1))^2)/sqrt((1-R*(1/r1(i,1))^2 -(1/r1(i,1))^2-(1/r2(i,1))^2 +...
%       (1/r1(i,1))^2*(1/r2(i,1))^2)^2 +4*(J^2)*((1/r1(i,1)) -(1/r1(i,1))*(1/r2(i,1)^2))^2) - y0;
  yd(i,1) = abs(1-(1/r2(i,1))^2)/(sqrt((1-R*(1/r1(i,1))^2 -(1/r1(i,1))^2-(1/r2(i,1))^2 +...
      (1/r1(i,1))^2*(1/r2(i,1))^2)^2 +4*(J^2)*((1/r1(i,1)) -(1/r1(i,1))*(1/r2(i,1)^2))^2)*y0);
end
Data = [r1 r2 yd];
end


