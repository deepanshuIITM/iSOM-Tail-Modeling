function [Cx, Cy, R] = curvature(ln_TPNT,yd)
for i = 1:4
    g = [1 100 175 200 300];
y = ln_TPNT(g(i):min(g(i+1),length(yd))); x = yd(g(i):min(g(i+1),length(yd)));
fcn = @(b) norm(((x(:)-b(1)).^2 + (y(:)-b(2)).^2) - b(3).^2);   % Objective Function
B = fminsearch(fcn, rand(3,1));                                 % Estimate Parameters
Cx = B(1);                                                      % Center X-Coordinate
Cy = B(2);                                                     % Center Y-Coordinate
R  = B(3);                                                    % Radius
% B = fsolve(fcn, rand(3,1))
plot(x, y, 'gs','markerSize',2.5)
hold on
plot(Cx, Cy, 'bs','markerSize',2.5)
plot([Cx x(1)], [Cy y(1)], 'b-.')
plot([Cx x(end)], [Cy y(end)], 'b-.')
end
end