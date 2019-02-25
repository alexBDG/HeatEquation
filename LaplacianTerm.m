function [fuxU] = LaplacianTerm(a,u,dx,cas)
% step forward the equation u_t + fu_x using finite differences on
% the grid given by dx, dt with speed f
% using Neumann bc
% using different choices of finite differences

% shift the vectors using the appropriate BC
[uB, uF, aB] = shiftNeumann(u,a);
% Now define the forward differences
fuxF = (uF-u)/dx;
% Backward differences
fuxB = (u-uB)/dx;
% centered differences
if cas==2
    fuxU = (a.*fuxF-aB.*fuxB)/(dx);
else
    fuxU = a.*(fuxF-fuxB)/(dx);
end
    


function [uB, uF, aB] = shiftNeumann(u,a)
u = u(:); n = length(u);
a = a(:);
% duplicating the first term leads to u_x = 0 at left boundary
uB = [u(1); u(1:n-1)];
uF = [u(2:n); u(n)];
aB = [a(1); a(1:n-1)];

function [uB, uF, aB] = shiftDirichlet(u,a)
u = u(:); n = length(u);
a = a(:);
% zero Dirichlet BC at endpoints
uB = [0; u(1:n-1)];
uF = [u(2:n);0];
aB = [a(1); a(1:n-1)];