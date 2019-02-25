% script to solve the variational problem
% min 1/n sum |y_i - u_i|^2 + \lambda \max_i |D_x u_i|
% with unstructured data points
clear; clc
nx = 200; 
x = linspace(-5,5,nx); dx = x(2)-x(1);
x = x(:);
%% Initial Function
%u0 = max(0*x, 1-x.^2);
u0 = max(0*x, 1-abs(x));



%% The Coeff Heat equation
a = ones(nx,1)*1.5;
abar = max(abs(a));
%% time step
Tf = .5;
dtmax = .1*(dx^2)/(abar);
Nt = Tf/dtmax;  Nt = ceil(Nt);
dt = .1*dx^2/abar;
%% Iterate
% Trick to keep number of plots small:
nplots = 6;
aa = max(1,floor(Nt/nplots));

figure(1), plot(x,u0,'k'); pause(0.5)
hold on
u1 = u0;
for jj = 1:Nt
    [fuxU] = LaplacianTerm(a,u1,dx,1);
    u1 = u1 + dt*fuxU;
    if mod(jj,aa) == 0 % whether to plot, only want to plot a few times
        figure(1), plot(x,u1,'b-'); pause(0.2)
    end
    % better to plot a few times
end
hold off


%% The Coeff Heat equation
a = 1 + .9*sin(pi*x);
abar = max(abs(a));
%% time step
Tf = .5;
dtmax = .1*(dx^2)/(abar);
Nt = Tf/dtmax;  Nt = ceil(Nt);
dt = .1*dx^2/abar;
%% Iterate
% Trick to keep number of plots small:
nplots = 6;
aa = max(1,floor(Nt/nplots));

figure(2), plot(x,u0,'k'); pause(0.5)
hold on
u2 = u0;
for jj = 1:Nt
    [fuxU] = LaplacianTerm(a,u2,dx,1);
    u2 = u2 + dt*fuxU;
    if mod(jj,aa) == 0 % whether to plot, only want to plot a few times
        figure(2), plot(x,u2,'b-'); pause(0.2)
    end
    % better to plot a few times
end
hold off


%% The Coeff Heat equation
a = 1 + .9*sin(pi*x);
abar = max(abs(a));
%% time step
Tf = .5;
dtmax = .1*(dx^2)/(abar);
Nt = Tf/dtmax;  Nt = ceil(Nt);
dt = .1*dx^2/abar;
%% Iterate
% Trick to keep number of plots small:
nplots = 6;
aa = max(1,floor(Nt/nplots));

figure(3), plot(x,u0,'k'); pause(0.5)
hold on
u3 = u0;
for jj = 1:Nt
    [fuxU] = LaplacianTerm(a,u3,dx,2);
    u3 = u3 + dt*fuxU;
    if mod(jj,aa) == 0 % whether to plot, only want to plot a few times
        figure(3), plot(x,u3,'b-'); pause(0.2)
    end
    % better to plot a few times
end
hold off


figure(4), plot(x,u1,'b',x,u2,'b',x,u3,'b');
