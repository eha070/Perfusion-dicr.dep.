% ------------------------------------------------------------------------- 
%                               E98_TestPressure
% 
% Checks, if p fulfills the relationship -K/mu \Delta p = Q
% 
% 
% 
%                                          (c)Constantin Heck, 06-Apr-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clear
close all;
clc


%% prepare data
% load data
[prm,Fmat]   = settings;
basenameflow = perfusion1c.providenameflow(prm.phiopt,prm.Kopt,prm.dim);
pathload = ['results/synt-createflowTPFA-' basenameflow '.mat'];
msg = ['Loading ' pathload];
disp(msg);
D = load(pathload);

%get pressure
u     = D.pmat;

%get physiological dimensions
m     = size(u);
h     = prm.physdim(1:2)./m;

%setup permeabiltiy K and viscosity mu
K = D.Kmat(1,1,1,1);
mu = prm.mu;

%K  = 5e-6;
%mu = 5e-6; %i don't find mu in the code, the paper says its 5e-6

%% build laplacian

%function handle for identity
id = @(i) speye(m(i));

%get finite difference matrix for first derivative with neumann boundary
%conditions
D           = @(i) spdiags(ones(m(1)+1,1)*[-1,1],[-1,0],m(1)+1,m(1));


D1          = D(1)./h(1);
D1(1,1)     = 0;
D1(end,end) = 0;

D2          = D(2)./h(2);
D2(1,1)     = 0;
D2(end,end) = 0;

%get discrete gradient
Di   = kron(id(2),D1);
Dj   = kron(D2,id(1));
GRAD = [Di;Dj];
DIV  = GRAD';
LAP  = DIV*GRAD;


%% calculate -K/mu*\Delta u and visualize results
f = K/mu*(LAP*u(:));
f = reshape(f,m);
% Convert f (which we call Q in our paper) from m^3/s/m^3 to m^3/s which are the units of Fmat 
f = f*prod(prm.h);

figure(1); clf;
subplot(1,3,1);
imagesc(f);
axis image;
title('f')

subplot(1,3,2);
imagesc(Fmat);
axis image;
title('Fmat')

subplot(1,3,3);
imagesc(f - Fmat);
axis image;
title('f - Fmat')

