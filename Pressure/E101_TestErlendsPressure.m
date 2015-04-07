% ------------------------------------------------------------------------- 
%                               E98_TestErlendsPressure
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
D = load(pathload);

%get pressure
u = D.pmat;

%get physiological dimensions
m = size(u);
h = prm.h;

%setup permeabiltiy K and viscosity mu
K  = 5e-6;
mu = prm.mu;

%% build laplacian

%function handle for identity
id = @(i) speye(m(i));

%setup finite difference matrix for laplacian in 1D
Li = @(i) spdiags(ones(m(1)+1,1)*[1,-2,1],[-1,0,1],m(1),m(1)); 


%get 1D-derivative-matrices with neumann boundary conditions
L1          = Li(1);
L1(1,1)     = -1;
L1(end,end) = -1;
L1          = L1./h(1)^2;

L2          = Li(2);
L2(1,1)     = -1;
L2(end,end) = -1;
L2          = L2./h(2)^2;

%kronecker to get the full Laplacian
LAP = kron(id(2),L1) + kron(L2,id(1));



%% calculate -K/mu*\Delta u and visualize results
f = -K/mu*(LAP*u(:));
f = reshape(f,m);


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