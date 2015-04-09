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
D = load(pathload);

%get pressure
p = reshape(D.pmat,[],1);
q = D.qmat;

%get physiological dimensions
m        = prm.dim(1:2);
h        = prm.h;
voxelvol = prod(prm.h);

K        = D.Kmat(1,1,1,1);
mu       = prm.mu;

%% build laplacian

%function handle for identity
id = @(i) speye(m(i));

%setup finite difference matrix for laplacian in 1D
Li = @(i) spdiags(ones(m(i)+1,1)*[1,-2,1],[-1,0,1],m(i),m(i)); 

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


%% get Gradient to setup the flow

%setup finite difference matrix for laplacian in 1D
Di = @(i) spdiags(ones(m(1)+1,1)*[-1,1],[-1,0],m(i)+1,m(i));


%get 1D-derivative-matrices with neumann boundary conditions
D1          = Di(1);
D1(1,1)     = 0;
D1(end,end) = 0;
D1          = D1./h(1);

D2          = Di(2);
D2(1,1)     = 0;
D2(end,end) = 0;
D2          = D2./h(2);

%kronecker to get the full Gradient
GRADy = kron(id(2),D1);
GRADx = kron(D2,id(1));



%% calculate -K/mu*\Delta u and visualize results
f = -K/mu*(LAP*p);
f = reshape(f,m);

% Convert from m^3/s/m^3 to m^3/s which are the units of Fmat
f = f*voxelvol;


%compare flow as well
qTest{1} = voxelvol*reshape(-GRADy*p,m(1)+1,m(2));
qTest{2} = voxelvol*reshape(-GRADx*p,m(1),m(2)+1);




return;

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