function [prm,Fmat] = settings

prm.reconflowIO = 1e-1;
prm.reconflowalpha = 0;
prm.reconflowregularization = '';
prm.reconflowroitype = 'full';
prm.reconflowcode = 'classic';
prm.reconflowdata = 'synt';
prm.reconflowdeconvmethod = 'circSVD';

% reduce time steps
prm.stepred = 100;

% type of porosity
prm.phiopt = 'flat';

% type of permeability
prm.Kopt = 'flat';

prm.aiftype = 'gamma';
% prm.aiftype = 'parker';
% prm.aiftype = 'delta';

prm.physdim = [10,10,1];
prm.physvol = prod(prm.physdim);

% dimension of problem
%prm.dim = [32,32,1];
prm.dim = [64,64,1];
% prm.dim = [20,20,1];
% prm.dim = [128,128,1];
% prm.dim = [256,256,1];
% prm.dim = [512,512,1];

% prm.dim = [4,4,1];
dim = prm.dim;

% voxel size in mm
prm.h = prm.physdim./prm.dim;

% voxel volume in ml
prm.voxelvol = prod(prm.h);

% viscosity blood in kPa*s
prm.mu = 5*1e-6;

% timeline and dt
prm.T = 90;
% prm.T = 480;
% prm.T = 40;
prm.dt = 0.002;
% prm.aifpeak = 7;
% prm.aifwidth = 0.3;
prm.timeline = (0:prm.dt:prm.T);


%
% Source and sink
%
a = zeros(dim);
a(1,1) = 1;
ind = find(a);
[cso(:,1),cso(:,2),cso(:,3)] = ind2sub(dim,ind);
prm.cso = cso;

% prm.cso = [1,1,1];
% prm.cso = [1,1,1;...
%     2,1,1;...
%     1,2,1;...
%     2,2,1];
% prm.cso(:,1) = ones(100,1);
% prm.cso(:,2) = (1:100)';
% prm.cso(:,3) = 1;
prm.nso = size(prm.cso,1);

a = zeros(dim);
a(end,end) = 1;
ind = find(a);
[csi(:,1),csi(:,2),csi(:,3)] = ind2sub(dim,ind);
prm.csi = csi;

%prm.csi = [dim(1),1,1];
% prm.csi = [dim(1),dim(2),1];
% prm.csi = [dim(1),1,1;...
%     dim(1)-1,1,1;...
%     dim(1),2,1;...
%     dim(1)-1,2,1];
% prm.csi(:,1) = 100*ones(100,1);
% prm.csi(:,2) = (1:100)';
% prm.csi(:,3) = 1;
prm.nsi = size(prm.csi,1);

% Perfusion of 5-70ml/min/100ml (From constantin)
% Go for a perfusion of 50ml/min/100ml 
qin = 50;

% This corresponds to 0.5ml/min/ml = 0.5mm^3/min/mm^3
qin = qin/100;

% Multiply by volume of our slab to get the total flow into the region to
% get ml/min
qin = qin*prm.physvol;

% Divide by 60 to get per second: mm^3/sec
qin = qin/60;

% this value gives a perfusion of around 50ml/min/100ml!!! Thats why we use
% it.
% qin = 0.0007;

% Actually, must multiply Q with a unit converter a = 1(ml/mm^3), becoming
% Q*a/rho
Fmat = zeros(dim);

qout = -qin;
for i = 1 : prm.nso
    Fmat(prm.cso(i,1),prm.cso(i,2),prm.cso(i,3)) = qin/prm.nso;
    Fmat(prm.csi(i,1),prm.csi(i,2),prm.csi(i,3)) = qout/prm.nsi;
end;

% fluid density map in units mg/mm^3
prm.rho = 1;

% scale away the density
Fmat = Fmat/prm.rho;

% for plotting Cmat
prm.scaling = 1;



