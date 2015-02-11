function [prm,Qmat] = settings

% reduce time steps
prm.stepred = 20;

% type of porosity
prm.phiopt = 'flat';

% type of permeability
prm.Kopt = 'flat';

prm.aiftype = 'gamma';
% prm.aiftype = 'parker';
% prm.aiftype = 'delta';

prm.physdim = [10,10,0.1];

% dimension of problem
%prm.dim = [32,32,1];
% prm.dim = [64,64,1];
prm.dim = [128,128,1];
% prm.dim = [256,256,1];
% prm.dim = [512,512,1];

% prm.dim = [4,4,1];
dim = prm.dim;

% voxel size in mm
prm.h = prm.physdim./prm.dim;
%prm.h = [1,1,1]*0.1;

% voxel volume in ml
prm.voxelvol = prod(prm.h);

% viscosity blood in kPa*s
prm.mu = 5*1e-6;

% timeline and dt
% prm.T = 240;
prm.T = 480;
% prm.T = 40;
prm.dt = 0.01;
% prm.aifpeak = 7;
% prm.aifwidth = 0.3;
prm.timeline = (0:prm.dt:prm.T);


%
% Source and sink
%

prm.cso = [1,1,1];
% prm.cso = [1,1,1;...
%     2,1,1;...
%     1,2,1;...
%     2,2,1];
% prm.cso(:,1) = ones(100,1);
% prm.cso(:,2) = (1:100)';
% prm.cso(:,3) = 1;
prm.nso = size(prm.cso,1);

%prm.csi = [dim(1),1,1];
prm.csi = [dim(1),dim(2),1];
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
% Comment: 50/60/100 is 0.0083ml/s/ml = 0.0083mm^3/s/mm^3 as average flow in the
% brain
% qin = 50;
% 
% % scale to ml/sec/100ml = mm^3/s/100mm^3
% qin = qin/60;
% 
% % scale to mm^3/s/mm^3
% qin = qin/100;

% this value gives a perfusion of around 50ml/min/100ml!!! Thats why we use
% it.
qin = 0.0007;

% Actually, must multiply Q with a unit converter a = 1(ml/mm^3), becoming
% Q*a/rho
Qmat = zeros(dim);

% % volume of THIS brain in mm^3: 
% brainvol = prod(dim)*prm.voxelvol;
% 
% % scale the flow to feed the brain in total (mm^3/s)
% qin = qin*brainvol;

% scaling it to the voxel volume to have a flow mm^3/s/mm^3 in the source
% shouldn't this be here????
% qin = qin/prm.voxelvol; 


qout = -qin;
for i = 1 : prm.nso
    Qmat(prm.cso(i,1),prm.cso(i,2),prm.cso(i,3)) = qin;
    Qmat(prm.csi(i,1),prm.csi(i,2),prm.csi(i,3)) = qout;
end;

% fluid density map in units mg/mm^3
prm.rho = 1;

% scale away the density
Qmat = Qmat/prm.rho;


