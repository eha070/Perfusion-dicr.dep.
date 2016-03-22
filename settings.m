function [prm,Fmat] = settings

% Make experiment folder
prm.experiment = '3D-flow';
prm.resultfolder = ['results-' prm.experiment];
[a,b] = mkdir(prm.resultfolder);
[a,b] = mkdir(['figs-' prm.resultfolder]);

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

% physical dimension in meter
if isequal(prm.experiment,'3D-flow')
    prm.physdim = [3,3,0.5]*milli*meter;
else
    prm.physdim = [10,10,1]*milli*meter;
end;
prm.physvol = prod(prm.physdim);


% dimension of problem
if isequal(prm.experiment,'3D-flow')
    prm.dim = [64,64,3];
    prm.capslice = 2;
else
    prm.dim = [64,64,1];
end;

% prm.dim = [4,4,1];
dim = prm.dim;

% voxel size
prm.h = prm.physdim./prm.dim;

% voxel volume
prm.voxelvol = prod(prm.h);

% viscosity blood in Pa*s
% Source: https://en.wikipedia.org/wiki/Viscosity
prm.mu = 3.5*1e-3;

% timeline and dt in seconds
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
if isequal(prm.experiment,'3D-flow')
    a(1:10,1:10,1) = 1;
else
    a(1,1) = 1;
end;
ind = find(a);
[cso(:,1),cso(:,2),cso(:,3)] = ind2sub(dim,ind);
prm.cso = cso;
prm.nso = size(prm.cso,1);


a = zeros(dim);
if isequal(prm.experiment,'3D-flow')
    a(end-9:end,end-9:end,3) = 1;
else
    a(end,end) = 1;
end;
ind = find(a);
[csi(:,1),csi(:,2),csi(:,3)] = ind2sub(dim,ind);
prm.csi = csi;
prm.nsi = size(prm.csi,1);

% Perfusion of 5-70ml/min/100ml (From constantin)
% Go for a perfusion of 50ml/min/100ml 
qin = 50;
prm.perfusionin = qin;

% This corresponds to 0.5ml/min/ml = 0.5m^3/min/m^3
qin = qin/100;

% Multiply by volume of our slab to get the total flow into the region to
% get m^3/min
qin = qin*prm.physvol;

% Divide by 60 to get per second: m^3/sec
qin = qin/60;

% Fmat = int_Omega_i Q(x) dx where Q = [m^3/s/m^3], so Fmat = [m^3/s],
% the absolute inflow
% Assume equal inflow in all sources and sinks
Fmat = zeros(dim);
qout = -qin;
for i = 1 : prm.nso
    Fmat(prm.cso(i,1),prm.cso(i,2),prm.cso(i,3)) = qin/prm.nso;
    Fmat(prm.csi(i,1),prm.csi(i,2),prm.csi(i,3)) = qout/prm.nsi;
end;

% % Fluid density in units kg/m^3
% % Value for blood according to http://hypertextbook.com/facts/2004/MichaelShmukler.shtml
% NB fluid density is already scaled away. The original source term is in
% mg/s/m^3. By dividing with rho we get m^3/s/m^3
% prm.rho = 1025;

% % scale away the density
% Fmat = Fmat/prm.rho;

% for plotting Cmat
prm.scaling = 1;



