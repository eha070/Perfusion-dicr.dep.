function [pmat,flowmat] = syntTPFA(Kmat,Fmat,prm)
%function [pmat,flowmat] = syntTPFA(Kmat,Fmat,prm)
%
%INPUT:
% Kmat - cell-centered diagonal tensor [mm^2]
% Fmat - cell-centered absolute in-/outflow in mm^3/s
% prm  - parameters with fields
%        dim - dimension
%        h   - voxel-sizes [mm]
%        mu  - viscosity  [Pa*s]
%
%OUTPUT:
% pmat    - pressure map [Pa]. Note that this is only defined up to a 
%           constant.
% flowmat - map with absolute flow [mm^3/s].
%
%
%MOTIVATION FOR TPFA:
%We want to avoid differntiation across voxels and write the flux between 
%voxel 1 and 2 as (see figure below for layout):
%
%   F = -|gamma_12|*lambda_1*(p12-p1)/(dx/2) = -2*c*lambda_1*(p12-p1)
%   F = -|gamma_12|*lambda_2*(p2-p12)/(dx/2) = -2*c*lambda_2*(p2-p12)
%
% where c:=(|gamma_12|/dx) and we impose continuity in flux and pressure 
%(cf. [1]). Substituting p12 = p2 + F/(c*lambda_2) we obtain
%
% F = 2*|gamma_12|/dx*(lambda_1*lambda_2)/(lambda_1+lambda_2)*(p1-p2)
%
%
% FIGURE: 
% |------|------|    |----------|----------|      |------|  |
% | p1  p12  p2 |    | lambda_1 | lambda_2 |      |      | dy
% |------|------|    |----------|----------|      |------|  |
%                                                  <-dx->
%    Pressure                 K/mu                 Spatial layout.
%    (cell-centered)       (cell-centered)         In 2D gamma_12 = dy 
%                                                  In 3D gamma_12 = dy*dz 
%
%
%
%[1] http://link.springer.com/article/10.1007%2Fs10596-007-9067-5
%[2] We make the code from Aarnes, Gimse and Lie for incompressible
%    single phase flow: "An introduction to the Numerivs of Flow in Porous
%    Media using Matlab"


% dimension of problem
dim = prm.dim;
ndim = numel(dim);
nvox = prod(dim);
h = prm.h;
% voxelvol = prod(h);

msg = ['This is ' mfilename ' using settings'];
disp(msg);
perfusion1c.printstructscreen(prm);

% make source/sink into a vector
F = Fmat(:);

%
% Compute transmissibilities by harmonic averaging
%

% K = K/mu
Kmat = Kmat/prm.mu;

% inverse
L = Kmat.^(-1);

% they can be represented as
%                       2|gamma_ij|
% t_ij = -------------------------------------------
%         delta x_i(lambda_i,ij^-1 + lambda_j,ij^-1)
%
% since delta x_i = delta x_j
%

% the term 2*|gamma|/h
tr = zeros(ndim,1);
trmat = cell(ndim,1);
for i = 1 : ndim
    v = dim;
    v(i) = v(i) + 1;
    trmat{i} = zeros(v);
    v = h;
    v(i) = [];
    tr(i) = 2*prod(v)/h(i);
end;
trmat{1}(2:dim(1),:,:) = tr(1)./(L(1:dim(1)-1,:,:,1) + L(2:dim(1),:,:,1));
trmat{2}(:,2:dim(2),:) = tr(2)./(L(:,1:dim(2)-1,:,2) + L(:,2:dim(2),:,2));
trmat{3}(:,:,2:dim(3)) = tr(3)./(L(:,:,1:dim(3)-1,3) + L(:,:,2:dim(3),3));

% Assemble TPFA discretization matrix
v = cell(6,1);
% x, above
v{1} = trmat{1}(1:dim(1),:,:);
% x, below
v{2} = trmat{1}(2:dim(1)+1,:,:);
% y, left
v{3} = trmat{2}(:,1:dim(2),:);
% y, right
v{4} = trmat{2}(:,2:dim(2)+1,:);
% z, above
v{5} = trmat{3}(:,:,1:dim(3));
% z, below
v{6} = trmat{3}(:,:,2:dim(3)+1);
sumv = zeros(nvox,1);
for i = 1 : numel(v)
    v{i} = v{i}(:);
    sumv = sumv + v{i};
end;

% value and position of bands in the matrix
vals = [-v{6}, -v{4}, -v{2}, sumv, -v{1}, -v{3}, -v{5}];
index = [-dim(1)*dim(2),-dim(1),-1,0,1,dim(1),dim(1)*dim(2)];

% make matrix
A = spdiags(vals,index,nvox,nvox);
A(1,1) = A(1,1) + sum(Kmat(:,1,1,1));

% solve linear system for the pressure p
%u = minres(A,F);

% b = voxelvol*Q = F;
b = F;

p = A\b;
pmat = reshape(p,dim(1),dim(2),dim(3));

% flow is gradient of pressure
% flowmat{1}(2:dim(1),:,:) = (pmat(1:dim(1)-1,:,:) - pmat(2:dim(1),:,:)).*trmat{1}(2:dim(1),:,:);
% flowmat{2}(:,2:dim(2),:) = (pmat(:,1:dim(2)-1,:) - pmat(:,2:dim(2),:)).*trmat{2}(:,2:dim(2),:);
% flowmat{3}(:,:,2:dim(3)) = (pmat(:,:,1:dim(3)-1) - pmat(:,:,2:dim(3))).*trmat{3}(:,:,2:dim(3));

%setup flowmat. Note that the values will be absolute, as the
%transmissibilities include the voxel-surface (see explanation above).
flowmat{1} = zeros(dim(1)+1,dim(2),dim(3));
flowmat{1}(2:end-1,:,:) = -(pmat(2:end,:,:) - pmat(1:end-1,:,:)).*trmat{1}(2:end-1,:,:);
flowmat{2} = zeros(dim(1),dim(2)+1,dim(3));
flowmat{2}(:,2:end-1,:) = -(pmat(:,2:end,:) - pmat(:,1:end-1,:)).*trmat{2}(:,2:end-1,:);
flowmat{3} = zeros(dim(1),dim(2),dim(3)+1);
flowmat{3}(:,:,2:end-1) = -(pmat(:,:,2:end) - pmat(:,:,1:end-1)).*trmat{3}(:,:,2:end-1);

