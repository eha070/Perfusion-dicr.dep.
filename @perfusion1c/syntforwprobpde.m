function [Cmat] = syntforwprobpde(phimat,flowmat,Fmat,aifval,prm)
%function [Cmat] = syntforwprobpde(phimat,flowmat,Fmat,aifval,prm)
%
% INPUT:
%   phimat - Cell-centered porosity [fractions]
%  flowmat - Absolute flow (staggered) [mm^3/s]
%            Positive values of flowmat mean that flow goes from
%            - top   -> bottom
%            - left  -> right
%            - front -> back
%    Fmat  - inflow (cell-centered) [mm^3/s]
%  aifval  - arterial input function [mmol/mm^3]
%
%
% OUTPUT
% Cmat - simulated values [mmol/mm^3]

dim = size(phimat);
if numel(dim) == 2
    dim = [dim 1];
end;
h = prm.h;
cso = prm.cso;
csi = prm.csi;

msg = ['This is ' mfilename ' using settings'];
disp(msg);
perfusion1c.printstructscreen(prm);

% Define a timeline in seconds
dt = prm.dt;
timeline = prm.timeline;
ntime = numel(timeline);
dt = dt*ones(size(timeline));

maxaifval = max(aifval);
lim = [0,prm.scaling*maxaifval];

% matrix of tracer concentrations
cmat = zeros([dim,ntime]);
voxelvol = prod(h);

% neigh = [-1,0,0;
%     1,0,0;
%     0,-1,0;
%     0,1,0;
%     0,0,-1;
%     0,0,1];
% source/sink term
% sstermini = Fmat./phimat;
Fmatini = Fmat;
% sstermini = sstermini/voxelvol;

time = 0;
for i = 2 : ntime    
    time = time + dt(i);
   
    %setup source-inflow for current timepoint
    Fmat = Fmatini;
    for j = 1 : size(cso,1)
        Fmat(cso(j,1),cso(j,2),cso(j,3)) = Fmat(cso(j,1),cso(j,2),cso(j,3))*aifval(i-1);
    end;
    
    %setup source-outflow for current timepoint
    cmath = cmat(:,:,:,i-1);
    for j = 1 : size(csi,1)
        Fmat(csi(j,1),csi(j,2),csi(j,3)) = Fmat(csi(j,1),csi(j,2),csi(j,3))*cmath(csi(j,1),csi(j,2),csi(j,3));
    end;
    

    % sum of in and out fluxes
    sumflux = zeros(dim);
        
    %
    % Flow INTO voxels, i direction
    %
    
    % flow from top
    a   = flowmat{1}(1:end-1,:,:);
    ind = a > 0;
    cj  = perfusion1c.transim(cmath,-1,0,0);
    sumflux(ind) = sumflux(ind) + abs(a(ind)).*cj(ind);
    % flow from bottom
    a   = flowmat{1}(2:end,:,:);
    ind = a < 0;
    cj  = perfusion1c.transim(cmath,1,0,0);
    sumflux(ind) = sumflux(ind) + abs(a(ind)).*cj(ind);

    %
    % Flow INTO voxels, j direction
    %

    % flow from left
    a   = flowmat{2}(:,1:end-1,:);
    ind = a > 0;
    cj  = perfusion1c.transim(cmath,0,-1,0);
    sumflux(ind) = sumflux(ind) + abs(a(ind)).*cj(ind);
    % flow from right
    a   = flowmat{2}(:,2:end,:);
    ind = a < 0;
    cj  = perfusion1c.transim(cmath,0,1,0);
    sumflux(ind) = sumflux(ind) + abs(a(ind)).*cj(ind);

    %
    % Flow OUT from voxels, i direction
    %
    
    % flow towards top
    a   = flowmat{1}(1:end-1,:,:);
    ind = a < 0;
    sumflux(ind) = sumflux(ind) - abs(a(ind)).*cmath(ind);
    % flow towards bottom
    a   = flowmat{1}(2:end,:,:);
    ind = a > 0;
    sumflux(ind) = sumflux(ind) - abs(a(ind)).*cmath(ind);

    %
    % Flow OUT from voxels, j direction
    %

    % flow towards left
    a   = flowmat{2}(:,1:end-1,:);
    ind = a < 0;
    sumflux(ind) = sumflux(ind) - abs(a(ind)).*cmath(ind);
    % flow towards right
    a   = flowmat{2}(:,2:end,:);
    ind = a > 0;
    sumflux(ind) = sumflux(ind) - abs(a(ind)).*cmath(ind);

    % sum flows results from flux-fields
    fluxterm = Fmat + sumflux;

    % adding up.
    % remember: dC/dt = fluxterm*appropriateUnits
    % where the fluxterm has units mmol/s
    del = dt(i)*fluxterm./(phimat*voxelvol); 
    a = cmat(:,:,:,i-1) + del;
    if nnz(a<0)
        warning('Some concentrations are below 0, going to keyboard.');
        keyboard
    end
    a(a < 0) = 0;
    
    cmat(:,:,:,i) = a;
    
    % plotting
    if round(i/500) == i/500        
        ci = cmat(:,:,:,i);
        figure(1);
            imagesc(ci,lim);
            colormap(gray); 
            axis image; 
            caxis(lim)
            drawnow
    end;
    
end;
Cmat = zeros(size(cmat));
for i = 1 : ntime
    Cmat(:,:,:,i) = cmat(:,:,:,i).*phimat;
end;




