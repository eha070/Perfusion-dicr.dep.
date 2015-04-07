function [Cmat] = syntforwprobpde(phimat,flowmat,Fmat,aifval,prm)

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
   
    cmath = cmat(:,:,:,i-1);
    Fmat = Fmatini;
    
    for j = 1 : size(cso,1)
        Fmat(cso(j,1),cso(j,2),cso(j,3)) = Fmat(cso(j,1),cso(j,2),cso(j,3))*aifval(i-1);
    end;
    
    for j = 1 : size(csi,1)
        Fmat(csi(j,1),csi(j,2),csi(j,3)) = Fmat(csi(j,1),csi(j,2),csi(j,3))*cmath(csi(j,1),csi(j,2),csi(j,3));
    end;
    

    % sum of in and out fluxes
    sumflux = zeros(dim);
        
    %
    % Flow into voxel x direction
    %

    % the face area
    % dA = h(2)*h(3);
    
    % NB the flux flowmat is already multiplied by the voxel area so
    % therefore we set dA = 1 always
%     dA = h(2)*h(3);
    dA = 1;

    % flow from top
    a = flowmat{1}(1:end-1,:,:);
    ind = a > 0;
    % tracer comes from top
    cj = perfusion1c.transim(cmath,-1,0,0);
    sumflux(ind) = sumflux(ind) + abs(a(ind)).*cj(ind)*dA;
    % flow from bottom
    a = flowmat{1}(2:end,:,:);
    ind = a < 0;
    cj = perfusion1c.transim(cmath,1,0,0);
    sumflux(ind) = sumflux(ind) + abs(a(ind)).*cj(ind)*dA;

    %
    % Flow into voxel y direction
    %

    % the face area
%     dA = h(1)*h(3);
    dA = 1; %CONSTANTIN
    
    % flow from left
    a = flowmat{2}(:,1:end-1,:);
    ind = a > 0;
    cj = perfusion1c.transim(cmath,0,-1,0);
    sumflux(ind) = sumflux(ind) + abs(a(ind)).*cj(ind)*dA;
    % flow from bottom
    a = flowmat{2}(:,2:end,:);
    ind = a < 0;
    cj = perfusion1c.transim(cmath,0,1,0);
    sumflux(ind) = sumflux(ind) + abs(a(ind)).*cj(ind)*dA;

    %
    % Flow out from voxel x direction
    %
    cj = cmath;
    
    % the face area
%     dA = h(2)*h(3);
    dA = 1; %CONSTANTIN
    
    % flow towards top
    a = flowmat{1}(1:end-1,:,:);
    ind = a < 0;
    sumflux(ind) = sumflux(ind) - abs(a(ind)).*cj(ind)*dA;
    % flow towards bottom
    a = flowmat{1}(2:end,:,:);
    ind = a > 0;
    sumflux(ind) = sumflux(ind) - abs(a(ind)).*cj(ind)*dA;

    %
    % Flow out from voxel y direction
    %

    % the face area
%     dA = h(2)*h(3);
    dA = 1; %CONSTANTIN



    % flow towards left
    a = flowmat{2}(:,1:end-1,:);
    ind = a < 0;
    sumflux(ind) = sumflux(ind) - abs(a(ind)).*cj(ind)*dA;
    % flow towards right
    a = flowmat{2}(:,2:end,:);
    ind = a > 0;
    sumflux(ind) = sumflux(ind) - abs(a(ind)).*cj(ind)*dA;

    % sum of fluxes
    fluxterm = Fmat + sumflux;

    % adding up
    del = dt(i)*fluxterm./(phimat*voxelvol);

    a = cmat(:,:,:,i-1) + del;
    a(a < 0) = 0;
    cmat(:,:,:,i) = a;
    
    % to plot
    Chere = cmat(:,:,:,i);
    
    
    if round(i/500) == i/500        
        figure(1);imagesc(Chere,lim);colormap(gray);axis image;drawnow
    end;
    
end;
Cmat = zeros(size(cmat));
for i = 1 : ntime
    Cmat(:,:,:,i) = cmat(:,:,:,i).*phimat;
end;




