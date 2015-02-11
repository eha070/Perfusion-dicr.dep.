function [phimat] = porosity(Cmat,aifval,timeline)


dt = timeline([2:end end]) - timeline;
dt(end) = dt(end-1);
dt = (dt(:))';

dim = size(Cmat);
ntime = dim(4);
dim3 = dim(1:3);

% sum across time (NB scale by delta t!! Remember, this is the integral over time)
v = 0;
a = aifval;
for i = 1 : ntime-1
%    v = v + a(i)*dt(i);
    v = v + 0.5*(a(i+1)+a(i))*dt(i);
end;
sumCin = v;

% for voxels, sum across time
v = zeros(dim3);
for i = 1 : ntime-1
%    v = v + Cmat(:,:,:,i)*dt(i);
    v = v + 0.5*(Cmat(:,:,:,i+1)+Cmat(:,:,:,i))*dt(i);
end;
sumCv = v;


phiin = 1;
% compute for porosity
phimat = phiin * sumCv / sumCin;

% fix noisy pixels
phimat(phimat > 1) = 1;
phimat(phimat < 1e-2) = 1e-2;