function [P,Pnorm] = flux2perf(X,lenim,Fmat,h)
% Convert flux values to perfusion values. Assuming the input X is the flux 
% with units mm^3/s/mm^2
%

% transfer flux to cell centered
X1 = (X{1}(2:end,:,:) + X{1}(1:end-1,:,:))/2;
X2 = (X{2}(:,2:end,:) + X{2}(:,1:end-1,:))/2;
X3 = (X{3}(:,:,2:end) + X{3}(:,:,1:end-1))/2;

X1 = X1/(h(2)*h(3));
X2 = X2/(h(1)*h(3));
X3 = X3/(h(1)*h(2));
absq = sqrt(X1.^2 + X2.^2 + X3.^2);

% fix the source/sinc later
absq(ne(Fmat,0)) = 0;

% perfusion in units mm^3/s/mm^3
P = absq./lenim;

% fix the source and sink values
ind = Fmat > 0;
len = lenim(ind);
% NB this area is an assumption of source along upper row, not valid
% elsewhere!!
val = abs(Fmat(ind))/(h(2)*h(3));
val = val./len;
P(ind) = val;

ind = Fmat < 0;
len = lenim(ind);
% NB this area is an assumption of sink along lower row, not valid
% elsewhere!!
val = abs(Fmat(ind))/(h(2)*h(3));
val = val./len;
P(ind) = val;

% NB: to recover the total perfusion use
% voxelvol = prod(h);
% totper = sum(P(:))*voxelvol

% Convert to mm^3_f/min/mm^3_v
Pnorm = P;
Pnorm = Pnorm*60;

% convert to mm^3_f/min/100mm^3_v
Pnorm = Pnorm*100;

return;

show(X{1},1);colorbar
show(X{2},2);colorbar

dim = size(X{1});
mid = dim(1)/2;
mid = round(mid);
a = X{1}(mid,:);
sum(a)
pause

% find perfusion as a sum of all positive (!) perfusion across the vertices
a1 = X{1}(1:end-1,:,:);
a2 = X{1}(2:end,:,:);
a3 = X{2}(:,1:end-1,:);
a4 = X{2}(:,2:end,:);
a5 = X{3}(:,:,1:end-1);
a6 = X{3}(:,:,2:end);

ind1 = a1 > 0;
ind2 = a2 < 0;
ind3 = a3 > 0;
ind4 = a4 < 0;
ind5 = a5 > 0;
ind6 = a6 < 0;

% total positive flow into the voxel: mm^3/s/voxel
P1 = abs(a1) .* ind1 + ...
     abs(a2) .* ind2 + ...
     abs(a3) .* ind3 + ...
     abs(a4) .* ind4 + ...
     abs(a5) .* ind5 + ...
     abs(a6) .* ind6;
 
% total negative flow into the voxel: mm^3/s/voxel
P2 = abs(a1) .* imcomplement(ind1) + ...
     abs(a2) .* imcomplement(ind2) + ...
     abs(a3) .* imcomplement(ind3) + ...
     abs(a4) .* imcomplement(ind4) + ...
     abs(a5) .* imcomplement(ind5) + ...
     abs(a6) .* imcomplement(ind6) ;

% This tricks fixes the problems at the source and sink of zero flux for
% either of them
% This quantity has units mm_f^3/s
P = (P1 + P2)/2;
dim = size(P);
mid = dim(1)/2;
a = P(mid,:);
sum(a)
pause

% convert mm_f^3/s -> mm_f^3/s/mm_v^3. This is the same as ml_f/s/ml_v
voxelvol = prod(h);
Pnorm = P/voxelvol;
mean(Pnorm(:))
mean(P(:))
show(Pnorm,3)
show(P,4)
pause
% Convert to ml_f/min/ml_v
Pnorm = Pnorm*60;

% convert to ml_f/min/100ml_v
Pnorm = Pnorm*100;


