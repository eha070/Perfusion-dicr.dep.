function [C,Cmask,mask] = preparedsc(im)
% PREPAREDSC  Prepares the DSC-MRI dataset for analysis, returning a mask
% of the brain as well as the GD concentration map (up to a multiplicative
% constant)
% Concentration map from "A control point interpolatio method for the
% non-parametric quantification of cerebral hemodynamics from dynamic
% susceptibility contrats MR", Mehdiratta 2013, Equation (7)
%

% Make a mask of the brain
[mask,mask4D] = createmask(im);

% Get GD concentration values
C = signal2conc(im,2,1);

% Mask out the background
Cmask = C.*mask4D;

%----------------------------------------

function [C] = signal2conc(im,TE,k)

im0 = mean(im(:,:,:,1:10),4);
th = 1e-1;
im0(im0 < th) = th;
C = im;

dim = size(im);
for i = 1 : dim(4)
    C(:,:,:,i) = -log(im(:,:,:,i)./im0)*k/TE;
end;

C(C < 0) = 0;
C(isinf(C)) = 0;
C(isnan(C)) = 0;

end

%----------------------------------------

function [mask,mask4D] = createmask(im)

%im = imcomplement(im);
d = mean(im,4);
mask = d > 450;
mask = bwkeep(mask,1,6);

mask4D = im;
dim = size(im);
for i = 1 : dim(4)
    mask4D(:,:,:,i) = mask;    
end;
end


end