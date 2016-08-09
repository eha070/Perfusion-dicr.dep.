function [F,phi,Irec,Crec,para] = fastPerfusionAnalysis(C,aif,timeline,thres,varargin) 
%function [F,phi,Irec,Crec,para] = fastPerfusionAnalysis(C,aif,timeline,thres,varargin) 
% Does CBV and CBF estimation in an efficient mode with a global threshold.
%
% INPUT:
%   C        - Tissue curves, size(C)=[ncurve,k];
%   aif      - The arterial input function
%   timeline - REGULARLY SPACED timeline in seconds.
%   thres    - threshold for the deconvolution (fraction with respect to
%              the larges eigenvalue)
%
% VARARGIN:
%   mode - Deconvolution Mode. Allowed is
%          bcSVD - Block-Circular SVD
%          sSVD  - standard deconvolution
% 
% OUTPUT:
%   F    - the estimated bloodflow in ml/s/ml
%   phi  - the estimated bloodvolume as a fraction
%   Irec - recovered residue functions
%   Crec - Estimated approximation to C
% 
%                                          (c)Constantin Heck, 09-Aug-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

%default parameters
mode = 'bcSVD';

%overwrites default parameter 
for j=1:2:length(varargin),  
	eval([varargin{j},'=varargin{',int2str(j+1),'};']); 
end 
 

%% data preprocessing

C        = reshape(C,[],numel(timeline));
ncurve   = size(C,1);
aif      = aif(:);
timeline = timeline(:);
k        = numel(timeline);
deltaT   = timeline(2)-timeline(1);




%% deconvolution analysis

fprintf('-----------------------------------------------\n\n');
fprintf('\t\tDeconvolution analysis: k=%i, ncurve=%i\n\n',k,ncurve);
fprintf('-----------------------------------------------\n\n\n');

%get deconvolution matrix
switch mode
    case 'bcSVD'
        M = perfusion1c.getCircularConvolutionMatrix(aif,deltaT);
    case 'sSVD'
        M = perfusion1c.getLinearConvolutionMatrix(aif,deltaT);
    otherwise
        error('mode is unknown');
end
     

%do the svd
fprintf('Starting SVD for a %i by %i matrix...',size(M,1),size(M,2)); tic;
    [U,S,V] = svd(M);
fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);


%prepare inverse of S
d       = diag(S);
dinv    = 1./d; 

%cut-off all eigenvalues which are too large
thresEig  = thres*max(d);
idx       = (d<=thresEig);

%get Sinv
dinv(idx) = 0;
Dinv      = diag(dinv,0);

%get cutOffS
d(idx)    = 0;
D         = diag(d,0);

%setup regularized M and its inverse
MregInv = V*Dinv*U';
Mreg    = U*D*V';

%zero-padded C
switch mode
    case 'bcSVD'
        fprintf('Zero-padding C...'); tic;
        CHat     = [C,zeros(ncurve,k)];
        fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);
    case 'sSVD'
        CHat     = C;
end

%recover Irec and Crec
fprintf('Deconvolving...');
tic;
Irec  = CHat*MregInv';
fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);

%get actual F
F    = max(Irec,[],2);
F    = F(:);


if nargout > 2
    
    fprintf('Reconstructing...');
    tic;
    Crec  = Irec*Mreg';
    fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);
    
    %reconstruction for the averaged curve
    fprintf('Reconstruction of averaged curve...');    
    tic;
    Cmean  = mean(CHat,1);
    Imean  = Cmean*MregInv';
    CmeanR = Imean*Mreg';
    Fmean  = max(Imean,[],2);
    para = struct('Cmean',Cmean(:),'Imean',Imean(:),'CmeanR',CmeanR(:),'Fmean',Fmean,'meanF',mean(F(:)));
    fprintf('...done. Elapsed time: %1.3fs.\n\n',toc);
    
end 



 




%% get CBV

fprintf('-----------------------------------------------\n\n');
fprintf('\t\tCBV-Estimation: k=%i, ncurve=%i\n\n',k,ncurve);
fprintf('-----------------------------------------------\n\n\n');


%get vector for integration
e  = timeline(2:end)-timeline(1:end-1);
Av = spdiags(1/2*ones(k,2),[0,1],k-1,k);
e  = e(:)'*Av; e=e(:);

%calculate phi
fprintf('Estimating phi...');
tic;
phi = (C*e)./(e'*aif(:));
fprintf('...finished. Elapsed time: %1.3fs.\n\n',toc);


%get phi
phi = phi(:);


 
end 
