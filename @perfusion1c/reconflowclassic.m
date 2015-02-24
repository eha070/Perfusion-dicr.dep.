function [qmat,perfmat,perfmatn] = reconflowclassic(im,timeline,organroi,aifval,prm)


dim = size(im);
dim3 = dim(1:3);
ntime = dim(4);
nvox = nnz(organroi);

dt = timeline(2:end) - timeline(1:end-1);
dt = dt(1);

% make sure its a row vector
aifval = (aifval(:))';

% indices of organ roi
indorganroi = organroi > 0;

% % image values for AIF
% imaifvec = im2vec4D(im,aifroi);
% imaifvec = imaifvec/phiin;

% image values within roi
imvec = perfusion1c.im2vec4D(im,organroi);
% imvec = imvec*1e3;
% aifval = aifval*1e3;

if isequal(prm.reconflowdeconvmethod,'ls')
    A = getlsmatrix(aifval,timeline);
    A = A + 0.0001*eye(size(A));
end;

aifval = aifval(:);

if isequal(prm.reconflowdeconvmethod,'circSVD')
    
    %get linear convolution matrix
    A = perfusion1c.getCircularConvolutionMatrix(aif,dt);

    % do the SVD
    [U,S,V] = svd(A);
elseif isequal(prm.reconflowdeconvmethod,'linearSVD')
    
    %get linear convolution matrix
    A = perfusion1c.getLinearConvolutionMatrix(aifval,dt);

    % do the SVD
    [U,S,V] = svd(A);
    
    
end;

msg = ['Using deconvolution method ' prm.reconflowdeconvmethod];
disp(msg);

% Deconvolution
perfvec = zeros(nvox,1);
for i = 1 : nvox

    % concentration curve of this voxel
    C = imvec(i,:);
    C = (C(:));

    % deconvolution
    if isequal(prm.reconflowdeconvmethod,'lucy')        
        I = deconvlucy(C,aifval,[],0.1);
        I = I/dt;
        perfvec(i) = max(I);
    elseif isequal(prm.reconflowdeconvmethod,'wiener')
        I = deconvwnr(C,aifval,0.01);
        I = I/dt;
        perfvec(i) = max(I);
    elseif isequal(prm.reconflowdeconvmethod,'ls')
        I = lsqlin(A,C,[],[]);
        perfvec(i) = max(I);
    elseif isequal(prm.reconflowdeconvmethod,'circSVD')        
        [F,~,~] = perfusion1c.circularDeconvolution(C,timeline,prm.reconflowIO,U,S,V);
        perfvec(i) = F;
    elseif isequal(prm.reconflowdeconvmethod,'linearSVD')        
        [F,~,~] = perfusion1c.linearDeconvolution(C,timeline,prm.reconflowIO,U,S,V);
        perfvec(i) = F;        
    end;
    
%     if i/500 == round(i/500)
%         [i,nvox]
%     end;
    
end;

% reshape as image
perfmat = zeros(dim3);
perfmat(indorganroi) = perfvec;

%
% Perfmat comes out as mm^3/sec/mm^3, convert to ml/min/100ml
%

% Convert to ml/min/ml
perfmatn = perfmat*60;

% convert to ml/min/100ml
perfmatn = perfmatn*100;

% dont know this quantity by this model
qmat = [];

%--------------------------------------------


function [A] = getlsmatrix(aifval,timeline)

dt = timeline(2)-timeline(1);

% make lower triagonal matrix A
ntime = numel(aifval);
A = zeros(ntime,ntime);
for i = 1 : ntime
    % deltat = t_i - t_j
    deltat = timeline(i) - timeline(1:i);
    
    % find the aif value at deltat
    val = interp1(timeline,aifval,deltat);
    
    % put into matrix
    A(i,1:i) = val;
end;
A = A*dt;

%----------------------------

% 
% function [I] = deconvls(C,V,Sinv,Ut)
% 
% % make column vector
% C = (rowvect(C))';
% 
% 
% % the solution to the linear system from SVD
% I = V*Sinv*Ut*C;
% 
% 
% %---------------------------------------------

