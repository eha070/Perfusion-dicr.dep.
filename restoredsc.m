function [CBF,CBFn,I] = restoredsc(Cim,brainmask,hdr,aifmask)

roi = 'voxel';
% roi = 'average';

% Thrash first three volumes due to magnetization issues
Cim = Cim(:,:,:,4:end);
dim = size(Cim);
ntime = dim(4);


% show4D(Cim)

% % smooth the data
% filt = fspecial('gaussian',5,3);
% for i = 1 : ntime
%     Cim(:,:,:,i) = imfilter(Cim(:,:,:,i),filt,'replicate');
% end;

% Extract AIF curve
aifval = im2vec4D(Cim,aifmask);
aifval = mean(aifval,1);


%get deconvolution matrix
deltaT  = hdr.dt;
A       = perfusion1c.getCircularConvolutionMatrix(aifval,deltaT);
% A       = perfusion1c.getLinearConvolutionMatrix(aifval,deltaT);
fprintf('Starting SVD...');
tic; [U,S,V] = svd(A);
fprintf('...done. Elapsed time: %1.3fs\n',toc);


% plot(imval,'r');
% hold on;
% plot(aifval,'b');
% hold off

% Make timeline
timeline = (0:deltaT:(ntime-1)*deltaT)';

% oscind = 0.01:0.005:0.1;
oscind = 0.01;

% scaling = [1,2,4];
scaling = [1];
CBF = cell(numel(scaling),1);
for i = 1 : numel(scaling)
                
    for j = 1 : dim(4)        
        a = imresize3d(Cim(:,:,:,j),scaling(i),'bilinear');
        dimhere = [size(a),ntime];
        if j == 1
            imhere = zeros(dimhere);
        end;
        imhere(:,:,:,j) = a;
    end
    brainmaskhere = imresize3d(brainmask,scaling(i),'nearest');

    % Extract average intensity value within mask
    imval = im2vec4D(imhere,brainmaskhere);
    if strcmp(roi,'average')
        imval = nanmean(imval,1);    
    end;
    nvox = size(imval,1);

    
    for j = 1 : nvox

        msg = ['Voxel ' int2str(j) ' out of ' int2str(nvox)];
        disp(msg);
    
        for k = 1 : numel(oscind)
            
            % oscilattion index
            OI = oscind(k);            

            msg = ['Osciallation index ' num2str(OI)];
            disp(msg);
            
            val = imval(j,:);
            val = val(:);
        
%             plot(timeline,val);
%             hold on;
%             plot(timeline,aifval)
%             hold off
%             pause

            % [CBF,I,C] = perfusion1c.linearDeconvolution(imval,timeline,OI,U,S,V);
            [CBF{i}(j,k),I,C] = perfusion1c.circularDeconvolution(val,timeline,OI,U,S,V);         
            
            perf = CBF{i}(j,k)*60*100;
            msg = ['Perfusion: ' num2str(perf) ' ml/min/100ml'];
            disp(msg);
            % mean(CBF{i})*60*100
%             plot(I)
%             pause
        end;
    end;
end;
% CBF{1}*60*100
% pause

% %do the maximum slope technique
% CBFrecMS = perfusion1c.maximumSlope(imval,timeline,aifval);

for i = 1 : numel(CBF)
    CBFn{i} = CBF{i}*60*100;    
end;
% CBFrecMS*60*100
% plot(oscind,CBFn)
% plot(Irec)

% %get CBV
% CBVrec = perfusion1c.cbvEstimation(Cav,timelinelow,AIFlow);

function [data] = im2vec4D(varargin)

f = varargin{1};
dim = size(f);
dim3 = dim(1:3);
if nargin == 1
    roi = ones(dim3);
else
    roi = varargin{2};
end;

nvox = sum(roi(:));

% rearrange the data into vectors n x time
ntime = size(f,4);
data = zeros(nvox,ntime);

for i = 1 : ntime    
    datahere = f(:,:,:,i);
    datahere = datahere(roi == 1);
    data(:,i) = datahere(:);
end;

