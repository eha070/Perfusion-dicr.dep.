function [] = createindicatorconv(prmin)
% This function uses estimated flow from synthetic model to generate a
% syntetic concentration map, to see if we can restore that one by
% deconvolution


% settings
[prm,Fmat] = settings;
prm = perfusion1c.mergestruct(prm,prmin);

stepred = prm.stepred;

% basename flow
basenameflow = perfusion1c.providenameflow(prm.phiopt,prm.Kopt,prm.dim);

% load flow and pressure field
pathload = ['results/synt-createflowTPFA-' basenameflow '.mat'];
msg = ['Loading ' pathload];
disp(msg);
D = load(pathload);

perfmat = D.perfmat;
phimat = D.phimat;

% Reduced data: aif curve, timeline and the concentration map to find TTP
basenameindicator = perfusion1c.providenameindicator(prm.phiopt,prm.Kopt,prm.dim,prm.aiftype,prm.T);
pathload = ['results/synt-createindicatorpde-' basenameindicator '.mat'];
msg = ['Loading ' pathload];
disp(msg);
E = load(pathload);
aifval = E.aifval;
timeline = E.prm.timeline;
% The loaded timeline from prm struct is wrong since we reduced it, 
% get the right one here for saving
prm.timeline = timeline;

% mean transit time, eq (18) in NMR Biomed 2013 Sourbron
T = phimat./perfmat;
T(isinf(T)) = 0;

dt = timeline([2:end end]) - timeline;
dt(end) = dt(end-1);
ntime = numel(timeline);

dim3 = size(perfmat);
if numel(dim3) < 3
    dim3 = [dim3,1];
end;
dim = [dim3,ntime];

% find time to peak map
[tilde,delay] = max(E.Cmat,[],4);
[tilde,delayaif] = max(E.aifval);
% account for the delay in the true concentration map field as wel
delay = delay - delayaif;
% delay = zeros(dim3);

clear E;

% move to zero (one, since matlab starts counting at 1) delay
% delay = delay - min(delay(:)) + 1;

aifmat = zeros(dim);
for i = 1 : dim(1)
    for j = 1 : dim(2)
        for k = 1 : dim(3)
            aifmat(i,j,k,:) = aifval;
            v = delay(i,j,k);
            n = ntime-v;
            aifmat(i,j,k,v+1:end) = aifmat(i,j,k,1:n);            
            aifmat(i,j,k,1:v) = 0;
        end;
    end;
end;

% make delay field on real time seconds to show
for i = 1 : dim(1)
    for j = 1 : dim(2)
        for k = 1 : dim(3)
            delay(i,j,k) = timeline(delay(i,j,k)+1);
        end;
    end;
end;

%method = 'summation';
method = 'conv';
if isequal(method,'conv')
    % forward model from convolution
    dim = [dim3 ntime];
    Cmat = zeros(dim);
    for i = 1 : dim(1)
        i
        for j = 1 : dim(2)
            for k = 1 : dim(3)
                val = squeeze(aifmat(i,j,k,:));
                IH = exp(-timeline./T(i,j,k));
                C = conv(IH,val);
                C = perfmat(i,j,k)*dt(1)*C(1:ntime);       
                Cmat(i,j,k,:) = C;
            end;
        end;
    end;
elseif isequal(method,'summation')
    % This is a slow version, but more transparent and give (almost) same answer
    dim = [dim3 ntime];
    Cmat = zeros(dim); 
    for i = 1 : ntime
        if i/100 == round(i/100)
            [i,ntime]    
        end;
        t = timeline(i);
        a = zeros(dim3);
        for j = 1 : i
            s = timeline(j);
            f = exp(-(t-s)./T);
            g = aifmat(:,:,:,j);
            a = a + f.*g*dt(j);
        end;
        Cmat(:,:,:,i) = perfmat.*a;
    end;
end;


savpaper = 1;
savdata = 1;
basenameindicator = perfusion1c.providenameindicator(prm.phiopt,prm.Kopt,prm.dim,prm.aiftype,prm.T);
lim = [0,prm.scaling*max(aifval)];

if savpaper
    
    % make panel
    start = 0;
    starti = find(timeline > start,1);
    n = round(numel(timeline)/2) - starti;
    step = floor(n/16);
    c = starti;
    timeplot = zeros(16,1);
    C = 0;
    for i = 1 : 4
        for j = 1 : 4            
            C = C + 1;
            A{i,j} = Cmat(:,:,:,c);
            clim{i,j} = lim*phimat(1,1,1);
            timeplot(C) = prm.timeline(c);
            c = c + step;
        end;
    end;
    
    H = perfusion1c.panelstruct(A,0.01,700,'lim',clim);
    pathsave = ['figs' '/' 'syntconv-' mfilename '-' basenameindicator '-panel.eps'];
    msg = ['Saving ' pathsave];
    disp(msg);
    print(H,pathsave,'-deps')
    
    % makemovieraw(numberim,1,1,pathsave);
end

if savdata
    pathsave = ['results/syntconv-' mfilename '-' basenameindicator  '.mat'];
    msg = ['Saving ' pathsave];
    disp(msg);
    save(pathsave,'aifval','delay','Cmat','prm','-v7.3')
end;



