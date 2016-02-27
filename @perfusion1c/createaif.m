function [aifval] = createaif(aiftype,timeline)
%function [aifval] = createaif(aiftype,timeline)
%
% INPUT:
% aiftype  - 'gamma', 'parker', 'delta'
% timeline - timeline in seconds
%


% Define an AIF curve
if isequal(aiftype,'parker')
    aifval = perfusion1c.getParkerAIF(timeline/60);
elseif isequal(aiftype,'gamma')
    aifval = perfusion1c.getGammaAIF(timeline/60);
elseif isequal(aiftype,'delta')
    aifval = zeros(ntime,1);
    ind = find(timeline > 5);
    ind = ind(1);
    aifval(ind) = 1e6;
end;

% The Parker AIF is 1mMol/liter = 1mMol/1000ml = 1e-3mMol/ml = 1e-6mMol/mm^3
aifval = aifval*1e-6;
aifval = (aifval(:))';

% aifval = delta(timeline-prm.aifpeak,prm.aifwidth);
% aifval = scale(aifval,0,1); 
plot(timeline,aifval,'-')
pause(2)

