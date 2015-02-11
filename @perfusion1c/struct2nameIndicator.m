function basename = struct2nameIndicator(prm,varargin)
%function basename = struct2nameIndicator(prm,varargin)
% Creates a basename from the fields of prm. 
% 
% INPUT: 
%   prm - a struct array
% 
% VARARGIN:
%  field1,...,fieldN - fields to create the basename. class(field)=char.
%   
% OUPUT:
%    basename - Has the form field1-value1-field2-value2-...-fieldN-valueN
%                                          (c)Constantin Heck, 30-Jan-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 


assert(isstruct(prm),'prm needs to be struct');

basename = '';
for i = 1:numel(varargin)
    field = varargin{i};
    assert(ischar(field),'fields need to be char');
    value = eval(['prm.',field]);
    
    %convert fieldname to identifier in fielname
    switch field
        case 'stepred'
            field = 'red'; 
        case 'aiftype'
            field = 'aif';
        case 'Kopt'
            field = 'K';
        case 'dim'
            value = value(1);
        case 'phiopt'
            field = 'phi';
    end
    
    
    basename = [basename,field,'-',num2str(value),'-'];
end

%remove last '-' sign
basename = basename(1:end-1);

end

