function [lenim] = arclength(q,Q,h)

ind = ne(Q,0);
dim = size(Q);
ndim = numel(dim);
if numel(dim) == 2
    dim = [dim,1];
end;

% transfer to cell centered grid
q1 = (q{1}(2:end,:,:) + q{1}(1:end-1,:,:))/2;
q2 = (q{2}(:,2:end,:) + q{2}(:,1:end-1,:))/2;
q3 = (q{3}(:,:,2:end) + q{3}(:,:,1:end-1))/2;
q{1} = q1;q{2} = q2;q{3} = q3;

a = perfusion1c.distimage(ind,4);
a = a == 1;
a = bwperim(a);
a = perfusion1c.constborder(a,1,0);
lenim = inf(dim);

indseed = find(a);
for i = 1 : numel(indseed)
    
    % this pixel
    [c1,c2,c3] = ind2sub(dim(1:3),indseed(i));
   
    
    clear indtraj len;
    mode = [1,-1];
    for j = 1 : 2
        pos = [c1,c2,c3];

        ctraj{j} = [];
        indtraj{j} = [];
        len{j} = 0;
        while 1

            % go one direction
            vect(1) = perfusion1c.interp3c(q{1},pos);
            vect(2) = perfusion1c.interp3c(q{2},pos);
            vect(3) = perfusion1c.interp3c(q{3},pos);

            % normalize
            absvect = norm(vect);
            vect = vect/absvect;

            % move
            dt = 0.30;
            deltapos = dt*vect;
            pos = pos + mode(j)*deltapos;
            for k = 1 : 3
                if pos(k) > dim(k)
                    pos(k) = dim(k);
                end;
                if pos(k) < 1
                    pos(k) = 1;
                end;
            end;

            posr = round(pos);
            ctraj{j} = [ctraj{j}; posr];

            % scale to real distance
            deltapos = deltapos .* h;

            % absolute value
            absdeltapos = norm(deltapos);

            % add to the integral
            len{j} = len{j} + absdeltapos;
            if ne(Q(posr(1),posr(2),pos(3)),0)
                break;
            end;

        end;
        
    end;
    lentot = len{1} + len{2};
    indtraj{1} = sub2ind(dim,ctraj{1}(:,1),ctraj{1}(:,2));
    indtraj{2} = sub2ind(dim,ctraj{2}(:,1),ctraj{2}(:,2));
    indtraj = [indtraj{1};indtraj{2}];
    indtraj = unique(indtraj);
    val = lenim(indtraj);
    val(val > lentot) = lentot;
    lenim(indtraj) = val;
    
    % show(lenim,1)
%    pause
end;

% remove isinf values
val1 = lenim([2:end end],:);val1 = val1(:);
val2 = lenim([1 1:end-1],:);val2 = val2(:);
val3 = lenim(:,[2:end end]);val3 = val3(:);
val4 = lenim(:,[1 1:end-1]);val4 = val4(:);
val = [val1,val2,val3,val4];
ind = isinf(lenim);
val(isinf(val)) = nan;
val = nanmean(val,2);
lenim2 = reshape(val,dim(1),dim(2));
lenim(ind) = lenim2(ind);
% there may be some more nan-values left if all neighbors are nan
lenim(isnan(lenim)) = max(lenim(:));
%lenim(end,1) = max(lenim(:));
%lenim(1,end) = max(lenim(:));

% find shortest distance between source and sink to account for numerical
% errors in curve length measurements
ind1 = find(Q > 0);
ind2 = find(Q < 0);
clear c1 c2;
[c1(:,1), c1(:,2),c1(:,3)] = ind2sub(dim,ind1);
[c2(:,1), c2(:,2),c2(:,3)] = ind2sub(dim,ind2);
minlen = inf;
for i = 1 : numel(ind1)
    v = c2 - repmat(c1(i,:),size(c2,1),1);
    v = v.*repmat(h,size(v,1),1);
    v = v.^2;
    v = sum(v,2);
    v = sqrt(v);
    v = min(v);
    minlen = min(minlen,v);
end;
delta = minlen - min(lenim(:));
lenim = lenim + delta;

% to the real boundary
lenim = lenim + h(1);


