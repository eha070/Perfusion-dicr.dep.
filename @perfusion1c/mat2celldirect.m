% MAT2CELLDIRECT Convert a matrix into a cell matrix
% A = CELL2MATDIRECT(A) Converts the matrix A into a cel matrix which is
% returned
%
function [A] = mat2celldirect(a)

% if already a cell then dont do anything
if iscell(a)
    A = a;
    return;
end;

[M N O] = size(a);
A = cell(M,N,O);
for i = 1 : M
    for j = 1 : N
        for k = 1 : O
            A{i,j,k} = a(i,j,k);
        end;
    end;
end;
