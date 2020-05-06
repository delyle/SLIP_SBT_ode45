function [ B ] = wrtsort( A, dim, i, mode ) 
%WRTSORT Sorts with respect to a single row or column.
%   B = Sorted array
%   A = Array to sort
%   dim = dimension along which to sort (1 = along columns (default),
%           2 = along rows)
%   i = row/column to sort with respect to (default is 1)
%   mode = 'ascend' <- Sort in ascending order (default)
%        = 'descend' <- Sort in descending order

narginchk(1,4)
if nargin == 1
    dim = 1;
    i = 1;
    mode = 'ascend';
elseif nargin < 3
    i = 1;
    mode = 'ascend';
elseif nargin ==3 
    if isstrprop(i,'alphanum')
        mode = i;
        i = 1;
    else
        mode = 'ascend';
    end
elseif isstrprop(i,'alpha')
    mode1 = i;
    i = mode;
    mode = mode1;
end

if ~(mod(dim,1)==0)
    txt = 'Sort dimension must be an integer (1 or 2)';
    error(txt)
end

s = size(A);
B = NaN(s);
if dim == 1
    [B(:,i), I] = sort(A(:,i),mode);
    B(:,[1:i-1,i+1:s(2)]) = A(I,[1:i-1,i+1:s(2)]);
elseif dim == 2
    [B(i,:), I] = sort(A(i,:),mode);
    B([1:i-1,i+1:s(1)],:) = A([1:i-1,i+1:s(1)],I);
elseif dim > 2
    txt = ['This function currently cannot handle arrays with dimension >2',char(10),...
        'Future functionality may be added'];
    error(txt)
end
end





