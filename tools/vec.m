function [ x ] = vec( X )
% ----------------------------------------------------------------------- %
% Vectorizes a matrix to a column vector.
%
% T. Reynolds -- RAIN Lab
% ----------------------------------------------------------------------- %

[m,n]   = size(X);
x       = reshape(X,m*n,1);

end

