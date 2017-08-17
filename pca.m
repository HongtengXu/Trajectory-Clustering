function [U,lambda,xc] = pca(X,L)
% pca   principal component analysis for input data set.
% 
% Syntax:
%        [U,lambda,xc] = pca(X,L)
%       
% Input:
%        X = a NxD matrix which contains one D dimension data in each
%            row, or in other words X =[x1';x2';...xn'].
%        L = the dimension to which the dataset will be reduced.
% 
% Output:
%        U = a DxL matrix which stores principle components in columns.
%   lambda = a vector containing min(N,D) eigenvalues.
%       xc = the centroid of input data set.
% 
% History:
%    Weiran Wang        October 19, 2008   created.
%    Weiran Wang        January 26, 2009   revised.
%    Weiran Wang        February 24, 2009  revised.
%    Weiran Wang        April 1, 2009      revised.

[N,D] = size(X);          % N data points with dimension D.

if ~exist('L','var') || isempty(L)
    L = min(N,D); 
end; 
% When latent dimension is not given, set default value for it.

xc = mean(X,1);           % Obtain the centroid of original dataset.
X = X-repmat(xc,N,1);     % Zero-mean data set.               

if D<=N
    S = (X'*X);           % DxD Covariance matrix.
else
    S = (X*X');           % NxN Gram matrix.
end

[U,lambda] = eig(S);      % Eigenvalues decomposition
lambda = diag(lambda);    % Extract the eigenvalues.
[lambda,I]=sort(lambda,'descend');  % Sort the eigenvalues in descending order.
U = U(:,I(1:L));          % Extract corresponding eigenvectors.

if D>N
    lv = lambda(1:L); z = find(lv>0); lv(z) = 1./sqrt(lv(z));
    U = (X'*U)*diag(sparse(lv));
end

lambda = lambda/N;

