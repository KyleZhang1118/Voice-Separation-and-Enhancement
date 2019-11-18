function [E,D] = PCA(vectors, firstEig, lastEig, rankTolerance)
%The step of PCA
%vectors:Num*samples


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Default values:
if nargin < 2, firstEig = 1; end
if(nargin >2 && lastEig>size(vectors,1))
    lastEig = size(vectors,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate PCA
% module = abs(max(max(vectors))); Normalization
module=1;
vectors = vectors./module;
%Calculate the covariance matrix.
% covarianceMatrix = cov(vectors.', 1);
covarianceMatrix = vectors*vectors'/size(vectors,2);
% Calculate the eigenvalues and eigenvectors of covariance matrix.
[E, D] = eig (covarianceMatrix);

%We reserve the eigenvalues which is larger than ranktolerance,to wipe out
%the noise component.
if nargin < 4
    rankTolerance = 1e-7;
end
if nargin < 3
    lastEig_ada = sum (diag (D) > rankTolerance);
    lastEig= lastEig_ada; 
end

% Sort the eigenvalues - decending.
[~,order] = sort(diag(D),'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See if the user has reduced the dimension enough

% if lastEig < maxLastEig
%     fprintf('检测到 %d个信号源 \n',...
%            lastEig);
%     num=lastEig;
% else
%      fprintf ('信号源与传感器数量一致.\n');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drop the smaller eigenvalues
if lastEig==0 
    sprintf('no component!');
    E = [];
    D = [];
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select the colums which correspond to the desired range
% of eigenvalues.
E = E(:,order);
D = D(order,order);
E = E(:,firstEig:lastEig);
D = D(firstEig:lastEig,firstEig:lastEig)*module^2;
return;