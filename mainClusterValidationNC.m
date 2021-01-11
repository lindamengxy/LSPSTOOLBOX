% Toolbox for Estimating the Number of Clusters (Version 2.0)
% based on the clustering results of PAM or K-means clustering algorithm
% Please read the help file "Readme.txt" before running this program
% Kaijun WANG: wangkjun@yahoo.com, April 2007.
function mainClusterValidationNC(X)
%clear;
alg = 2;                     % 1 --- PAM, 2 --- K-means
newp=1;                   % opening a new figure window ?
subp=0;                    % plotting in which sub-window
pcolor=1;                  % color plotting ?
pc=1;                        % plotting data by Principal Component Analysis
staz=0;                     % standardization to [0 1] when mixed metric
Hsep=0;                   % threshold for System Evolution method
nk=1;                        % skip computation of error rate if NC is unknown
N2=10;                     % searching limit is max(N2,nk+6)
type = 1;                   % 1 - using Euclidean distances for general data;
                                  % 2 - Pearson correlation coefficients for gene data
                                 % it is preseted in row 64: if id > 20  type = 2; end

data = X;
[nrow, dim] = size(data);
N1=2;                          % a forward search starting at k = 2
N=10;%floor(max([N2,nk+6])/2);     % stopping at k = N
truelabels = ones(nrow,1);


% calculating dis-similarity/distance matrix of a data set
if type == 2
  % Pearson similarity [-1,1] is normalized to Pearson distance [0,1]
  Dist = 1-(1+similarity_pearson(data'))/2; 
  for j = 1:nrow
     Dist(j,j) = 0;
  end
  dc = 2;                         % a sign to call energy_pearson.m
  dmax = 1;                    % max value of data
else
   if staz 
       data = standarz(data); end            % columns are standardized within [0,1]
   [Dist,dmax]= similarity_euclid(data); % Euclidean distances between rows
   dc = 1;                       % a sign to call energy_euclid.m
end
dissim = [];
if  alg == 1
  for i = 2:nrow
   dissim = [dissim Dist(i,1:i-1)];  % dissimilarity vector
  end
end

% Part 2: Running PAM or K-means clustering algorithm
classlabel = ones(nrow,N);
vtype = 4*ones(1,dim);
Rd='euclidean';

opts = statset('Display','final');
for i = 1:N %N1 : N
  fprintf('\n  = = running on K-means clustering at k= %i', i);
  classlabel(:,i) = kmeans(X, i, 'distance', 'sqeuclidean','Replicates',5,'Options',opts);
  Q = ind2cluster(classlabel(:,i)); 
  % less than 4 data points in one cluster, stop !
  ns = [];
  for j =1:numel(Q)
     ns(j) = numel(Q{j});
  end
  if min(ns) < 4
     N = i-1;
     break;
  end
end

% Part 3:  Estimating the number of clusters by validity indices
validity_Index
