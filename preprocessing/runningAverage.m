function [runAvgX,runAvgY]=runningAverage(x,y,intervalWidth,stepSize,fillVal)
% ** function [runAvgX,runAvgY]=runningAverage(x,y,intervalWidth,stepSize,fillVal)
% computes the 'running average' of data points y with abscissa values x
% using a window intervalWidth wide (same unit as x) and resolution of
% stepSize (same unit as x). x must be a column vector (m x 1), y must be a
% vector of the same size or an array of size m x n. Gaps will be filled
% with fillVal (any number including NaN, which is the default, or 'last').
% Output arguments runAvgX and runAvgY are the resulting abscissa and
% ordinate values, respectively.
% The function performs operations similar to those of conv (with a
% rectangular window) but on unevenly spaced points (including multiple
% values at single abscissa values).

if nargin<=4
  fillVal=nan;
end

% sanity check
if stepSize>intervalWidth
  error('check stepSize and intervalWidth')
end

[n1, n2]=size(y);
if n1~=size(x,1)
  error('check sizes of x and y')
end

% intervals
avgIntv=(min(x)-intervalWidth/2:stepSize:max(x)-intervalWidth+intervalWidth/2+stepSize)';
avgIntv(:,2)=avgIntv(:,1)+intervalWidth;
nAvIntv=size(avgIntv,1);
runAvgX=mean(avgIntv,2);
% preallocate
runAvgY=nan(nAvIntv,n2);
for k=1:nAvIntv
  ix=x>=avgIntv(k,1) & x<avgIntv(k,2);
  if any(ix)
    runAvgY(k,:)=mean(y(ix,:));
  else
    if ischar(fillVal) && strcmp(fillVal,'last') && k>1
      runAvgY(k,:)=runAvgY(k-1,:);
    else
      runAvgY(k,:)=fillVal;
    end
  end
end

% below is an implementation using accumarray which is, unfortunately and
% strangely, much slower

% binEdges=(min(x):stepSize:max(x)+stepSize)';
% runAvgX=binEdges+stepSize/2;
% binCounts=discretize(x,binEdges);
% % note usage of new (R2016b) arithmetic
% multiBinCounts=binCounts+[1:stepIntvRatio]-ceil(stepIntvRatio/2);
% multiBinCounts=multiBinCounts(:);
% multiBinCounts=max(multiBinCounts,1);
% multiBinCounts=min(multiBinCounts,numel(binEdges));
% 
% runAvgY=accumarray(multiBinCounts,repmat(y,stepIntvRatio,1),size(binEdges),@mean,nan);
