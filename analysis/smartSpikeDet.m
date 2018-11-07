function [tsl,thresh]=smartSpikeDet(d,si,varargin)
% ** function tsl=smartSpikeDet(d,si,varargin)
% ...'smart' implies that i) the threshold can be computed automatically,
% based on the distribution of the data, and ii) some of the computations
% may be run on the GPU, potentially speeding things up.
% It is assumed that data d are properly filtered (highpass, possibly
% notch filters applied).
% All optional input parameters must be specified as parameter/value pairs,
% e.g. as in
%          smartSpikeDet(d,si,'thresh'4);
%
%                         >>> INPUT VARIABLES >>>
% NAME        TYPE/DEFAULT          DESCRIPTION
% d           1d array              time series
% si          double                sampling interval (µs)
% thresh      double, NaN           absolute threshold; if specified, will
%                                   override normThresh (below)
% normThresh  double, 4.0           threshold, expressed as distance from 
%                                   base line in multiples of 'robustly'
%                                   computed standard deviations (sign is
%                                   irrelevant and will be computed based
%                                   on skew)
% detMode     char, 'cross'         optional input arg to tcd.m, see there
% periSpxWin  double, [-10 10]
% doUseGPU    logical, false        if true, computations will be partly 
%                                   run on GPU (likely faster than on CPU)
%
%                         <<< OUTPUT VARIABLES <<<
% NAME        TYPE/DEFAULT         DESCRIPTION
% tsl         1d array             time stamp list
% thresh      scalar               computed threshold


% ----- default values & varargin -----
thresh=NaN;
normThresh=4.0;
detMode='cross';
periSpxWin=[-10 10];
doUseGPU=false;
pvpmod(varargin,{'thresh','normThresh','periSpxWin','detMode','doUseGPU'});

if isfinite(thresh) 
  doDetThresh=false;
  if isfinite(normThresh)
    warning('absolute threshold was specified - cancelling threshold estimation')
  end
else
  doDetThresh=true;
end

if doDetThresh
  periSpxWinPts=round(periSpxWin/(si/1000));
  convKernel=ones(diff(periSpxWinPts),1);
end

tic
% convert, if requested
if doUseGPU
  d=gpuArray(d);
  if doDetThresh
    convKernel=gpuArray(convKernel);
  end
end

if doDetThresh
  % sign of threshold: if there are any spikes in the data, they should
  % skew the distribution
  s=sign(skewness(d));
  % compute median & iqr/2
  [m,v]=detbaseline(d,'meth','median');
  % threshold for removal of spikes/artifacts
  thresh=s*m+v*normThresh;
  % detect points beyond range demarcated by [-threshold threshold] 
  beyondIx=abs(d)>=abs(thresh);
  % convolute with boxcar kernel of appropriate width to extend intervals
  % of detected points by user-defined amounts
  convD=conv(beyondIx,convKernel,'same');
  beyondIx=convD>0;
  % Setting a hard threshold-based cut-off distorts the estimate of the
  % interquantile range (iqr) of the remaining data. Assuming that they are
  % normally distributed, this can be corrected for by reimplanting into d
  % data points above |thresh| and below -|thresh| that would exist had
  % they not just been removed (along with spikes). However, with
  % reasonable values of normThresh (3 and above) the proportion of these
  % data points is so tiny that this correction can be neglected.
  % Nonetheless, here is the number of points above |thresh| that would
  % have to be reimplanted:
  % halfNumPts=round((1-diff(normcdf([-1 1]*normThresh)))*numel(d)/2)
  
  % warn if the proportion of points above threshold is suspiciously large
  if sum(beyondIx)/numel(d)>.5
    warning('threshold is estimated from less than half of data points - consider decreasing periSpxWin and/or increasing normThresh')
  end
  
  % compute median & iqr/2 of points below thresh
  [m,v]=detbaseline(d(~beyondIx),'meth','median');
  % from these, compute final threshold
  thresh=s*m+v*normThresh;
end
% detect spikes
tsl=tcd(d,si,thresh,'detMode',detMode);

tsl=tsl{1};
if doUseGPU
  tsl=gather(tsl);
  if doDetThresh
    thresh=gather(thresh);
  end
end

toc