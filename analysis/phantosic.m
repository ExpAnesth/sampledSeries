function [base,dev,phas,varargout]=phantosic(d,nPts,bin,varargin)
% ** function [base,dev,phas,varargout]=phantosic(d,nPts,bin,varargin)
%  computes the base line, variability of the base line (noise) and the
%  integral of phasic deviations from the base line from a data trace.
%  Classically, the data trace is the continuously recorded current in a
%  voltage-clamped neuron; in this case the computed parameters correspond
%  to the tonic holding current, the noise of the tonic current, and the
%  integral of the phasic postsynaptic currents. The data trace will be
%  subdivided into smaller segments of length defined by the user; all
%  parameters will be computed individually from these. The principal idea
%  of the algorithm is from Glykys and Mody, J Physiol 582.3 (2007), pp
%  1163–1178: first, the base line and the base line noise from those parts
%  of the data not containing PSCs is assessed. Then, the phasic
%  components' integral is computed. The first task can be accomplished by
%  fitting Gaussians, in which case the base line is the peak of the fitted
%  Gaussian and the noise its standard deviation. This is the (computation-
%  intensive) default. In case the data points in the base line are not
%  normally distributed, or if a faster method is required, input parameter
%  'method' should be set to 'peak'. In this case, the base line and its
%  noise will be determined directly from the amplitude distribution of the
%  data points.
%
% All input parameters except d, nPts and bin are optional and must be
% specified as parameter/value pairs, e.g. as in
%          phantosic(d,2000,-500:2:100,'method','Gauss');
%
%                         >>> INPUT VARIABLES >>>
% NAME         TYPE/DEFAULT          DESCRIPTION
% d            column array          data (sampled series)
% nPts         scalar integer        length of elementary segment in
%                                    data points
% bin          scalar                bin width of amplitude histogram
%               OR 1D array          bins of histogram
% method       char arr, 'peak'      method for base line assessment:
%                                    - 'peak' (robust)
%                                    - 'Gauss' 
% polarity     char arr, 'neg'       direction into which phasic
%                                    deflections go ('pos' or 'neg')
% prc          scalar, 15.866        method 'peak': percentile to compute 
%                                    from base line as measure of noise
%                                    (15.866 results in quasi-standard 
%                                    deviation)
% frame        scalar, 0             if a positive integer each frame-th
%                                    data excerpt will be plotted along
%                                    with histogram (and fit)
% pau          scalar, 0.1           pause (s) between each frame
% si           scalar, NaN           sampling interval in µs (only required
%                                     for plots)
%
%                     <<< OUTPUT VARIABLES <<<
% NAME         TYPE/DEFAULT           DESCRIPTION
% base         column array           base line, one datum per interval
% dev          column array           variation of data around base line,
%                                     one datum per interval
% phas         column array           average phasic current, one datum per
%                                     interval (same unit as dev and phas)
% gof          column array, optional method 'Gauss': goodness of fit 
%                                     (adjusted R2)

% to do: 
% - maybe incorporate externally determined base line

% ----- default values & varargin -----
method='peak';
polarity='neg';
prc=normcdf(-1)*100; % 15.8655253931457 
frame=0;
pau=.1;
si=NaN;
pvpmod(varargin);

% -----------------------------------------------------------------------------
%                          I. PRELIMINARIES
% -----------------------------------------------------------------------------
[n1, n2]=size(d);
if n1<nPts
  error('interval length in points must be smaller than total number of data points');
end

if ~ismember(method,{'peak','Gauss'})
  error('input var ''method'' must be ''peak'' or ''Gauss''');
end

maxNBin=400;

% --- some basics of data
switch polarity
  case 'neg'
    d=d*-1;
    bin=bin(:);
    bin=flipud(bin*-1);
  case 'pos'
    % do nothing
  otherwise
    error('input var ''polarity'' must be ''pos'' or ''neg''');
end

if frame>=1
  figH=figure;
  tmpScrSz=get(0,'Screensize');
  tmpScrSz=round([tmpScrSz([3 4]).*[.05 .4]  tmpScrSz([3 4]).*[.4 .4]]);
  set(gcf,'position',tmpScrSz);
  if isfinite(si)
    % time axis in ms
    t=(0:nPts-1)*(si/1000);
    tLabel='time (ms)';
  else
    t=0:nPts-1;
    tLabel='time (points)';
  end    
end

% --- take care of bins
if numel(bin)==1
  % ** input arg 'bin' is the bin width
  % rename
  binW=abs(bin);
  % bin borders: it is not important to capture all points of the
  % positively skewed part of the histogram, so we could theoretically cut
  % off quite a bit more than on the unskewed part (in the latter:
  % maximally 0.1% per segment). However, if there is a trend in the data,
  % this will not work, so we have to cover pretty much the whole range.
  bibo=prctile(d,[.1/(n1/nPts) 99]);
  bin=(bibo(1):binW:bibo(2))';
  % number of bins
  nBin=numel(bin);
else
  % ** input arg 'bin' contains bins
  % make sure bin is a col array
  bin=bin(:);
  nBin=numel(bin);
  binW=diff(bin(1:2));
end

if nBin>maxNBin
  warning(['number of bins is pretty large (' int2str(nBin) ') - consider purging outliers, increasing bin width or defining limited range of bins']);
end

% --- deal with time intervals
% cut off any points at end
nIntv=floor(n1/nPts);
tmp=nIntv*nPts;
disp(['given length of data and choice of interval the last ' int2str(n1-tmp) ' data points will be ignored']); 
d=d(1:tmp,:);
% reshape d
d=reshape(d,nPts,nIntv);

% -----------------------------------------------------------------------------
%                          II. COMPUTATIONS
% -----------------------------------------------------------------------------
% preallocations
base=nan(1,nIntv);
dev=base;
phas=base;
gof=base;
% amplitude histogram
n=histc(d,bin);
% peaks of individual histograms are needed for both methods
r=evdeal(n,'idx',{'minmaxpeak','allpeaks'});
% if peaks of histograms could not be determined sth is really foul
if any(~isfinite(r.maxPeakT))
  error('please check raw data and/or bins - most likely the base line shifts strongly or the bins are completely off');
end

switch method
  
  case 'peak'
    % .maxPeakT points to the lower bin border, but as the true base line
    % let's pick the middle of the interval
    base=(bin(r.maxPeakT))'+binW/2;
    % nth percentile of data points in unskewed part of histogram: pick all
    % data points smaller than lower border of (original) peak bin + half
    % of all data pts in it; then compute the 2*prc-th percentile
    for g=1:nIntv
      % take all data pts smaller than the lower border of the peak bin...
      ix1=find(d(:,g)<bin(r.maxPeakT(g)));
      % ...and the ones in the peak bin..
      ix2=find(d(:,g)>=bin(r.maxPeakT(g)) & d(:,g)<bin(r.maxPeakT(g)+1));
      % ..but retain only half of the latter
      ix2=ix2(1:2:end);
      % d(union(ix1,ix2),g) are all points in the current data segment
      % below the base line plus half of those on the base line (as defined
      % by the bins and their width). These points are putatively devoid of
      % PSCs; computing 2*prc from these is the same as computing prc from
      % the assumed underlying full (unskewed) distribution
      dev(g)=base(g)-prctile(d(union(ix1,ix2),g),2*prc);
      if ~mod(g,frame)
        plotDetail(figH, method, t, tLabel, g, d, base, bin, [], n, [])
        pause(pau);
      end
      progressbar(g/nIntv);
    end
    
  case 'Gauss'
    % start by guessing initial fit parameters for each interval
    % individually 
    tmpR=evdeal(n,'idx','minmaxpeak');
    mn=bin(tmpR.maxPeakT);
    ampl=tmpR.maxPeak;
    % roughly determine sigma from histograms
    tmp=cumsum(n);
    tmp=tmp./tmp(end);
    % set sigma to the distance between the first percentile and the mode
    [~,tmpIx]=min(abs(tmp-.01));
    sigma=mn-bin(tmpIx);
    % average
    sigmaAv=nanmedian(sigma);
    % set up fit
    ft = fittype(['a/s*exp(-0.5*((x-m)/s).^2)'] ,...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a','m','s'});
    fo = fitoptions('method','NonlinearLeastSquares');

    for g=1:nIntv
      % set identical sigma but individual amplitude and mean (baseline),
      % which is essential for long recordings with drifts
      set(fo,'Startpoint',[ampl(g)*sigmaAv mn(g) sigmaAv]);      
      % pick the first peak which attains at least 33 % of the maximal peak
      ix=find(r.posPeak{g}/r.maxPeak(g)>.33,1);
      % take all histogram bins up to that bin
      binIx=(1:r.posPeakT{g}(ix))';
      % ** note that the original Glykys & Mody paper method used the max
      % peak and more bins to the right of it, but this has been found to
      % be too conservative (the base line was shifted too much in the
      % direction of the PSC peaks)
      % binIx=(1:min(nBin,r.maxPeakT(g)+2))';
      try
        [ff,curGof]=fit(bin(binIx),n(binIx,g),ft,fo);
      catch
        warning('Bibo Binbin had a problem')
      end
      tmpPar=coeffvalues(ff);
      base(g)=tmpPar(2);
      dev(g)=tmpPar(3);
      if nargout>3
        gof(g)=curGof.adjrsquare;
      end
      if ~mod(g,frame)
        plotDetail(figH, method, t, tLabel, g, d, base, bin, binIx, n, ff)
        pause(pau);
      end
      progressbar(g/nIntv);
    end
    % same story here: as the true base line
    % let's pick the middle of the interval
    base=base+binW/2;
end

if frame>=1
  close(figH);
end

% sanity check of base line values
alpha=0.005;
if numel(base)>=5
  ix=find_outliers_Thompson(base,alpha,'biweight');
  % don't heed results of outlier correction if a third or more of the
  % intervals are labeled outliers
  if any(ix) && numel(ix)/nIntv<=1/3
    warning([int2str(numel(ix)) ' interval(s) with outlying base line values - all computed parameters in this interval will be set to NaN']);
    base(ix)=nan;
    dev(ix)=nan;
  end
end

% phasic current: mean of base line-corrected signal
phas=mean(d)-base;

% aftermath
if strcmpi(polarity,'neg')
  base=base*-1;
  phas=phas*-1;
end

base=base';
dev=dev';
phas=phas';
if nargout>3
  varargout{1}=gof';
end

function plotDetail(figH, method, t, tLabel, g, d, base, bin, binIx, n, ff)
figure(figH);
% raw data
subplot(2,1,1)
ah=area(t,d(:,g),base(g));
set(ah,'facecolor','c');
xlabel(tLabel)
ylabel('amplitude (orig. units)')
% histogram
subplot(2,1,2)
bar(bin,n(:,g),1.0,'k');
hold on
% in case of Gaussian fit, plot it
if strcmp(method, 'Gauss')
  bar(bin(binIx),n(binIx,g),1.0,'r');
  plot(bin,ff(bin),'b-');
elseif strcmp(method, 'peak')
  % base line in red
  ph=plot(base(g),0,'rd');
  set(ph,'markersize',8,'markerfacecolor','r');
  niceyuax;
end
xlabel('amplitude (orig. units)')
ylabel('N')
hold off
drawnow

