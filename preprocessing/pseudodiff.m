function [dd,pdWin]=pseudodiff(d,varargin)
% ** function [dd,pdWin]=pseudodiff(d,varargin)
% performs 'pseudo-differentiation' of time series d with the aim of
% producing a new time series in which events with a fast rise time and a
% slower decay (e.g. postsynaptic potentials) can be easily detected even
% with fluctuating or sloping baselines (see Banks at al, Neuron
% 25:449-257, 2000). Figuratively speaking, we run, sampling point by
% sampling point, a pair of short identical windows (input arg win) across
% the data, average the data points in each of them and subtract the left
% window's average ('baseline') from the right window's ('peak') one. If
% input arg deltaT is set equal to the number of points in win, the windows
% are adjacent (gap-free). The procedure is implemented as a convolution;
% the average delay of the resulting trace dd is approximately taken care
% of. For acceptable detection of e.g. cortical GABAergic IPSCs at ~35°C
% the window could be a boxcar window covering on the order of 1 ms.
% However, note that pseudodifferentiation with a boxcar window is
% essentially a convolution with a step function (that is, a funny bandpass
% filtering procedure); other shapes may be more appropriate for more
% sensitive detection of events. The user can specify an alternative window
% shape, but for more sophisticated shapes convolution should be used
% directly instead of this function.
%       ** all time units in (sampling) points **
%
% All input parameters except d are optional and must be specified as
% parameter/value pairs, e.g. as in 
%          d=pseudodiff(d,'win',ones(7,1)/7);
%
%                         >>> INPUT VARIABLES >>>
%
% NAME          TYPE/DEFAULT          DESCRIPTION
% d             array (up to 3D)      time series 
% win           vector, [1 1 1]'/3    window used for convolution
% deltaT        scalar, 3             distance between windows 
%
%                         <<< OUTPUT VARIABLES <<<
%
% NAME          TYPE/DEFAULT           DESCRIPTION
% dd            array                  pseudo-differentiated time series
%                                      with close to zero phase shift and
%                                      equal length as d
% pdWin         vector                 filter coefficients resulting from
%                                      pseudodifferentiation

% (C) Harald Hentschke, University of Tübingen
% Last code change Oct 14 2016

% -----------------------------------------------------------------------------
%                          I. PRELIMINARIES
% -----------------------------------------------------------------------------
% ----- default values & varargin -----
% convolution window
win=[1 1 1]'/3;
% delta t for computation of difference (=time difference between baseline
% and peak windows, measured between identical points) in pts
deltaT=3;
si=[];

pvpmod(varargin);
disp(['** ' mfilename]);

% --- checks of input
[n1 n2 n3]=size(d);
if n1<=1
  error('input array ''d'' does not contain a sufficient number of rows');
end

[tmp1 tmp2]=size(win);
if tmp1*tmp2==1
  error('''win'' is a scalar but must be an 1D array of at least 2 elements');
elseif all([tmp1 tmp2]>1)
  error('''win'' is a 2D array but must be a 1D array');
end
% now make sure it's a column array
win=win(:);
nWinPts=numel(win);

if deltaT<1
  error('''deltaT'' must be positive and equal to or greater than nWinPts');
elseif deltaT<nWinPts
  error('a value of deltaT smaller than nWinPts does not make sense - the desired effect can be achieved more efficiently by shortening win');
end

% -----------------------------------------------------------------------------
%                          II. COMPUTATIONS
% -----------------------------------------------------------------------------
% construct the full window for convolution 
pdWin=[win; zeros(deltaT-nWinPts,1); -win];
nPdWinPts=numel(pdWin);
% compute group delay (§ could interpolate trace to account for group delay
% more precisely)
delay=floor(mean(grpdelay(pdWin)));

% preallocate
dd=zeros([n1 n2 n3]);
for s=1:n3
  for c=1:n2
    % convolute; compared to d tmpd is longer by nPdWinPts-1 points
    tmpD=conv(d(:,c,s),pdWin);
    % so, to get rid of delay 
    dd(:,c,s)=tmpD(1+delay:end-(nPdWinPts-1-delay));
  end
end


% % ------ historical, cumbersome code below (changed Oct 14 2016) -------
% % compute length difference between d(:,c,s) and tmp after computation of
% % difference of convoluted traces
% lDiff=nWinPts-1-deltaT;
% % preallocate
% dd=zeros([n1 n2 n3]);
% for s=1:n3
%   for c=1:n2
%     % convolute
%     tmp=conv(d(:,c,s),win);
%     % compute difference
%     tmp=tmp(1+deltaT:end)-tmp(1:end-deltaT);
%     % embed: the center of vector tmp should be aligned to center of
%     % vector d(:,c,s), but tmp differs in length by nWinPts-1-deltaT (as
%     % computed above). So,
%     if lDiff<=0
%       % tmp is shorter than or as long as current slot in dd
%       dd((ceil(abs(lDiff)/2):end-floor(abs(lDiff)/2)-1)+rem(abs(lDiff)+1,2),c,s)=tmp;
%     else
%       % tmp is longer than current slot in dd
%       dd(:,c,s)=tmp(1+floor(lDiff/2):end-ceil(lDiff/2));
%     end
%   end
% end
% 



% ------ historical, problematic code (in terms of alignment) below -------
% (code was changed 8. Oct 2016)
%
% winHalfPts=(nWinPts-1)/2;
% ...
% % cut off additional points at either end of array resulting from
% % convolution
% tmp=tmp(1+winHalfPts:end-winHalfPts);
% % compute difference
% tmp=tmp(1+deltaT:end)-tmp(1:end-deltaT);
% % correct for the rightward shift of the first term above by deltaT
% % points
% dd(1+deltaT:end,c,s)=tmp;
