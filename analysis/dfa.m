function [fluc,win,alpha,stats]=dfa(d,si,varargin)
% ** function [fluc,win,alpha,stats]=dfa(d,si,varargin)
% performs detrended fluctuation analysis (DFA) on a single-channel time
% series with the aim of detecting/determining power-law scaling behavior
% as an expression of 'self-organized criticality' in nonstationary
% signals. NOTE: input data d should be the envelope of a bandpass-filtered
% signal, otherwise the algorithm makes little sense.
% Brief description: 
% - divide time series into segments of equal length (window size)
% - for each segment subtract the local trend (linear fit) 
% - compute fluctuation (root mean square)
% - average values from all segments
% - repeat for all window sizes
% The next (optional) step is to plot the mean fluctuations against the
% corresponding window sizes on a double logarithmic scale and compute a
% linear fit. Within the framework of 'self-organized criticality' the
% slope of this fit (and, foremost, the question whether it is a linear
% slope at all) informs us about....
%
% Optional input parameters listed below (= all except d and si) must be
% specified as parameter/value pairs, e.g. as in 
%   [fluc,win]=dfa(d,si,'nWin',30);
%
% --------------- INPUT VARIABLES ------------------------------------
% NAME      TYPE, DEFAULT      DESCRIPTION
% d         column array       time series
% si        scalar             sampling interval in microseconds
% nWin      scalar, 40         number of windows
% win       2el-array,         minimal and maximal window length in points;
%             [10 numel(d)/2]   mininal length should be substantially
%                               longer than one period of the dominant 
%                               frequency; maximal window length should not
%                               longer than half of data length
% fh        figure handle, []  if a valid figure handle, plot will be generated
% ---------------- OUTPUT VARIABLES ----------------------------------
% fluc      row array          mean fluctuation values for each window size
% win       row array          window lengths
% alpha     scalar             slope of linear fit
% stats     struct             'stats' output of regress function

%
% For details and background of DFA see:
% 1. Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE and Goldberger
%    AL. Mosaic organization of DNA nucleotides. Physical Review E 49:
%    1685-1689, 1994.
% 2. Linkenkaer-Hansen K, Nikouline VV, Palva JM and Ilmoniemi RJ.
%    Long-range temporal correlations and scaling behavior in human brain
%    oscillations. J Neurosci 21: 1370-1377, 2001.

% check time series
[n1,n2]=size(d);
if n1>1 && n2>1
  error ('d must be a column vector')
elseif n2>n1
  d=d(:);
  [n1,n2]=size(d);
end

% defaults
nWin=40;
win=[10 floor((length(d)/2))];
fh=[];
% optional inputs 
pvpmod(varargin);

% check inputs
if ~isequal(nWin,floor(nWin))
  error('nWin must be an integer');
end
if win(1)<3
  error('min window length must be > 3');
end
if win(2)>n1
  error('max window length must be <= data length')
end

% first step is some data pre-processing
d=d-mean(d);
d=cumsum(d);

% array of logarithmically spaced window lengths
win=round(logspace(log10(win(1)),log10(win(2)),nWin));
% preallocation of fluc
fluc=zeros(size(win));

% loop through window lengths
for g=1:numel(win)
  % define start and stop times of excerpts overlapping by half 
  [intv]=mkintrvls([1 n1],'ilen',win(g),'olap',ceil(win(g)/2));
  % number of data segments
  nSeg=size(intv,1);

  % preallocation of variable holding rms value for each segment
  rmsd=zeros(nSeg,1);

  % prepare linear regression
  x=(1:win(g))';
  x=[ones(size(x)) x];

  % deal with segments
  for h=1:nSeg;
    % regression
    [rcoeff,confI,r]=regress(d(intv(h,1):intv(h,2)),x);
    % function regress provides the residuals as third output arg
    rmsd(h)= sqrt(mean(r.^2));
  end
  fluc(g)=mean(rmsd);
end

% final task is to fit a line on a loglog plot:
% convert time to s and log-transform both variables
logWin=(log10(win*(si/1e6)))';
logFluc=(log10(fluc))';

% prepare regression
logWin=[ones(size(logWin)) logWin];
% regress
[rcoeff,confI,r,rint,stats]=regress(logFluc,logWin);
alpha=rcoeff(2);
disp(['alpha=' num2str(rcoeff(2))]);

if ~isempty(fh) && ishandle(fh)
  figure(fh)
  % don't clear axes to allow multiple plots
  subplot(2,1,1), hold on
  plot(d);
  niceyax;
  % we are particularly interested in the slope of the fitted line
  disp(['linear regression y=a+b*x: a=' num2str(rcoeff(1)) '; slope of the line b=' num2str(rcoeff(2)) ]);

  % now we want to plot the regression line on a loglog plot together with
  % the data; for this, values of the line must be generated and transformed
  line4log=10.^(rcoeff(2)*logWin(:,2) + rcoeff(1));
  subplot(2,1,2), hold on
  % plot the results on double logarithmic scale
  plot(win,fluc,'o',win,line4log, 'k');
  axis tight
  set(gca,'xscale','log');
  set(gca,'yscale','log');

  xlabel('{\tau} (points)');
  ylabel('mean fluctuation');

  disp(['DFA computation finished, total data length t=' num2str(length(d)*(si*1e-6)) ' sec;'])
  disp (['number of windows=' num2str(nWin) '; window size range from ' num2str(win(1)*(si*1e-6)) ' to ' num2str(win(end)*(si*1e-6)) ' sec'])

end