function [xe,fz,c,poolSD]=xapen(d,m,varargin)
% ** function [xe,fz,c,poolSD]=xapen(d,m,varargin) computes
% cross-approximate entropy (XApEn) of all non-redundant pairs of data
% columns ('channels') in input variable d and places the values in output
% variable xe. Channel pairs in xe follow the order 1-1, 1-2, ..., 2-2,
% 2-3,..., which is listed in output variable c. By default, 'auto'
% approximate entropy will be computed, but this option can be turned off.
% All input parameters except d and m are optional and must be specified as
% parameter/value pairs, e.g. as in
%          xe=xapen(d,2,'tol',0.4);
%
%               >>> INPUT VARIABLES >>>
% NAME        TYPE/DEFAULT          DESCRIPTION
% d           array                 data (channels in columns)
% m           scalar, 1             subsequence length
% tol         scalar or array, 0.2  tolerance: if tolType is 'sd', tol must
%                                    be a scalar; otherwise it may be a
%                                    scalar or an array with the elements
%                                    corresponding to the channel pairs
% tolType     char array, 'sd'      'sd': tolerance specified in terms of
%                                    standard deviations (pooled for each
%                                    channel pair separately)
%                                   'abs': tolerance specified in absolute
%                                     units
% doApen      logical, true         if true, will include 'auto' approxima-
%                                    to entropy (XApen of a time series
%                                    with itself) for all data columns
% doNorm      logical arr, [0 0]    if first element is true, each data
%                                    column will be mean-subtracted
%                                   if second element is true, each data
%                                    column will be divided by its sd (in
%                                    which case tolType 'abs' is not
%                                    allowed)
%                                    
%             <<< OUTPUT VARIABLES <<<
% NAME        TYPE/DEFAULT           DESCRIPTION
% xe          column array           XApEn, as many elements as there are
%                                     non-redundant channel pairs
% fz          column array           the fraction of segments of length m+1
%                                     for which no single match was found
%                                     (one value for each channel pair)
% c           N by 2 array           list of channel pairs (which can also 
%                                     be obtained outside this function via
%                                     combin(size(d,2))
% poolSD      column array           pooled standard deviations (one for
%                                     each channel pair)
%
% *** PLEASE NOTE: this implementation is reasonably fast because for each
% pair of channels it computes all possible pairwise differences between
% elements in channel1 and channel2 only once. The drawback of this
% approach is that the intermediate matrix holding these values requires a
% lot of memory (e.g. ~750 MB for data with 10000 points per channel). If
% your data are much larger than this you should consider using a different
% algorithm. Alternatively, consider converting input data d to single
% number type.
%
% -------------------------------------------------------------------------
% Version 1.0, Oct. 2012. Code by Harald Hentschke (University of Tübingen)
% -------------------------------------------------------------------------



% to do/improvements
% - more elaborate checks of tol
% - compute memory demand before running computations and exit if
% insufficient memory available
% - parallel computing!

% defaults
tol=.2;
tolType='sd';
doApen=true;
doNorm=[false false];
pvpmod(varargin);

% -------------------------------------------------------------------------
% ------- PART I: CHECKS OF INPUT ARGS & PREPARATORY COMPUTATIONS ---------
% -------------------------------------------------------------------------
[n1 n2]=size(d);
if n1<n2
  error('data channels must be along columns');
end
if n2<2 && ~doApen
  error('d must contain at least two data columns for strict XApEn analysis');
end
if ~any(strcmp(tolType,{'abs','sd'}))
  error('illegal value for ''tolType''');
end
if strcmp(tolType,'abs') && doNorm(2)
  error('specifying absolute tolerance with data to be normalized does not make sense');
end

% pairs of channels
[c,nc]=combin(n2,'autoC',doApen);

% subtract mean?
if doNorm(1)
  d=d-repmat(mean(d),n1,1);
end

% variance of channels
dVar=var(d);
% pooled sd = square root of [mean of variances]
poolSD=sqrt(sum(dVar(c),2)/2);

% tol is expanded here, one entry for each nonredundant channel pair
switch tolType
  case 'abs'
    if numel(tol)==1
      tol=repmat(tol,1,nc);
    elseif numel(tol)==nc
      tol=tol(:)';
    else
      error('tol must be a scalar or an array with as many elements as channel pairs');
    end
  case 'sd'
    if numel(tol)~=1
      error('tol must be a scalar if tolType is set to ''sd''');
    end
    if doNorm(2)
      % normalize channels by division by their sd
      d=d./repmat(sqrt(dVar),n1,1);
      % expand tol
      tol=repmat(tol,1,nc);
    else
      % tolerance, scaled to pooled sd
      tol=tol*poolSD';
    end
end

% number of subsegments of length m
nSeg=n1-m+1;
% preallocate xe
xe=nan(1,nc);

% -------------------------------------------------------------------------
% ------- PART II: XAPEN COMPUTATION ---------
% -------------------------------------------------------------------------

% In the implementation here, first all non-redundant pairwise differences
% between elements in channel1 and channel2 are computed and placed in a
% square matrix. This matrix is converted into logical matrix cm,
% indicating whether the given difference is within tolerance r (true) or
% not (false). In this matrix, similar sequences of length m in the
% original data are represented as diagonal sequences of logical ones of
% length m. So, these have to be identified and counted, which is achieved
% by iterative logical ANDing subsets of cm, and then just summing up the
% remaining column entries.

% loop over first channel within the pairs to be compared
for chIx1=unique(c(:,1))'
  ch1=d(:,chIx1);
  % loop over second channel
  for chIx2=c(c(:,1)==chIx1,2)'
    ch2=d(:,chIx2);
    % index into c to current channel pair
    cIx=find(chIx1==c(:,1) & chIx2==c(:,2));
    % tolerance for current channel pair
    curTol=tol(cIx);
    % now construct matrix of logicals indicating whether differences
    % between all non-redundant pairwise combinations of points in channel1
    % and channel2 are smaller than tol
    cm=abs(bsxfun(@minus,ch1',ch2))<curTol;
    % counting blocks of m occurrences of logical(1) along diagonals:
    % iteratively ANDing 'opposing corner'-subsets of cm
    for mIx=1:m-1
      cm=cm(1:end-1,1:end-1) & cm(2:end,2:end);
    end
    % count entries in the columns of the remaining cm and divide the
    % counts by N-m+1, the number of comparisons possible between one
    % segment of length m from channel1 and all segments of length m of
    % channel2. If there is any segment in channel1 which is not similar to
    % any other segment in channel2, the column mean is zero and XApEn in
    % the strict sense is undefined as in the next step the logarithm of
    % that mean count is taken. As a workaround, we here simply ignore the
    % segment so that it does not spoil the entire XApEn value for the
    % current channel pair, and count the fraction of instances this
    % happens. Finally, compute the average of such comparisons possible
    % (which is also N-m+1, because that's the number of possible segments
    % in channel 1).
    cmMn=sum(cm)/nSeg;
    simIx=find(cmMn>0);
    ctm=sum(log(cmMn(simIx)))/numel(simIx);
    % repeat steps above for segment length m+1... 
    cmMn=sum(cm(1:end-1,1:end-1) & cm(2:end,2:end))/(nSeg-1);
    simIx=find(cmMn>0);
    ctm1=sum(log(cmMn(simIx)))/numel(simIx);
    % ...and now keep the fraction of segments without a match (because
    % this number will be higher than that for segment length m)
    fz(cIx)=1-numel(simIx)/(nSeg-1);
    % finally, XApEn
    xe(cIx)=ctm-ctm1;
  end
end

