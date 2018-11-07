function xapen_demo(d,m,varargin)
% ** function xapen_demo(d,m,varargin) is a variant of xapen.m
% intended to illustrate graphically the computation of cross-approximate
% entropy (XApEn) between all non-redundant pairs of data columns
% ('channels') in input variable d. For further information see xapen.m.
% All input parameters except d and m are optional and must be specified as
% parameter/value pairs, e.g. as in
%          xe=xapen_demo(d,2,'tol',0.4);
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
    % variance of channels
    dVar=var(d);
    if doNorm(2)
      % normalize channely by division by their sd
      d=d./repmat(sqrt(dVar),n1,1);
      % expand tol
      tol=repmat(tol,1,nc);
    else
      % tolerance, scaled to pooled sd
      % (** note: pooled sd = square root of [mean of variances])
      tol=tol*sqrt(sum(dVar(c),2)/2)';
    end
end

% number of subsegments of length m
nSeg=n1-m+1;

% -------------------------------------------------------------------------
% ------- PART II: XAPEN COMPUTATION & GRAPHICS ---------
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

% set up figure for plotting
fh=figure; clf
colormap(coma('bluered','ncols',300));
set(gcf,'DefaultLineMarkersize',5);
scs=get(0,'screensize');
marg=round(scs(4)/40);
set(fh,'position',[scs(1)+marg floor(scs(4)/2)-marg scs(3)-2*marg floor(scs(4)/2)-2*marg]);


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
    dm=bsxfun(@minus,ch1',ch2);
    cm=abs(dm)<curTol;
    
    % plot this matrix
    ma=max(dm(:));
    mSph=subplot(1,3,3);
    hold on
    imagesc(dm);
    set(mSph,'clim',[-ma ma]);
    spy(cm,'k.');
    xlabel('channel 1');
    ylabel('channel 2');
    title('difference matrix')
    axis equal
    axis tight
    % now plot the raw data
    pos=get(mSph,'position');
    % assume margins of 0.05
    rSph=subplot('position',[.05 pos(2) 1-pos(3)-.2 pos(4)]);
    hold on
    ph=plot(d,'-');
    set(ph(1),'color',[.4 .4 .6]);
    set(ph(2),'color',[.6 .4 .4]);
    ph=plot(d,'k.');
    set(ph,'markersize',10);
    nicexyax(40);
    
    % counting blocks of m occurrences of logical(1) along diagonals:
    % iteratively ANDing 'opposing corner'-subsets of cm
    for mIx=1:m-1
      cm=cm(1:end-1,1:end-1) & cm(2:end,2:end);
    end
    % repeat, looking for m+1, but save as different var
    cm1=cm(1:end-1,1:end-1) & cm(2:end,2:end);
    
    % work through columns of cm (=subsegments of channel1)
    for gIx=1:nSeg-1
      mIx=[];
      % index to similar subsegments of length m of channel2 
      mIx{1}=find(cm(:,gIx));
      % index to similar subsegments of length m+1 of channel2 
      mIx{2}=find(cm1(:,gIx)); 
      ph1m={};
      ph2m={[]};
      subplot(rSph);
      for k=0:1
        if m+k==1
          lineSty='none';
        else
          lineSty='-';
        end
        % current subsegment of channel 1
        d1Ix=gIx:gIx+m-1+k;
        ph1m{k+1}=plot(d1Ix,d(d1Ix,1),'bo-');
        set(ph1m{k+1},'linestyle',lineSty);
        if ~isempty(mIx{k+1})
          % matching segments of channel 2
          d2Ix=ones(m+k,1)*mIx{k+1}'  +  (0:m+k-1)'*ones(1,numel(mIx{k+1}));
          ph2m{k+1}=plot(d2Ix,reshape(d(d2Ix,2),[m+k numel(mIx{k+1})]),'mo-');
          set(ph2m{k+1},'linestyle',lineSty);
        end
      end
      if numel(ph1m)>1
        set(ph1m{2},'markersize',8);
      end
      if numel(ph2m)>1
        set(ph2m{2},'markersize',8);
      end
      th=smarttext(...
        [int2str(numel(mIx{1})) ' similar segments (m=' int2str(m) '); ', ...
        int2str(numel(mIx{2})) ' similar segments (m=' int2str(m+1) ')'] ...
        );
      title('button to continue');
      pause
      for k=1:numel(ph1m)
        delete(ph1m{k});
      end
      for k=1:numel(ph2m)
        delete(ph2m{k});
      end
      delete(th);
    end
    
  end
end

