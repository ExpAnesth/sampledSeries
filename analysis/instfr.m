function [ifr,time,krnl,kscal]=instfr(tsl,varargin)
% ** function [ifr,time]=instfr(tsl,varargin)
% generates from time stamp list (discrete points in time, e.g. spike
% train) an estimate of the instantaneous firing rate (= a quasi-continuous
% equivalent)
%  
%   -- function half-finished -->> see ievrate.m for a successor
%                        >>> INPUT VARIABLES >>>
%
% NAME         TYPE/DEFAULT          DESCRIPTION
%
% tsl          1d array              time stamp list (unit: ms)
% resol        scalar, 1000          time resolution of the ifr (unit: us)
% meth         string, 'hist'        the method to use for conversion:
%                                    'hist' - spikes are binned regularly, 
%                                      producing a simple histogram
%                                    'krnl' - convolute spike train with 
%                                      kernel function 
% muetze       char arr,'triang'     type of kernel function:
%                                    'triang' - triangular window
%                                    'gauss' - Gaussian window
% sigma        scalar, 20            standard width of kernel in ms
% Int          2-element arr         ms, if specified, ifr will cover 
%                                     exactly this interval (even if this 
%                                     amounts to putting out many or all zeros)
%
%                         <<< OUTPUT VARIABLES <<<
% NAME         TYPE/DEFAULT          DESCRIPTION
% ifr          1d (column) array     instananeous firing rate. If time 
%                                     stamp list is empty, an empty array 
%                                     will be put out, unless 'Int' had 
%                                     been specified
% time         1d (column) array     time (for 'hist' method these are the 
%                                     centers of the bins)
% krnl         2-column array        the kernel function and its time axis 
%                                      (empty if method other than krnl chosen)
% kscal        scalar                the factor by which the kernel will be divided to have unit area

% default values
resol=1000;
meth='hist';
muetze='triang';
sigma=20;
muetzenplotz=0;
Int=[];
ifr=[];
time=[];

pvpmod(varargin);
disp(['**** ' mfilename ':']);

% first thing to do: convert resol from us to ms - more convenient
resol=resol*.001;

% ----- preparing kernel func
switch muetze
  case 'gauss'
    % total half-length of muetze in ms - for Gaussian, set to 3*sigma
    mHalfLen=3*sigma;
    % same in ticks
    mHalfLen_pts=round(mHalfLen/resol);
    % by definition, krnl should contain an odd number of points and be
    % symmetrical 
    krnl=exp(-.5*(([-mHalfLen_pts:mHalfLen_pts]*resol).^2)/sigma^2);
  case 'triang'
    mHalfLen=sqrt(6)*sigma;
    % same in ticks
    mHalfLen_pts=round(mHalfLen/resol);
    krnl=(triang(2*mHalfLen_pts+1))';
  otherwise
    error('no valid muetze specified');
end;

% scale krnl such that area==1
kscal=sum(krnl);
krnl=krnl/kscal;
if muetzenplotz
  mpfh=figure;
  ph=plot([-mHalfLen_pts:mHalfLen_pts]*resol,krnl);
  set(ph,'linewidth',2);
  xlabel('ms');
  % close(mpfh);
end;



compute=1;
% check Int and tsl and decide whether further processing makes sense given
% the choice of Int
if isempty(Int)
  if isempty(tsl),
    warning('time stamp list is empty and ''Int'' is not specified - returning empty arrays');
    compute=0;
    % this is redundant
    time=[];
    ifr=[];
  else
    % further processing makes sense, make Int cover all spikes in tsl
    %§
    if length(tsl)>1, Int=[tsl(1) tsl(end)+resol];
    else Int=[tsl tsl+resol];
    end
  end
else
  % perform simple checks on Int
  if sum(size(Int))~=3, error('Int must be a 2-element array'); end;
  if diff(Int)<=0.0, error('time interval for ifr is zero or bounds are reversed in order'); end;
  % at this stage it is safe to generate time array. remember: centers of
  % bins
  % §
  time=Int(1)+resol/2:resol:Int(2);
  if isempty(tsl),
    compute=0;
    warning('time stamp list is empty - returning array of zeros');
    ifr=zeros(size(time));
  else
    % further processing makes sense unless..
    if Int(1)>tsl(end) | Int(2)<tsl(1), 
      compute=0;
      warning('time interval for ifr does not cover any of the spikes in tsl - returning array of zeros');  
      ifr=zeros(size(time));
    end
  end
end

if compute,
  switch meth
  case 'hist'
    % remember: centers of bins
    % §
    time=Int(1)+resol/2:resol:Int(2);
    % histc interprets elements of 'time' as edges (not as centers like
    % 'hist'), so shift time by resol/2 
    ifr=histc(tsl,time-resol/2);
    % kill last bin because it counts only values lying exactly on the 
    % border
    ifr(end)=[];
    time(end)=[];
    ifr=ifr/resol*1e3;
  case 'krnl'
    % cut down tsl
    tsl=tsl(find(tsl>=Int(1) & tsl<Int(2)));
    % convert time stamps in ms to time stamps in ticks
    %§
    tsl_pts=round(tsl/resol);
    % if any two spikes are assigned the same time stamp (expressed in
    % ticks) time resolution is not adequate to resolve the difference
    % between them In the algorithm below this is not tragic, but when
    % convolution is applied make sure the HISTOGRAM is convoluted
    tmp1=diff(tsl_pts);
    if any(~tmp1),
      warning(['time resolution for kernel-computed ifr is not adequate: ' num2str(100*2*length(find(~tmp1))/length(tsl)) ' % spx assigned to same bin']);
    end
%    [time,time_pts,ndistance,ndistance_pts]=regspace(Int,'resol',resol,'distance',resol);
    [time,time_pts,ndistance,ndistance_pts]=mkintrvls(Int,'resol',resol,'ilen',resol,'border','skip');
    % centers of bins, please - only needed for time
    time=makerow(time(:,1)+ndistance/2-1);
    time_pts=makerow(time_pts(:,1));
    spTrain_pts=zeros(size(time_pts));
    % find spikes too close to border
    borderTsl_pts=tsl_pts(find(tsl_pts>(time_pts(end)-mHalfLen_pts) | tsl_pts<time_pts(1)+mHalfLen_pts));
    % 
    for i=1:length(borderTsl_pts)
      idx=[-mHalfLen_pts:mHalfLen_pts]+borderTsl_pts(i);
      idx_idx=find(idx>=time_pts(1) & idx<=time_pts(end));
      spTrain_pts(idx(idx_idx)-time_pts(1)+1)=spTrain_pts(idx(idx_idx)-time_pts(1)+1)+krnl(idx_idx);
    end;
    % 
    tsl_pts=setdiff(tsl_pts,borderTsl_pts);
    for i=1:length(tsl_pts)
      idx=[-mHalfLen_pts:mHalfLen_pts]+tsl_pts(i)-time_pts(1)+1;
      spTrain_pts(idx)=spTrain_pts(idx)+krnl;
    end;
    ifr=spTrain_pts;
    
    % convolute
    % convolution is stupid for 'sparse' spike trains; 
    % adding a kernel per spike should be more effective
%    ifr=conv(spTrain_pts,krnl);
%    time=(tsl_pts(1)-mHalfLen_pts:tsl_pts(end)+mHalfLen_pts)*resol;
  end
end
ifr=ifr(:);
time=time(:);
krnl=[[-mHalfLen_pts:mHalfLen_pts]*resol; krnl]';