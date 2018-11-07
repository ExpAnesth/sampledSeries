function nd=halign(d,intv,alignMeth,ipFac)
% ** function nd=halign(d,intv,alignMeth,ipFac)
% aligns data traces (e.g. single action potentials) horizontally according 
% to distinctive characteristics in trace (e.g. humps or flanks) 
%
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% d                2D array              data traces (along columns)
% intv             2 element array       interval (points) in which to look
%                                         for charactteristics (max slope,
%                                         peak, etc)
% alignMeth        char arr              'posSlope','negSlope','posPeak','negPeak',
% ipFac            scalar                interpolation factor
%                     
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT           DESCRIPTION
%


% § checks 

[n1,n2]=size(d);

% the maximal number of points any trace will be shifted given interval
% above
maxShift=diff(intv([1 end]));

x=(1:n1)';
xi=(1:1/ipFac:n1)';
nXi=numel(xi);
% indices into xi to points whithin which to determine slope
iSlopIx = find(xi>=x(intv(1)) & xi<=x(intv(end))) ;

% preallocate array into which to place the aligned and downsampled traces
nd=repmat(nan,n1+maxShift,n2);
% interpolate
interpD=interp1(x,d,xi,'spline');
% determine parameters in given interval
tmpr=evdeal(interpD(iSlopIx,:),'idx',{'minmaxslope','minmaxpeak'});

% these are the indices to the steepest points in the interpolated traces
% corrected for the offset imposed by interval above
switch alignMeth
  case 'posSlope'
    % iSlopIx(1)-1 to account for the fact that indices in matlab start
    % with 1, +1 to define the second (of the two) point defining a slope
    % as the one
    sMaxT=(tmpr.slopemaxT+iSlopIx(1)-1+1)';
  case 'negSlope'
    sMaxT=(tmpr.slopeminT+iSlopIx(1)-1+1)';
  case 'posPeak'
    sMaxT=(tmpr.maxPeakT+iSlopIx(1)-1+1)';    
  case 'negPeak'
    sMaxT=(tmpr.minPeakT+iSlopIx(1)-1+1)';    
  otherwise
    error('illegal alignment method')
end
clear tmpr;
% in case there are events with no slope or peak at all (flat lines) impose
% arbitrary shift 
sMaxT(isnan(sMaxT))=min(sMaxT);

% as the traces will be downsampled compute the index to the steepest slope
% of the downsampled traces under the assumption that the point of steepest 
% slope will be in the reading raster
sync2=floor((sMaxT-1)/ipFac)+1;

% downsample & align
for g=1:n2
  % pre-alignment point points
  iIx=fliplr(sMaxT(g):-ipFac:1);
  % all
  iIx=[iIx sMaxT(g)+ipFac:ipFac:nXi];
  % embed
  nd((1:numel(iIx))-sync2(g)+intv(end),g)=interpD(iIx,g);
end

% prune nd: maxShift at either end
nd([1:maxShift end-maxShift+1:end],:)=[];


