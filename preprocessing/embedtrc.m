function [idx1,idx2]=embedtrc(nPts1,sync1,nPts2,sync2)
% ** function [idx1,idx2]=embedtrc(nPts1,sync1,nPts2,sync2)
%   Calculates indices into traces (possibly of differing length) which
%   shall be aligned in relation to a common synchronization point. trace2 
%   is embedded into trace1, which means that points in trace2 extending 
%   beyond the borders of trace1 will be sacrificed after alignment
%                        - UNITS ARE POINTS -
%
%                         >>> INPUT VARIABLES >>>
%
% NAME              TYPE/DEFAULT          DESCRIPTION
% nPts1/nPts2       scalar                number of points in trace1/trace2
% sync1/sync2       scalar OR col arr     synchronization (=trigger, if you will) points
% 
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% idx1/idx2        2-element arrays      start and stop indices into trace1 & trace2 after alignment
%
 

% ----- default values & varargin -----
verbose=0;

%pvpmod(varargin);
if verbose, disp(['**** ' mfilename ':']); end;


sz1=size(sync1);
sz2=size(sync2);
if sz1(2)>sz1(1), error('sync1 must be a column vector'); end;
if sz2(2)>sz2(1), error('sync2 must be a column vector'); end;

n=unique([sz1(1) sz2(1)]);
if size(n)~=[1 1],
  error('sync1,sync2 must have same size');
end


idx1=[ones(n,1) repmat(nPts1,n,1)];
idx2=[ones(n,1) repmat(nPts2,n,1)];

% shift of trace2 relative to trace1 required for synchronization (positive values: to the right)
shift=sync1-sync2;
% these are the indices into trace1 into which trace2 would fit
idx1=idx2+repmat(shift,1,2);

% if none of the cases below apply, this is the situation (* is the synchronization sample):
% trace1:   ----------------*-----------
% trace2:      -------------*-------             
% trace2 fits completely into trace1, no manicure required, idx1 is changed, idx2 is unchanged

% right border
cutidx=find(idx1(:,2)>nPts1);
if ~isempty(cutidx)
% trace1:   ~-----*---------
% trace2:   ~-----*--------.-------             (only points up to & including . will be taken from trace2)
  dif=nPts1-idx1(:,2);
  idx1(:,2)=idx1(:,2)+dif;
  idx2(:,2)=idx2(:,2)+dif;
end
% left border
cutidx=find(idx1(:,1)<1);
if ~isempty(cutidx)
% trace1:       -----------*---
% trace2:   ----.----------*-------             (only points starting at & including . will be taken)
  dif=idx1(:,1)-1;
  idx1(:,1)=1;
  idx2(:,1)=idx2(:,1)-dif;
end
      