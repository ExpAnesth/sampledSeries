function [b,v,varargout]=detbaseline(d,varargin)
% ** function [b,v,varargout]=detbaseline(d,varargin)
% computes value of base line of sampled series d
% All input parameters except d are optional and must be specified as
% parameter/value pairs, e.g. as in 
%          detbaseline(d,'meth','xxx');
%
%                         >>> INPUT VARIABLES >>>
%
% NAME          TYPE/DEFAULT          DESCRIPTION
% d             array                 sampled series (up to 2D, columns)
% meth          char, 'mode'          method of base line detection:
%                                       'mode' - most frequent value in d
%                                       'median'
% relBinWid     scalar, .002          bin width, relative to amplitude range,
%                                      to be used in 'mode' method
%
%                         <<< OUTPUT VARIABLES <<<
%
% NAME          TYPE/DEFAULT           DESCRIPTION
% b             array                  values of base line (mode or median)
% v             array                  variability of base line
% varargout{1}  array                  base line-corrected time series 
%
% 


% ----- default values & varargin -----
meth='mode';
relBinWid=.002;
pvpmod(varargin);

[n1,n2,n3]=size(d);

% preallocate
b=nan(1,n2,n3);
v=nan(1,n2,n3);

switch meth
  case 'mode'
    for g=1:n2*n3
      % assume that the most frequent value in d is within the 1st and 99th
      % percentile
      p=prctile(d(:,g),[1 99]);
      bin=linspace(p(1),p(2),round(1/relBinWid));
      [~,ix]=max(histcounts(d(:,g),bin));
      b(g)=mean(bin([ix ix+1]));
      % § no idea for v yet
    end
  case 'median'
      p=prctile(d,[15.8655253931457 50 84.1344746068543],1);
      b=p(2,:,:);
      % interquantile range divided by half
      v=diff(p([1 3],:,:),1,1)/2;
  otherwise
    error('illegal base line detection method')
end

if nargout>1
  % varargout{1}=d-repmat(b,[n1 1 1]);
  varargout{1}=d-b;
end
