function [fd,diffFIR,scaleFac]=differfi(d,si,fo,passbf,stopbf,varargin)
% ** function [fd,diffFIR,scaleFac]=differfi(d,si,fo,passbf,stopbf,varargin)
% designs a differentiator FIR filter given the input parameters, passes
% time series d through it, corrects for the lag imposed by the filter and
% removes transients at the edges of the filtered time series, setting them
% to zero.
%                         >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT          DESCRIPTION
% d              numeric array         time series
% si             scalar                sampling interval (µs)
% fo             scalar                filter order (the higher, the
%                                       steeper the filter characteristics)
% passbf         scalar                lower end of passband (Hz)
% stopbf         scalar                lower end of stopband (Hz)
% scaleFac       char or scalar, 1     if 'auto', output fd will be scaled
%                                       such that its amplitudes roughly
%                                       match those of input d (see scaleFac)
%                                      if a valid double, output fd will be
%                                       multiplied by this value
%
%                         <<< OUTPUT VARIABLES <<<
%
% NAME           TYPE/DEFAULT           DESCRIPTION
% fd             numeric array          filtered time series
% diffFIR        digital filter object  filter parameters
% scaleFac       scalar                 scaling factor


% ----- default values & varargin -----
scaleFac=1;
pvpmod(varargin);

if ischar(scaleFac) && strcmpi(scaleFac,'auto')
  doAutoScale=true;
else
  doAutoScale=false;
end
% design differentiator filter 
diffFIR=designfilt(...
  'differentiatorfir',...
  'FilterOrder',fo, ...
  'PassbandFrequency',passbf,...
  'StopbandFrequency',stopbf, ...
  'SampleRate',1e6/si);


delay=floor(mean(grpdelay(diffFIR)));
fd=filter(diffFIR,d); %/(si/1000);
% circshift to get rid of delay
fd=circshift(fd,-delay);

% substitute samples to get rid of transient
fd([1:ceil(fo/2) end-floor(fo/2):end],:)=0;

if doAutoScale
  % match 1-99 percentile range in d to .1-99.9 range in fd
  dRng=diff(prctile(d,[1 99]));
  fdRng=diff(prctile(fd,[.1 99.9]));
  scaleFac=dRng/fdRng;
end
fd=fd*scaleFac;