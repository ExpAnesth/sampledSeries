function d=elim_hum(d,si,f,varargin)
% **function d=elim_hum(d,si,f,varargin)
% eliminates typical line pickup ('hum') in data d by capping, in the
% frequency domain, the magnitude of the signal in frequency range
% f +/- sExt Hz (default: 0.2) to values in the surround (f +/- estHalfWid Hz;
% default: 0.4). si is the sampling interval of the data in µs; f is the
% array of hum center frequencies to be treated this way. 
%
%                    >>> INPUT VARIABLES >>>
% NAME            TYPE/DEFAULT           DESCRIPTION
% d               array (up to 3D)       sampled data, time runs along column
% si              scalar                 sampling interval in µs
% f               array                  frequencies to be eliminated
% substHalfWid   scalar, 0.2             substitution half width (Hz)
% estHalfWid     scalar, 1.0             estimation half width (Hz)
% searchHalfWid  scalar, 1.0             peak search half width (Hz); set
%                                         to [] to prevent fine-tuning of f
%
%                  <<< OUTPUT VARIABLES <<<
% NAME        TYPE/DEFAULT            DESCRIPTION
% d           same as input           cleaned raw data

% defaults
substHalfWid=0.2;
estHalfWid=1.0;
searchHalfWid=1.0;
pvpmod(varargin);

[n1 n2 n3]=size(d);
% sampling rate in Hz
fs=1e6/si;		
% number of frequencies to deal with
nF=numel(f);

% generate vector of frequencies of DFT series:
% - a is 0 for odd values, 1 for even
a=mod(n1+1,2);
% - freq values of one-sided spectrum
freq=fs/2*linspace(0,1,floor(n1/2)+1);
% - indexes to those for second half 
freqIx=floor(n1/2)+1-a:-1:2;
% append
freq=[freq freq(freqIx)];

% fft
ft=fft(d,[],1);
pha=angle(ft);
mag=abs(ft);

% means across both 2D and 3D
if n2==1 && n3==1
  magMn=mag;
else
  if n2>1
    magMn=mean(mag,2);
  end
  if n3>1
    magMn=mean(mag,3);
  end
end

% fine-tune frequencies to be eliminated by searching for maxima
if ~isempty(searchHalfWid) 
  fSearch=f(:)*[1 1]+repmat(searchHalfWid*[-1 1],nF,1);
  for g=1:nF
    ix=find(freq>=fSearch(g,1) & freq<=fSearch(g,2));
    % get rid of indexes to redundant frequencies in second half of array
    % (important, because due to limited precision the values of the FFT
    % transform in both halves are not exactly the same, so we may end up
    % with an index which is totally off)
    ix(ix>floor(n1/2)+1)=[];
    [~,mIx]=max(magMn(ix));
    f(g)=freq(ix(1)+mIx-1);
  end
end

% boundaries of frequencies to be 'eliminated' by interpolation
fSubst=f(:)*[1 1]+repmat(substHalfWid*[-1 1],nF,1);
fEst=f(:)*[1 1]+repmat(estHalfWid*[-1 1],nF,1);

for g=1:nF
  % - index to frequencies to be substituted
  fSubstIx=freq>=fSubst(g,1) & freq<=fSubst(g,2);
  fEstIx=(freq>=fEst(g,1) & freq<=fEst(g,2)) & ~fSubstIx;
  if ~any(fEstIx)
      error('fEstIx is all zeroes - estHalfWid will have to be widened')
  end
  mag(fSubstIx,:,:)=repmat(mean(mag(fEstIx,:,:)),[numel(find(fSubstIx)), 1, 1]);
end
% convert back
d=real(ifft(mag.*exp(1i*pha),[],1));