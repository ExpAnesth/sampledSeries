function [P,F,avP,krnlfft]=ac2fspecp(ac,resol,varargin)
% ** function [P,F,avP,krnlfft]=ac2fspecp(ac,resol,varargin)
%  computes power spectral density (psd) estimate from an autocorrelation (AC) sequence.
%  This function is most useful for obtaining spectral estimates of spike trains:
%  simply compute their unbiased autocorrelation sequence (see tslcc.m for an 
%  explanation) and put it in this function.
%  ** NOTES **
%  (A) The frequency resolution of P is 
%         [sampling frequency] divided by [total number of points in AC sequence];
%      thus both the length of 'ac' and the value of 'resol' determine the spectral resolution
%  (B) The windows: 
%      windowing can be done in two stages, on the AC and on the fft of the AC. 
%      1. The first ('wintype') serves to reduce the variance of the spectrum
%      (see Oppenheim '99, p. 744) stemming from the fact that border bins in AC always 
%      have a larger variance than more central bins. 
%      2. The second ('muetze') will be applied to the Discrete Fourier Transform (DFT) of the AC;
%      this window must itself be the (symmetric) DFT of a 'kernel' function, like a Gaussian or 
%      triangular window. If windowing is done in this place, the spectrum is the same as the one 
%      that would be obtained from a DFT of the convolution of the original series (spike train) 
%      with this kernel function. 
% 
%                    >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% ac               1d array              AC sequence
% resol            scalar                the bin width (sampling interval) of AC in ms
% wintype          char array, 'rect'    window to apply to AC (size is given by length of AC)
%                                        'rect' - rectangular window
%                                        'triang' - triangular (or Bartlett) window
% muetze             char array, 'none'  fft of this window type (or 'kernel') will be multiplied 
%                                        with the fft of the AC sequence; choose among 
%                                        'none', 'triang' and 'gauss'
% sigma            scalar, 10            the standard width of the kernel (ms)
% limFreq          2element-arr,[0 60]   lower & upper limit of freqs of spectra (Hz) 
%
%                    <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT      DESCRIPTION
% P                column array      power spectral density in power per Hz (e.g. uV^2/Hz) 
% F                column array      vector of frequencies corresponding to P
% avP              scalar            average power over total Nyquist interval
% krnlfft          2-col array       1st column: frequency, 2nd column: abs(fft(kernel))
%                                    (useful if you want to evaluate choice of kernel)

wintype='rect';
limFreq=[0 60];
muetze='none';
sigma=10;
verbose=1;

pvpmod(varargin);
disp(['**** ' mfilename ':']);

% ********** preparatory stuff *************************************
ac_pts=length(ac);
krnlfft=[];

if ~mod(ac_pts,2)
  error('ac must have an uneven number of elements');
end

switch wintype
case 'rect'
  window=ones(ac_pts,1);
  scalfac=1;
case 'triang'
  window=triang(ac_pts);
  scalfac=ac_pts/sum(window.^2);  
otherwise
  error('no valid window specified');
end;

% apply window to AC
ac=ac.*window;
% fft
P=fft(ac);

% kernel funcs - must be column arrays
if ~strcmpi(muetze,'none')
  switch muetze
    case 'gauss'
      % total half-length of krnl in ticks
      mHalfLen_pts=(ac_pts-1)/2;
      krnl=(exp(-.5*(([-mHalfLen_pts:mHalfLen_pts]*resol).^2)/sigma^2))';
    case 'triang'
      mHalfLen=sqrt(6)*sigma;
      % same in ticks
      mHalfLen_pts=round(mHalfLen/resol);
      kitt=zeros((ac_pts-(2*mHalfLen_pts+1))/2,1);
      krnl=[ kitt; triang(2*mHalfLen_pts+1); kitt];
    otherwise
      error('no valid muetze specified');
  end;
  % scale krnl such that area==1
  krnl=krnl/sum(krnl);
  krnlfft=fft(krnl);
  P=P.*krnlfft;
  % after job is done 1. make one-sided 2. pre-pend freq 
  krnlfft=[makecol(1000*[0:(ac_pts-1)/2-1]/(resol*ac_pts)) abs(krnlfft(1:(ac_pts-1)/2))];    
end

% power spectral density = magnitude of (windowed) AC sequence divided by sampling freq
P=scalfac*abs(P)*resol*.001;
% to obtain power spectrum multiply by fs; to obtain average power over 
% entire Nyquist interval multiply by fs and divide by (2*number of bins)
% (*2 because in one-sided psds the total power is contained in half the 
% points a two-sided psd is composed of)
% So, average power=
avP=sum(P)*1000/resol/(2*length(P));

% frequency resolution is fs/number of pts, here expressed as Hz
F=((0:length(P)-1)*(1/resol)*(1/ac_pts)*1000.0)';
if verbose, disp(['frequency resolution: ' num2str(diff(F(1:2))) ' Hz']); end
% cut down results to desired frequency range
limFreqInd=find(F(:,1)<=limFreq(2) & F(:,1)>=limFreq(1));
P=P(limFreqInd,:);
F=F(limFreqInd);

