function [P,F,avP,varargout]=fspecp(d,si,varargin)
% ** function [P,F,avP,varargout]=fspecp(d,si,varargin)
%   computes power spectral density estimate of a sampled series
%   (real-valued) with the averaged modified periodogram method (division
%   of data into segments, windowing them & averaging ffts of these). The
%   psd is one-sided. Reason for using this function instead of the
%   matlab-provided functions:
%   (i) no need for the signal processing toolbox if method 'fft' is
%       specified
%   (ii) if a range of acceptable windows is specified, actual window  
%      length will optimize (=minimize) left-over data pieces
%   (iii) if method 'fft' is chosen amplitude, phase and power can be 
%      obtained for each segment individually (in variable output arguments)
% Optional input parameters listed below (= all except d and si) must be
% specified as parameter/value pairs, e.g. as in 
%   [P,F,avP]=fspecp(d,si,'win',[1000 2000]);
% *** Copy-and-paste line:
% [P,F,avP]=fspecp(d,si,'meth',meth,'win',win,'wintype',wintype,'olap',olap,'limFreq',limFreq);
%                    >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% d                1d array              sampled series
% si               scalar                sampling interval in us
% meth             string,'amp'          the method to use
%                                        'amp' - Averaged Modified Periodogram method going back to Welch
%                                        'psd' - same as above using psd.m
%                                        'fft' - same as above using fft directly
% win              scalar OR             ms, length of subsegment from which to compute single spectrum.
%                  2 element arr         A range of acceptable values may be specified. If omitted, win is
%                                        set to the length of sampled series d (which amounts to no averaging at all)
% wintype          char arr, 'hanning'   type of window; in addition to the matlab-provided windows 'rect' 
%                                        may be specified
% olap             scalar, 0             ms, overlap between segments
% limFreq          2element-arr,[0 60]   Hz, lower & upper limit of freqs of spectra
% detr             char arr, 'constant'  detrending method: 
%                                        'constant' - subtract mean
%                                        'linear' - subtract straight-line fit (maybe useful for intracellular data)
%
%                    <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT      DESCRIPTION
% P                column array      power spectral density in power per Hz (e.g. uV^2/Hz) 
% F                column array      frequencies corresponding to P
% avP              scalar            average power over total Nyquist interval
% varargout{1}     2D array          phase of segments
% varargout{2}     2D array          amplitude of segments
% varargout{3}     2D array          power of segments

% ************ default values *************************************
win=1000;
olap=0;
limFreq=[0 60];
meth='amp';
detr='constant';
wintype='hanning';
verbose=1;
varargout={};

% *************** check input arguments ***************************
pvpmod(varargin);
disp(['**** ' mfilename ':']);
if nargout > 3 && ~strcmpi(meth,'fft')
  error('segment-wise estimates of spectral values can only be obtained with method ''fft''');
end

% ********** preparatory stuff *************************************
% not all of it is needed for all methods
% sampling frequency in Hz
fs=1e6/si;
% moving window methods: 
% variables in points
win_pts=round(win/si*1e3);
olap_pts=round(olap/si*1e3);
dlen_pts=length(d);
% length of data segment in ms
dlen=dlen_pts*si/1000.0;

% ****************************************************************************
% attention: win_pts may be modified after call to mkintrvls
% intrvls_pts holds the start and stop times (in points) of the intervals in its 1st and 2nd column, respectively
% ****************************************************************************
[intrvls,intrvls_pts,win,win_pts]=mkintrvls([0 dlen],'resol',si*.001,'ilen',win,'olap',olap,'verbose',1);
% only after this point can variables depending on win be calculated 
if strcmpi(wintype,'rect'), window=ones(win_pts,1); 
else eval(['window=' wintype '(win_pts);']);
end
if verbose, disp(['frequency resolution without zero-padding: ' num2str(fs/win_pts) ' Hz']); end
% for length of window, find the next larger power of 2 
nfft_pts=2^(ceil(log2(win_pts)));
if (length(nfft_pts)-length(window))/length(nfft_pts)>0.5, disp('extensive zero-padding)'); end;
if verbose, disp(['frequency resolution after zero-padding: ' num2str(fs/nfft_pts) ' Hz']); end  
% detrend to get rid of frequency leakage at low frequencies 
d=detrend(d,detr);

switch meth
case 'amp'
  [P,F]=pwelch(d,window,olap_pts,nfft_pts,fs);
case 'psd'
  warning(['in the former version of ' mfilename  ' ''psd'' called pwelch.m while now it calls psd.m - make sure this is what you want']);
  [P,F]=psd(d,nfft_pts,fs,window,olap_pts,'none');
  % this scales P (factor 2 for one-sidedness):
  P=2*P(:,1)/fs;
case 'fft'
  nInts=size(intrvls_pts,1);
  % psd and pwelch automatically scale the window - here it has to be done manually
  % (condition: sum of squared elements of window divided by (scalfac*win_pts) must be==1)
  scalfac=win_pts/sum(window.^2);
  % ********** do it  *************************************
  % P is the (complex) fft, one column per interval, F the corresponding frequencies
  for i=nInts:-1:1
    % pick interval, apply window
    tmpd=d(intrvls_pts(i,1):intrvls_pts(i,2)).*window;
    P(:,i)=fft(tmpd,nfft_pts);
  end
  % first thing to do: cut down P to half+1 points (one-sided fft)
  P=P(1:nfft_pts/2,:);
  % corresponding frequency vector
  F=fs/nfft_pts*(0:nfft_pts/2-1)';
  % compute magnitude and phase (unwrapped) for individual segments (but only 
  % if requested because it is computationally intensive)
  % - multiply by scaling factor for window
  % (@ scaling fac for magnitude not yet done, probably different scaling factor)
  if nargin>=4
    varargout{1}=angle(P);
  end
  if nargin>=5
    varargout{2}=abs(P);
  end
  % compute power (individual segments) & scale to psd:
  % - multiply by 2 because we need the one-sided power spectrum
  % - multiply by scaling factor for window, see above
  % - divide by sampling freq fs and number of points in original (NOT zero-padded) data segment
  P=scalfac*2*(P.*conj(P))/win_pts*si*1e-6;
  % power of individual segments
  if nargin>=6
    varargout{3}=P;
  end
  % averaged power
  P=mean(P,2);
  
otherwise
  error('no valid method for power spectral density estimation specified');
end;
% to obtain power spectrum multiply by fs; to obtain average power over 
% entire Nyquist interval multiply by fs and divide by (2*number of bins
% PRIOR TO ZERO-PADDING)
% (*2 because in one-sided psds the total power is contained in half the 
% points a two-sided psd is composed of)
% So, average power=
avP=sum(P)*fs/(2*win_pts);

% on the occasion, to obtain the power in a specific frequency band from
% the PSD: 
% - multiply the PSD by fs and average the values in the bins in
% question
% which amounts to the same as 
% - multiplying the PSD with the bin width and summing the bins in question

% ** concluding section **
% cut down results to desired frequency range
limFreqInd=find(F(:,1)<=limFreq(2) & F(:,1)>=limFreq(1));
P=P(limFreqInd,:);
F=F(limFreqInd);
for g=1:length(varargout)
  varargout{g}=varargout{g}(limFreqInd,:);
end


% ------------------------------------------------------------------------------------------
% ** A short compendium of some matlab methods for power spectrum estimates ****
% a) nonparametric
% Generally, all functions 
% - are capable of generating plots
% - allow specification of a window (or even do so automatically) and scale it properly
% - return one-sided psd for real input (equal to 2* two-sided psd evaluated at positive frequencies)
%   and two-sided psd for complex input. THIS DOES NOT APPLY TO psd.m!
%
% psd.m is a little funny - there is no help in the help window, the code looks a little
% fussed-with by several programmers and it does no scaling but detrending.
% 
% function name  units/scaling                       averaging of segments?      detrending possible
% periodogram    power per Hz (one-and twosided)     no!                         nay
%                OR power per radians per sample 
%                (fft divided by 2pi or fs)
% 
% psd            none, must divide by fs and         yes                         yes
%                multiply by 2 (two-sided)
% pwelch         power per Hz (one-and twosided)     yes                         no
%                OR  power per radians per sample 
%                (fft divided by 2pi or fs)
%                
% pmtm           - multitaper, not yet dealt with -


% note below is taken from psd.m:

%   NOTE 2: The Power Spectral Density of a continuous-time signal,
%           Pss (watts/Hz), is proportional to the Power Spectral 
%           Density of the sampled discrete-time signal, Pxx, by Ts
%           (sampling period). [2] 
%       
%               Pss(w/Ts) = Pxx(w)*Ts,    |w| < pi; where w = 2*pi*f*Ts

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice hall, 1997, pg, 15
%     [2] A.V. Oppenheim and R.W. Schafer, Discrete-Time Signal
%         Processing, Prentice-Hall, 1989, pg. 731
%     [3] A.V. Oppenheim and R.W. Schafer, Digital Signal
%         Processing, Prentice-Hall, 1975, pg. 556




% trash ------------------------------ trash ------------------------------ trash 

%   % cut down results to desired frequency range - doing this here before calculating power saves ressources
%   limFreqInd=find(F(:,1)<=limFreq(2) & F(:,1)>=limFreq(1));
  % compute power & average & scale to psd:
  % - multiply by 2 because we need the one-sided power spectrum
  % - multiply by scaling factor for window, see above
  % - divide by sampling freq fs and number of points in original (NOT zero-padded) data segment
%   P=scalfac*2*mean(abs(P(limFreqInd,:)).^2,2)/win_pts/fs;
%   F=F(limFreqInd);

