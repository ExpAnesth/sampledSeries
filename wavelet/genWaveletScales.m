function [scales,f,scale2Freq,scaleTickFun]=genWaveletScales(typeWavelet,freq,varargin)
% ** function scales=genWaveletScales(typeWavelet,freq,varargin)
% generates wavelet scales given a range of frequencies of interest, the
% type of wavelet and the frequency overlap (see below). Scales are
% determined according to Jordan et al, Rev.Sci.Instrum. 68:1484-1494,
% 1997 (all equation numbers listed here refer to this paper).
% Input variable overlap is optional and must be specified as
% parameter/value pairs, e.g. as in 
%          genWaveletScales('morl',[1 100],'overlap',.8)
%
%                         >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT      DESCRIPTION
% typeWavelet    char,'morl'       type of wavelet (currently, only Morlet)
% freq           array             if of the form [x1 x2] scales spanning 
%                                  frequencies x1 and x2 will be computed
%                                  if [x1 x2...xn] scales corresponding
%                                  to exact values will be computed
% overlap        scalar, .9        definition see eq. 34
%                     
%                         <<< OUTPUT VARIABLES <<<
%
% NAME           TYPE/DEFAULT       DESCRIPTION
% scales         array              the wavelet scales
% f              array              the corresponding frequencies (Hz)
% scale2Freq     scalar             factor by which 1/scales must be multi-
%                                   plied to obtain frequency values (Hz)
% scaleTickFun   function handle    compute e.g. tickVal=scaleTickFun([100 50 25 10])
%                                   to obtain ordinate values in plots when e.g.
%                                   imagesc(t,1:nWavelet,...) is used

overlap=.9;
pvpmod(varargin)

switch typeWavelet
  case 'morl';
    % - omega0 (sic!! see help for cwtft)
    w0=6.0;
    % - analytical conversion factor from scales to freq
    % *** don't use function centfrq with cwtft! ***
    scale2Freq=(w0+sqrt(2+w0^2))/(4*pi);
    % - define initial scale in dimensional time (s)
    a1=scale2Freq/freq(2);
    % - proportionality constant, eq. 36
    wd=-sqrt(-2*log(overlap));
    % - K, eq. 41
    K=w0/(w0+wd);
    if numel(freq)>2
      scales=scale2Freq./freq;
    elseif numel(freq)==2
      % - number of wavelets, given frequency resolution and initial scale (eq.
      % 43)
      nWavelet=ceil(log(freq(2)/freq(1))/log(K)+1);
      % - finally, the scales:
      scales=a1*K.^(0:nWavelet-1);
    else
      error('freq must contain two or more values');
    end
    f=scale2Freq./scales;
    scaleTickFun=@(x) 1+log((scale2Freq./x)/a1)./log(K);
  otherwise
    error('sorry, so far only ''morl'' implemented')
end
