function w=awfm(d,fs,wf)
% ** function w=awfm(d,fs,wf)
% 'arbitrary waveform frequency modulation' does principally the same as
% vco.m, the difference being that an arbitrary waveform instead of a sine
% will be modulated in frequency 
%  Notes: 
% 1. this is a potentially memory-demanding task. If the load is too high, 
%  split d up into smaller chunks
% 2. no safeguard whatsoever against nonsense values in d (especially 
%  negative values which could result in negative frequencies)
% 
%                    >>> INPUT VARIABLES >>>
% NAME    TYPE/DEFAULT          DESCRIPTION
% d       2D or 3D array        the frequency-modulating data (time running 
%                                down the columns)
% fs      scalar                sampling frequency (Hz) of d AND output w
%                     
%                    <<< OUTPUT VARIABLES <<<
% NAME     TYPE/DEFAULT         DESCRIPTION
% w        array (same as d)    output waveform


% ----- default values & varargin -----

% number of points in waveform
wnpts=length(wf);
% this is the carrier frequency (at which the original waveform could be
% put out by just playing a concatenation of it with fs)
fc=fs/wnpts;
% map d to fc: scale [0 1] in d to [0 fc]
d=d*fc;
% ..and calculate indices into wv
d=mod(round(cumsum(d/fc)),wnpts)+1;

% voila
w=wf(d);

