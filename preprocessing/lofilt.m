function [d, newsi, cfreq, ord]=lofilt(d, si, b, a, varargin)
% ** function [d, newsi, cfreq, ord]=lofilt(d, si, b, a, varargin)
%   Passes data d through butterworth lowpass filter using double reverse 
%   filtering routine 'filtfilt'. Double reverse filtering has the advantage 
%   of zero phase shift of the resulting signal. 
% 
%                         >>> INPUT VARIABLES >>>
% NAME              TYPE/DEFAULT         DESCRIPTION
% d                 array                data to be filtered (along columns!)
% si                scalar               sampling interval in us
% b                 array                filter coefficients 
% a                 array                filter coefficients (**Note: if b and a are empty,
%                                         they will be computed from cfreq, rs, and rp)
% cfreq             scalar               cutoff (also termed corner, or -3 dB) frequency 
%                                         of filter in Hz
% rs                scalar, 30           minimal attenuation in stopband (dB) - see buttord
%                                         also, minimal steepness of attenuation (dB/octave) 
%                                         for frequencies between passband and stopband
% rp                scalar, 0.5          max allowed ripples in passband (dB) - see buttord
% pickf             scalar (int), 1      downsampling factor
% verbose           scalar, 1            if 1, numerical results of the calculations will be 
%                                         printed on screen. If any other nonzero value, 
%                                         in addition the filter's frequency response will 
%                                         be plotted
%
%                         <<< OUTPUT VARIABLES <<<
% NAME              TYPE/DEFAULT         DESCRIPTION
% lfd               1d- or 2d-array      filtered data 
% newsi             scalar               the sampling interval of domnsampled data in us
% cfreq             scalar               cutoff frequency in Hz resulting from computations
% ord               scalar               length of filter parameter arrays (order of filter)


% default values
rs=30;
rp=0.5;
pickf=1;
verbose=0;

pvpmod(varargin);
if verbose
    disp(['**** ' mfilename ':']);
end

% compute filter parameters if they were not specified
if isempty(a) || isempty(b)
    [a, b, cfreq, ord] = lofi_par(si, cfreq, 'rs', rs, 'rp', rp, 'verbose', verbose);
end
% filter
d=filtfilt(b,a,d);
newsi=si;
% downsample?
if pickf > 1
    p = 1:pickf:size(d,1); 
    d = d(p,:,:);
    newsi = pickf*si;
    if verbose
        disp(['New sampling interval: ' num2str(newsi) ' µs']);
    end
end