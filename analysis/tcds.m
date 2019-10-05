function tsl=tcds(d,si,thresh,varargin)
% ** function tsl_=tcd(d,si,thresh,varargin)
% 'threshold crossing detector single' - a threshold-based spike detector
% producing time stamp lists of single-column arrays. 
%
%                    >>> INPUT VARIABLES >>>
% NAME       TYPE/DEFAULT               DESCRIPTION
% d          column array               raw data
% si         a) scalar                  sampling interval in µs
%            b) string 'idx'            sampling interval is set to 1; time
%                                       stamps will be given as indices
% thresh     scalar                     threshold; if negative, assume that 
%                                       negative-going crossing is looked
%                                       for
% detMode    char, 'cross'              - 'cross': time stamp is set to
%                                       crossing of threshold
%                                       - 'peak': each peak beyond
%                                       threshold will be assigned a time
%                                       stamp
%                  <<< OUTPUT VARIABLES <<<
%
% NAME       TYPE/DEFAULT               DESCRIPTION
% tsl        1d arr                     time stamp list (ms or indices)
%
% See also: tcd.m, the precursor function, which operates on up to 3D
% arrays.

% (c) Harald Hentschke 2019

detMode='cross';
pvpmod(varargin,{'detMode'});

validateattributes(d, {'numeric'}, {'column'})

if ischar(si)
    if strcmpi(si,'idx')
        si=1;
    else
        error('si must be a scalar or string ''idx''');
    end
else
    % convert to ms
    si=si*.001;
end

switch detMode
    case 'cross'
        % if thresh is negative, assume that negative-going crossing is looked for
        if thresh<0.0
            tc=d<=thresh;
        else
            tc=d>=thresh;
        end
        tc=diff(tc);
        % thus, first tick above thresh is defined as time of crossing
        tsl=(find(tc==1)+1)*si;
    case 'peak'
        % dd: positive-going peaks = -1, negative-going peaks = 1
        dd=diff(diff(d)>0);
        % if thresh is negative, assume that negative-going crossing is looked for
        if thresh<0.0
            % indexes to negative-going peaks
            ix=find(dd==1)+1;
            % pick those below threshold
            tsl=ix(d(ix)<thresh)*si;
        else
            % indexes to positive-going peaks
            ix=find(dd==-1)+1;
            % pick those above threshold
            tsl=ix(d(ix)>thresh)*si;
        end
    otherwise
        error('bad detMode')
end
