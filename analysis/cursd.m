function c=cursd(d,varargin)
% ** function c=cursd(d,varargin)
% computes current source density from data d (very simple difference
% algorithm, no smoothing).
% All input parameters listed below except d must be specified as
% parameter/value pairs, e.g. as in 
%       c=cursd(d,'dx',0.1);
% 
%                         >>> INPUT VARIABLES >>>
% NAME              TYPE/DEFAULT          DESCRIPTION
% d                 2d array              data (time runing down columns)
%                                         (unit: mV)
% dx                scalar                spacing between recording sites 
%                                         (unit: mm)
% imp               scalar                (assumed) impedance of extracellular medium
%                                         (unit: Ohm*mm)
%    ** if BOTH dx and imp are provided, output c will be scaled (unit: mA/mm^3)
%
%                         <<< OUTPUT VARIABLES <<<
% NAME             TYPE/DEFAULT           DESCRIPTION
% c                2d array               current source density
%

% ----- default values & varargin -----
dx=[];
imp=[];
pvpmod(varargin);


[n1 n2 n3]=size(d);
% checks
if n2<3
  error('a minimum of 3 voltage traces must be given')
end

% preallocate
c=repmat(0,[n1 n2-2 n3]);
% do it 
for ix=2:n2-1
  c(:,ix-1)=2*d(:,ix)-d(:,ix-1)-d(:,ix+1);
end
if length(dx)==1 && length(imp)==1
  c=c/dx^2/imp;
end

% to obtain real units divide c by square of dx and impedance of
% extracellular medium (unit of latter: ohm*mm, one estimate is 300 ohm*cm
% = 3000 ohm*mm = 30000 ohm*um): (mV/mm^2)/(Ohm*mm) = mA/mm^3
% See Brankack et al '93