function d=elim_artefact(d,si,varargin)
% ** function d=elim_artefact(d,si,varargin)
% is a collection of routines for detecting and removing artifacts in
% sampled series d (electrophysiological data), based on typical properties
% of artifacts such as high amplitude, high slope, or one-sample spikes. 
% Input parameters except d, and si are optional and must be specified
% as parameter/value pairs, e.g. as in
%          d=elim_artefact(d,100,'thresh',5);
%
%    -- UNDER CONSTRUCTION --
%
%                         >>> INPUT VARIABLES >>>
% NAME       TYPE/DEFAULT        DESCRIPTION
% d          array               time series data (time down columns)
% si         scalar              sampling interval in µs
% afType     char, 'needle-1s'   type of artifact:
%                                 - 'needle-1s': 1-sample 'needles' 
%                                 - 'hi-amp': large artifacts with much
%                                 higher amplitude than neuronal signals
%                                 - 'needle': spikes best removed by hampel
%                                   filter
% thresh     array, [3.5 2.5]    threshold(s): exact meaning depends on
%                                 afType
% sIntv      2-element arr       peri-artifact substitution interval (ms; 
%                                 see etslexcsubst.m)
% eIntv      2-element arr       peri-event estimation interval (ms) from 
%                                 which substitute data will be estimated;
%                                 see etslexcsubst.m
% neighHood  scalar, 1           number of neighbors for Hampel filter
%                                 (specified in ms, will be converted to
%                                 points; see hampel.m)
% doPlot     logical, false      if true, original and corrected data trace 
%                                 are plotted together in a single graph
%
%                         <<< OUTPUT VARIABLES <<<
%
% NAME           TYPE/DEFAULT           DESCRIPTION
%
%

% to do 
% - help, precise explanation of vars

afType='needle-1s';
thresh=[3.5 2.5];
sIntv=[];
eIntv=[];
neighHood=1;
doPlot=false;
pvpmod(varargin)

% some 'constants'
numPlotCutout=500;

switch afType
  case 'needle-1s'
    if numel(thresh)<2
      error('for artifacts of type ''needle-1s'' thresh must contain 2 values');
    end
  case 'hiAmp'
  case 'needle'
    if numel(thresh)~=1
      error('for artifacts of type ''needle'' thresh must contain exactly one value');
    end
    numNeigh=cont2discrete(neighHood,si/1000,'intv',1);
  otherwise
    error('undefined artifact type')
end

[n1,n2,n3]=size(d);
% catch a typical problem
if n3==1 && n1==1 && n2>n1
  warning('d is a row array - reshaping');
  d=d';
  [n1,n2,n3]=size(d);
end

if doPlot
  fh=figure;
  sphArr=gobjects(2,1);
  sphArr(1)=subplot(3,1,1);
  nPlotPts=min(1e6,n1)-1;
  t=(0:nPlotPts-1)'*(si/1000);  
  % plot only first column, otherwise it'll be too crowded
  plot(t,d(1:nPlotPts,1),'m');
  hold on
end

cutout=[];
sumNAf=0;
switch afType
  case 'needle-1s'
    % ** artifacts consisting of a single stray sample point **
    % differentiate, order one, down columns
    diffD=diff(d,1,1);
    % determine mean and sd
    mn=mean(diffD);
    sd=std(diffD);
    % *** be aware that this loops through all columns of what may be a 3D
    % array via 'linear column indexing' ***
    for colIx=1:n2*n3
      % find indices to points whose derivative is above...
      posIx=diffD(:,colIx)-mn(colIx)>sd(colIx)*thresh(1);
      % ...or below thresh(1)*sd; these are the potential offenders
      negIx=diffD(:,colIx)-mn(colIx)<-sd(colIx)*thresh(1);
      % find indices to points whose derivative is beyond thresh(2)*sd in
      % either direction (see below)
      beyondIx=diffD(:,colIx)-mn(colIx)>sd(colIx)*thresh(2) | ...
        diffD(:,colIx)-mn(colIx)<-sd(colIx)*thresh(2);
      % artifacts are two adjacent points with extreme and opposite slopes
      % (=abs of slopes of both consecutive points are beyond
      % thresh(1)*sd), flanked by points without extreme slopes (whose abs
      % slopes are below thresh(2)*sd)
      afIx=((posIx & [false; negIx(1:end-1)]) | (negIx & [false; posIx(1:end-1)])) ...
        & ~[false; false; beyondIx(1:end-2)] & ~[beyondIx(3:end); false; false];
      afIx=find(afIx);
      % by definition, the artifacts detected here cannot be right at the
      % beginning or end of the sampled series, so don't worry about border
      % effects
      nAf=numel(afIx);
      sumNAf=sumNAf+nAf;
      if doPlot && size(cutout,2)<numPlotCutout && colIx==1 
        cutout=cat(2,cutout,tsl2exc(d(:,colIx),'idx',{afIx},'win',[-10 10]));
      end
      % now substitute
      for ai=1:nAf
        d(afIx(ai),colIx)=mean(d(afIx(ai)+[-1 1],colIx));
      end
    end
    
  case 'hiAmp'
    % !! CASE NOT TESTED YET!!
    % *** be aware that this loops through all columns of what may be a 3D
    % array via 'linear column indexing' ***
    for colIx=1:n2*n3
      % detect artifacts as 'bursts': a critical parameter is the
      % time interval below which two threshold-surpassing portions
      % of data (=signal excursions of a biphasic artifact) will be
      % merged and treated as one. For now, assume that the signal
      % excursions of an artifact are separated by 1 ms at most
      afEtsl=tbt(abs(d(:,colIx)),'idx',thresh,...
        'minActive',0,'minInactive',round(1000/si),...
        'elimOrder','inactive');
      nAf=size(afEtsl,1);  
      sumNAf=sumNAf+nAf;
      if doPlot && size(cutout,2)<numPlotCutout && colIx==1 
        cutout=cat(2,cutout,tsl2exc(d(:,colIx),'idx',{afEtsl},'win',[-10 10]));
      end
      % substitution of values:
      d(:,colIx)=etslexcsubst(d(:,colIx),'idx',afEtsl,...
        cont2discrete(sIntv,si/1000),cont2discrete(eIntv,si/1000));
    end
    
  case 'needle'
    for colIx=1:n2*n3
      [d(:,colIx),ix]=hampel(d(:,colIx),numNeigh,thresh);
      afIx=find(ix);
      % to determine the locations of the artifacts, combine
      afEtsl=etslburstf(afIx,0.5,'recLen',afIx(end)+1);
      nAf=size(afEtsl,1);  
      sumNAf=sumNAf+nAf;
      if doPlot && size(cutout,2)<numPlotCutout && colIx==1
        cutout=cat(2,cutout,tsl2exc(d(:,colIx),'idx',{afEtsl},'win',[-10 10]));
      end
    end
end

disp([mfilename ': eliminated ' int2str(sumNAf) ' artifacts of type ' afType]);


if doPlot
  subplot(sphArr(1))
  plot(t,d(1:nPlotPts,1),'k');
  axis tight
  grid on
  if strcmp(afType,'needle-1s')
    sphArr(2)=subplot(3,1,2);
    plot(t,diffD(1:nPlotPts,1));
    hold on
    % plot first threshold (of first channel only)
    line(t([1,nPlotPts])*[1 1],(mn(1)+sd(1)*thresh(1)*[1;1]*[-1 1]),...
      'linestyle','--','color','r')
    % plot second threshold (of first channel only)
    line(t([1,nPlotPts])*[1 1],(mn(1)+sd(1)*thresh(2)*[1;1]*[-1 1]),...
      'linestyle',':','color','k')
    axis tight
    grid on
    xlabel('time (ms)')
    linkaxes(sphArr,'x');
  end
  % plot cutouts of eliminated artifacts
  subplot(3,2,5)
  plot(cutout(:,1:min(numPlotCutout,sumNAf)),'m')
  xlabel('sample')
  axis tight
end
