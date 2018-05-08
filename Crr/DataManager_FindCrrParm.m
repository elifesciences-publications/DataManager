function [pinfo, data] = DataManager_FindCrrParm(pinfo, data)
%%set up all parameters for later on calculations
timebinsize = 0.01; maxlag = 100*timebinsize; onesideintbin = 5;
firstPsearchwin = [-25 25]*timebinsize;
secondPsearchwin = [-25 25]*timebinsize;

nspike = numel(pinfo.general.animalname);
pinfo.parm.timeunit = ones(1, nspike);   %timeunit = 1s
pinfo.parm.timebin = timebinsize*ones(1, nspike);     %time binsize in second
pinfo.parm.maxlag = maxlag*ones(1, nspike);     %maximum timelag in second
pinfo.parm.setbacktime = 2*ones(1, nspike);     %setbacktime added at the begining/end of the event times in second
pinfo.parm.intbin = cell(1, nspike);  %time bins to integrate for computing correlations 
pinfo.parm.crrmode = cell(1, nspike);     %mode to compute crr = either of rate, count, normcount

pinfo.parm.sesskeyword = cell(1, nspike);     %session to select for computing crr
pinfo.parm.sesskeytype = cell(1, nspike);     %session to select for computing crr
pinfo.parm.evtkeyword = cell(1, nspike);     %event to select for computing crr
pinfo.parm.evtkeytype = cell(1, nspike);     %event to select for computing crr

pinfo.parm.smoothmode = cell(1, nspike);     %mode to compute crr = either of rate, count, normcount
pinfo.parm.smoothbin = 8*ones(1, nspike);

pinfo.parm.sessType = cell(1, nspike);
pinfo.parm.eventtype = cell(1,nspike);
%%specify parameters for particular analysis (samples below):

for (i = 1:nspike)
    pinfo.parm.intbin{i} = [-1 1]*onesideintbin*timebinsize;
    
    pinfo.parm.P1stSearchWin{i} = firstPsearchwin; pinfo.parm.P2ndSearchWin{i} = secondPsearchwin; 
    pinfo.parm.searchPmode{i} = 'peak'; %%%%%%either 'peak', 'trough', or 'both'
    
    pinfo.parm.crrmode{i} = 'rate'; %%%rate, count or normcount
    pinfo.parm.smoothmode{i} = 'no';
    
    pinfo.parm.sesskeyword{i} = 'track';
    pinfo.parm.sesskeytype{i} = 'linear';
    pinfo.parm.evtkeyword{i} = 'sleep';
    pinfo.parm.evtkeytype{i} = 'sws';
    
    sessnow = pinfo.general.sessionname{i}; 
    for (j = 1:numel(sessnow))
        pinfo.parm.sessType{i}{j} = findsesstype(sessnow{j});
    end
    evn = pinfo.general.eventname{i};
    for (j = 1:numel(evn))
        pinfo.parm.eventtype{i}{j} = findeventtype(evn{j});
    end
end

function sessType = findsesstype(sessnow)
sessType = 'others';
if (~isempty(strfind(sessnow, 'open'))) | (~isempty(strfind(sessnow, 'field'))) | (~isempty(strfind(sessnow, 'OF'))) ...
        | (~isempty(strfind(sessnow, 'OP'))) | (~isempty(strfind(sessnow, 'platform'))) | (~isempty(strfind(sessnow, 'circ')))...
        | (~isempty(strfind(sessnow, 'square'))) | (~isempty(strfind(sessnow, 'Open')))...
        | (~isempty(strfind(sessnow, 'QS'))) | (~isempty(strfind(sessnow, 'QW'))) | (~isempty(strfind(sessnow, 'QN')))...
        | (~isempty(strfind(sessnow, 'QE'))) | (~isempty(strfind(sessnow, 'QX')))
    sessType = 'open';
elseif (~isempty(strfind(sessnow, 'track'))) | (~isempty(strfind(sessnow, 'linear'))) | (~isempty(strfind(sessnow, 'LT')))...
        | (~isempty(strfind(sessnow, 'Linear'))) | (~isempty(strfind(sessnow, 'Track'))) | (~isempty(strfind(sessnow, 'TRK')))
    sessType = 'linear';
elseif (~isempty(strfind(sessnow, 'sleep'))) | (~isempty(strfind(sessnow, 'Sleep'))) | (~isempty(strfind(sessnow, 'rest')))...
        | (~isempty(strfind(sessnow, 'Rest'))) | (~isempty(strfind(sessnow, 'SLP')))
    sessType = 'sleep';
end

function evType = findeventtype(evn)
evType = 'others'; 
if ( (~isempty(strfind(evn, 'cw'))) | (~isempty(strfind(evn, 'acw'))) | (~isempty(strfind(evn, 'CW'))) | (~isempty(strfind(evn, 'ACW')))...
        | (~isempty(strfind(evn, 'leftright'))) | (~isempty(strfind(evn, 'rightleft'))) | (~isempty(strfind(evn, 'Leftright')))...
        | (~isempty(strfind(evn, 'Rightleft'))) | (~isempty(strfind(evn, 'ncw'))) | (~isempty(strfind(evn, 'NCW')))...
        | (~isempty(strfind(evn, 'RL'))) | (~isempty(strfind(evn, 'LR'))) | (~isempty(strfind(evn, 'run')))...
        | (~isempty(strfind(evn, 'lefttoright'))) | (~isempty(strfind(evn, 'righttoleft'))) )...
    & ( (isempty(strfind(evn, 'ripple'))) & (isempty(strfind(evn, 'spindle'))) & (isempty(strfind(evn, 'full'))) ...
        & (isempty(strfind(evn, 'session'))) )
   evType = 'run';
elseif (~isempty(strfind(evn, 'stop'))) | (~isempty(strfind(evn, 'Stop'))) | (~isempty(strfind(evn, 'eat'))) | (~isempty(strfind(evn, 'Eat')))...
        | (~isempty(strfind(evn, 'leftleft'))) | (~isempty(strfind(evn, 'rightright'))) | (~isempty(strfind(evn, 'Leftleft')))...
        | (~isempty(strfind(evn, 'Rightright')))
    evType = 'stop';
% elseif (~isempty(strfind(evn, 'sleep'))) | (~isempty(strfind(evn, 'Sleep')))...
%         | (~isempty(strfind(evn, 'rest'))) | (~isempty(strfind(evn, 'Rest')))
%     evType = 'sleep';
elseif (~isempty(strfind(evn, 'ripp')))
    evType = 'ripple';
elseif (~isempty(strfind(evn, 'REM'))) | (~isempty(strfind(evn, 'rem')))
    evType = 'rem';
elseif (~isempty(strfind(evn, 'SWS'))) | (~isempty(strfind(evn, 'sws')))
    evType = 'sws';    
end