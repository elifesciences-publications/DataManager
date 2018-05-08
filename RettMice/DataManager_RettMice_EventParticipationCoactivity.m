function DataManager_RettMice_EventParticipationCoactivity
%%identifying cell pairs with overlapping place fields
hf = gcbf; dbtype = getappdata(hf, 'dbtype'); celllistvar = getappdata(hf, 'celllistvar'); celllisttitle = getappdata(hf, 'celllisttitle');
%%%%%%%added options for all types of databases
ok = 1;
if ~strcmp(dbtype, '.spikedb')
    msgbox('Not a spikedb. Quit'); ok = 0;
end
if ok
    hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
    if ~isfield(pinfo, 'firing')
        msgbox('Firing properties not computed. Quit'); ok = 0;
    end
end
% if ok  %Are parameters needed?
%     inputIII=inputdlg({'Threshold for overlap (% of smaller field)'; 'New group name:'}, 'Pair Selection', 1, {'25'; 'FOverlap'}, 'on');
%     if (~isempty(inputIII))
%         [overlapthreshold, ok] = str2num(inputIII{1});
%     else
%         ok = 0;
%     end
% end
if ok
    hgroup = getappdata(hf, 'hgroup');hfield = getappdata(hf, 'hfield');
    plotparm = getappdata(hf, 'plotparm'); ow = plotparm.overwrite; %if ow=1, plot into the current database
    [writefilename, ok] = getoutputfile(hf, ow);
    %%%group selection
    groupselection = getappdata(hgroup, 'selection'); 
    %%%cell selection
    cellind = []; grpind = find(groupselection == 1); 
    if numel(grpind) == 0
        msgbox('No group(s) selected. Quit'); ok = 0;
    end
end
if ok
    for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
    pinfo = FindRippleSpikingProperties(pinfo,data,cellind);
    %%%update plots
    save(writefilename, 'pinfo', 'data');
    if (ow)
       iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
       fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
       newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
    else
        hmain = DataManager_DataManager_Callback; plotparm.linkspike = 1; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
    end
    setappdata(hmain,'plotparm', plotparm);
    DataManager_PlotSpikeDatabase(hmain, pinfo, data, 'Cell', 'clname', '.spikedb');
    set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
    setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); pinfo = []; data = [];
end

function pinfo = FindRippleSpikingProperties(pinfo,data,cellind)
%variable to assign
if (~isempty(cellind))
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo.firing, 'evtPartRate')) pinfo.firing.evtPartRate = cell(1, nspike); end
   if (~isfield(pinfo.firing, 'evtMeanNSpike')) pinfo.firing.evtMeanNSpike = cell(1, nspike); end
end
for (jjjk = 1:numel(cellind))
    i = cellind(jjjk); 
    %%%%load computing parameters 
    spiketime = sort(pinfo.parm.timeunit(i)*data.spike.spiketime{i}); [a,b] =size(spiketime); if a==1 spiketime = spiketime'; end %%%all column vectors
    evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i};        
    nev = numel(evName);
    pinfo.firing.evtParteRate{i} = NaN*ones(1, nev); pinfo.firing.evtMeanNSpike{i} = NaN*ones(1, nev);
    for (tt = 1:nev)
        [rate, meanN] = findpartrateNspike(spiketime, evTime{tt});
        pinfo.firing.evtPartRate{i}(tt) = rate; pinfo.firing.evtMeanNSpike{i}(tt) = meanN;
    end
end

function [rate, meanN] = findpartrateNspike(spiketime, ep)
rate = NaN; meanN = NaN;
if ~isempty(ep)
nev = numel(ep.start); Nspike = zeros(1, nev); partic = zeros(1, nev);
for kk = 1:nev
    idnow = find( (spiketime>=ep.start(kk)) & (spiketime<=ep.ent(kk)) );
    if (~isempty(idnow))
        Nspike(kk) = numel(idnow); partic(kk) = 1; 
    end
end
rate = numel(find(partic))/nev; meanN = mean(Nspike);
end

function [writefilename, okk] = getoutputfile(hf, ow)
okk = 1; writefilename = [];
if (ow == 0)
   [fname, pname] = uiputfile(fullfile(cd, '*.spikedb'), 'Write the new spike database to:');
   if (numel(fname)>1)
      writefilename = fullfile(pname, fname);
   else
      okk = 0;
   end
else
   %input = questdlg('The current database will be altered and overwritten. Are you sure?', 'Overwrite?', 'Yes');
   %if (strcmp(input, 'Yes'))
      fname = get(hf, 'Name'); ftt = strfind(fname, '__'); writefilename = fname(ftt+2:numel(fname));
   %else
   %   okk = 0;
   %end
end

