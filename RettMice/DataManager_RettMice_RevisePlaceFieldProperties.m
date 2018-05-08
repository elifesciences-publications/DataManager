function DataManager_RettMice_RevisePlaceFieldProperties
%%identifying cell pairs with overlapping place fields
hf = gcbf; dbtype = getappdata(hf, 'dbtype'); 
%%%%%%%added options for all types of databases
ok = 1;
if ~strcmp(dbtype, '.spikedb')
    msgbox('Not a spikedb. Quit'); ok = 0;
end
if ok
    hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
    if ~isfield(pinfo, 'field')
        msgbox('Place fields not computed. Quit'); ok = 0;
    end
end
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
    pinfo = RevisePlaceFieldProperties(pinfo,data,cellind);
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
disp('*********************');

function pinfo = RevisePlaceFieldProperties(pinfo,data,cellind)
%variable to assign
if (~isempty(cellind))
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo.field, 'TrajName')) pinfo.field.TrajName = cell(1, nspike); end
   if (~isfield(pinfo.field, 'TrajNField')) pinfo.field.TrajNField = cell(1, nspike); end
   if (~isfield(pinfo.field, 'TrajMeanRate')) pinfo.field.TrajMeanRate = cell(1, nspike); end
   if (~isfield(pinfo.field, 'AvgTrajName')) pinfo.field.AvgTrajName = cell(1, nspike); end
   if (~isfield(pinfo.field, 'AvgTrajSI')) pinfo.field.AvgTrajSI = cell(1, nspike); end
end
for (jjjk = 1:numel(cellind))
    i = cellind(jjjk); 
    disp(['-----> Revise place field measures (', num2str(jjjk), ' out of ', num2str(numel(cellind)), '): ', pinfo.general.parmfile{i}]);
    %%%% compute average SI across sessions
    iev = find( (pinfo.field.runSptlInfo{i}>-10) & (strcmp(pinfo.parm.eventtype{i}, 'run')) & (pinfo.firing.evtmeanrate{i}>0.5) );
    nev= numel(iev); evName = pinfo.general.eventname{i}(iev); evSI = pinfo.field.runSptlInfo{i}(iev);
    ss = cell(1, nev); ee = cell(1, nev);
    for (tt = 1:nev) [ss{tt}, ee{tt}] = strtok(evName{tt}, '_'); end
    [trajname, III] = unique(ee); pinfo.field.AvgTrajName{i} = trajname; pinfo.field.AvgTrajSI{i} = NaN*ones(numel(III),1);
    for (tt = 1:numel(III))
        jj = find(strcmp(ee, trajname{tt})); pinfo.field.AvgTrajSI{i}(tt) = mean(evSI(jj));
    end
    %%%%First, need to remove the fields on trajs with meanrate < 0.5
    fieldtraj = pinfo.field.PF1Devt{i}; ileft = zeros(1, numel(fieldtraj));
    AllevName = pinfo.general.eventname{i}; Allevrate = pinfo.firing.evtmeanrate{i};
    for (tt=1:numel(fieldtraj))
        III = find(strcmp(AllevName, fieldtraj{tt}));
        if (Allevrate(III) >= 0.5) ileft(tt) = 1; end
    end
    JJJ = find(ileft);
    pinfo.field.PF1Devt{i} = pinfo.field.PF1Devt{i}(JJJ); pinfo.field.PF1DNfield{i} = numel(JJJ);
    pinfo.field.PF1DBoundStart{i} = pinfo.field.PF1DBoundStart{i}(JJJ); pinfo.field.PF1DBoundEnd{i} = pinfo.field.PF1DBoundEnd{i}(JJJ);
    pinfo.field.PF1DLocPeakX{i} = pinfo.field.PF1DLocPeakX{i}(JJJ); pinfo.field.PF1DLocComX{i} = pinfo.field.PF1DLocComX{i}(JJJ);
    pinfo.field.PF1DLocStartX{i} = pinfo.field.PF1DLocStartX{i}(JJJ); pinfo.field.PF1DLocEndX{i} = pinfo.field.PF1DLocEndX{i}(JJJ);
    pinfo.field.PF1DLength{i} = pinfo.field.PF1DLength{i}(JJJ); pinfo.field.PF1DSize{i} = pinfo.field.PF1DSize{i}(JJJ);
    pinfo.field.PF1DSkewness{i} = pinfo.field.PF1DSkewness{i}(JJJ); pinfo.field.PF1DKurtosis{i} = pinfo.field.PF1DKurtosis{i}(JJJ);
    pinfo.field.PF1DInMeanrate{i} = pinfo.field.PF1DInMeanrate{i}(JJJ); pinfo.field.PF1DInPeakrate{i} = pinfo.field.PF1DInPeakrate{i}(JJJ);
    pinfo.field.PF1DBaseRate{i} = pinfo.field.PF1DBaseRate{i}(JJJ); 
    %%%% compute number of fields per traj
    fieldtraj = pinfo.field.PF1Devt{i}; evnamenow = unique(fieldtraj);
    pinfo.field.TrajName{i} = evnamenow; 
    pinfo.field.TrajNField{i} = NaN*ones(numel(evnamenow), 1); pinfo.field.TrajMeanRate{i} = NaN*ones(numel(evnamenow), 1);
    for (tt = 1:numel(evnamenow))
        III = find(strcmp(fieldtraj, evnamenow{tt}));  pinfo.field.TrajNField{i}(tt) = numel(III);
        III = find(strcmp(AllevName, evnamenow{tt}));  pinfo.field.TrajMeanRate{i}(tt) = Allevrate(III);
    end

    
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

