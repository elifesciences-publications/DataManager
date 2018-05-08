function DataManager_GenerateCrrDatabase_Callback
%generate a new crr database from a spikedb

%get the spikedb first
hf = gcbf; fname = get(hf, 'Name');
mmind = strfind(fname, '__'); extfile = fname(mmind+2:numel(fname)); [~, ~, ee] = fileparts(extfile);
spikedbfile = []; ok = 1; 
hfig = []; %%%%the figure handle input for selecting only selected groups of cells in the crr computation
disp('Create new cross-/auto-correlation database'); 
if (strcmp(ee, '.spikedb'))
    spikedbfile = extfile;
else
    ok = 0; disp('---------> not a spikedb; cannot compute crr');
end
spikepinfo = []; spikedata = [];
if ok
    spikepinfo = getappdata(hf, 'pinfo'); spikedata = getappdata(hf, 'data');
    hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); 
    grpind = find(groupselection == 1); cellind = cell(1, numel(grpind)); 
    for (kk = 1:numel(grpind)) cellind{kk} = spikedata.grouplist.groupindex{grpind(kk)}; end
        %spikeselection(cellind) = ones(1, numel(cellind));
end

% %%%assign the selection index to 1 for newly added spikes, but keep the old spikes selection the same
% subf = fieldnames(spikepinfo.general); 
% allspikeselection = ones(1,numel(spikepinfo.general.(subf{1})));
% allspikeselection(1:ntotalspike) = spikeselection;    

if ok && (numel(spikepinfo.general.parmfile)>1)
   sel = {'Allgroups'; 'Withingroups'; 'Crossgroups'};
    sss = listdlg('ListString', sel, 'PromptString', 'Spike selection');
    if ~isempty(sss)
       selectiontype = sel{sss};
    else
       ok = 0;
    end
end
if ok 
    cpathname = fullfile(cd, '*.crrdb'); 
    [fname, pname] = uiputfile(cpathname, 'Select a crr database file to save:');
    if (fname ~= 0)
       writefilename = fullfile(pname, fname);
    else
        ok = 0;
    end
    %[pp,nn,ee] = fileparts(spikedbfile);
    %writefilename = fullfile(pp, strcat(nn, '.crrdb'));
   if ok
    disp('-----> Get all spike pairs');
    %%%gather all the general information
    pinfo = []; pinfo.work = struct([]); data = [];
    %[pinfo, data] = DataManager_FindCrrAllSpikePairs(pinfo, spikepinfo, spikedata, allspikeselection);
    [pinfo, data] = DataManager_FindCrrAllSpikePairs(pinfo, spikepinfo, spikedata, cellind, selectiontype);
    %get parameters for all the calculations
    disp('-----> Set parameters for computing crr');
    [pinfo, data] = DataManager_FindCrrParm(pinfo, data);
    %%%compute initial variables
    vv = 0; cellindnow = 1:numel(pinfo.general.clname); 
    disp('-----> Compute intitial variables of the Crr database');
    [pinfo,data] = DataManager_ComputeCrrInit_Callback(pinfo, data, cellindnow, vv);
    %%%%set group list
    if (~isempty(pinfo.general.parmfile))
          data.grouplist.groupname{1} = 'List0'; data.grouplist.groupindex{1} = 1:numel(pinfo.general.clname);
          data.grouplist.grouptype{1} = 'Manual'; data.grouplist.groupcrit{1} = [];
          
          SSS = whos('pinfo');
          TTT = whos('data');
          disp(['--------> later pinfo size: ', num2str(SSS.bytes)]);
          disp(['--------> later data size: ', num2str(TTT.bytes)]);
          
          %%%Save the result
          save(writefilename, 'pinfo', 'data');
          %%plot reult  
          if (isempty(extfile))
             hmain = hf;
          else
             hmain = DataManager_DataManager_Callback;
          end
          DataManager_PlotSpikeDatabase(hmain, pinfo, data, 'CellPair', 'clname', '.crrdb'); 
          set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
          setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); pinfo = []; data = [];
    else
          disp(['-------------> no cell pairs found!']);
    end
   end
else
    disp(['-------------> action cancelled or no more than 2 cells found!']);
end
disp('**********************');

function [eegev, nev] = geteegevent(dbfile, cat, var)
nev = 0; eegev = [];
S = load(dbfile, '-mat'); %load the file
[~,~,ee] = fileparts(dbfile);
if (~isempty(strfind(ee, 'db')))
    if (strcmp(ee, '.eegdb'))
       eeg = S.eeg; eegdata = S.eegdata;
    elseif (strcmp(ee, '.behavdb'))
       eeg = S.behav; eegdata = S.bhdata;
    else
       eeg = S.pinfo; eegdata = S.data; 
    end
    S = []; 
    for (i = 1:numel(var))
        if isfield(eeg, cat{i}) && isfield(eeg.(cat{i}), var{i})
           [eegevnow, nevnow] = getev(eeg, eegdata, cat{i}, var{i}, ee);
           nev = nev + nevnow; eegev = combineeegev(eegev, eegevnow);
        else
           disp(['-------------> warning: ', cat{i}, '.', var{i}, ' is not a valid variable!']);
        end
    end
end

function eegev = combineeegev(eegev, eegevnow)
if (~isempty(eegevnow))
    if (isempty(eegev))
        eegev = eegevnow;
    else
        nev = numel(eegevnow.clname);
        if (nev >0)
            subf = fieldnames(eegev);
            for (i = 1:numel(subf))
                [A,B] = size(eegev.(subf{i}));
                if (A ==1)
                    if (isfield(eegevnow, subf{i}))
                        eegev.(subf{i}) = [eegev.(subf{i}) eegevnow.(subf{i})];
                    else
                        ttt = cell(1, nev); eegev.(subf{i}) = [eegev.(subf{i}) ttt];
                    end
                elseif (B ==1)
                    if (isfield(eegevnow, subf{i}))
                        eegev.(subf{i}) = [eegev.(subf{i}); eegevnow.(subf{i})];
                    else
                        ttt = cell(nev,1); eegev.(subf{i}) = [eegev.(subf{i}); ttt];
                    end
                end
            end
        end
   end
end

function [ev, nev] = getev(eeg, eegdata, cat, var, ee)
nev = 0; ev = [];  %%%%%single field now: eg. car.var=ripple.sessStartT
neeg = numel(eeg.general.finaldir); 
%%%define checkstr: this is to deal with the facts that cell lists in some datadb only span a session: need to create a str shared by the sessions on the same day
checkstr = cell(1, neeg); idstr = [];
if (strcmp(ee, '.eegdb'))
    idstr = 'eegTT';
elseif (strcmp(ee, '.spikedb'))%%%spikedb already span whole day
    idstr = 'clname';     
end %%%for behavdb or others, all items belong to the same finaldir are grouped 
if (~isempty(idstr))
   for (i = 1:neeg) checkstr{i} = strcat(eeg.general.finaldir{i}, filesep, eeg.general.(idstr){i}); end
else
   for (i = 1:neeg) checkstr{i} = eeg.general.finaldir{i}; end
end
%%%combine different event files (during different sessions) together to make a spike file for the whole day
unistr = unique(checkstr);
for (i = 1:numel(unistr)) 
    iii = find(strcmp(checkstr, unistr{i})); 
    datanow = [];
    for (j = 1:numel(iii))
        %disp(['-----> var now: ', cat, '.', var]);
        %disp(['-----> var size: ', num2str(size(eeg.(cat).(var)))]);
        %disp(['-----> index now: ', num2str(iii(j))]);
         [A,B] = size(eeg.(cat).(var){iii(j)});
         if (A ==1)
             datanow = [datanow eeg.(cat).(var){iii(j)}]; 
         elseif (B == 1)
             datanow = [datanow; eeg.(cat).(var){iii(j)}]; 
         end
    end
    if (~isempty(datanow)) %%%if variable is empty, it will not be added to the events; This poses a problem in some situation: some variables may be computed but results are empty: eg. no ripples found
        nev = nev + 1; ev.data{nev} = datanow; ev.eventtimes{nev} = cell(0,1); 
        ev.sessionname{nev} = cell(0,1); ev.eventname{nev} = cell(0,1); 
        ev.sessionstartT{nev} = []; ev.sessionendT{nev} = []; ev.sessionlength{nev} = [];
        %%%%assign all the general variables required by .spikedb
        for (j = 1:numel(iii))
            ev.eventname{nev} = [ev.eventname{nev} eeg.general.eventname{iii(j)}];
            ev.eventtimes{nev} = [ev.eventtimes{nev} eegdata.event.eventtimes{iii(j)}];
        end
        if (strcmp(ee, '.eegdb')) || (strcmp(ee, '.behavdb'))
            for (j = 1:numel(iii))
            ev.sessionname{nev} = [ev.sessionname{nev}; eeg.general.sessname{iii(j)}];
            ev.sessionstartT{nev} = [ev.sessionstartT{nev}; eeg.general.sessstartT{iii(j)}];
            ev.sessionendT{nev} = [ev.sessionendT{nev}; eeg.general.sessendT{iii(j)}];
            ev.sessionlength{nev} = [ev.sessionlength{nev}; eeg.general.sesslength{iii(j)}];
            end
            
        else
            for (j = 1:numel(iii))
            ev.sessionname{nev} = [ev.sessionname{nev}; eeg.general.sessionname{iii(j)}];
            ev.sessionstartT{nev} = [ev.sessionstartT{nev}; eeg.general.sessionstartT{iii(j)}];
            ev.sessionendT{nev} = [ev.sessionendT{nev}; eeg.general.sessionendT{iii(j)}];
            ev.sessionlength{nev} = [ev.sessionlength{nev}; eeg.general.sessionlength{iii(j)}];
            end
        end
        [datename, tok] = strtok(eeg.general.datedir{iii(1)}, '_');
        if (strcmp(ee, '.eegdb')) || (strcmp(ee, '.behavdb'))
            ev.clname{nev} = strcat(datename, '_', eeg.general.eegTT{iii(1)}, '_', cat, '_', var); %{str};
        else
            if (isfield(eeg.general, 'clname'))
                ev.clname{nev} = strcat(datename, '_', eeg.general.clname{iii(1)}, '_', cat, '_', var); %{str};
            else
                ev.clname{nev} = strcat(datename, '_', cat, '_', var); %{str};
            end
        end
        if (~isempty(idstr))
            ev.TTname{nev} = eeg.general.(idstr){iii(1)};
        else
            ev.TTname{nev} = [];
        end
        ev.parmfile{nev} = unistr{i}; ev.wavefile{nev} = []; ev.gain{nev}=[]; 
        subf = fieldnames(eeg.general);
        for (j = 1:numel(subf))
            if (~(isfield(ev, subf{j}))) || (numel(ev.(subf{j})) < nev) %%%if field not defined yet as above
               ev.(subf{j}){nev} = eeg.general.(subf{j}){iii(1)};
            end
        end
    end
end

function [cat, var] = geteegvar(ff)
nvar = 0; cat = cell(1, nvar); var = cell(1, nvar);
ok = 1; ttt = ff{1};
if (~isempty(ttt))
    while ok
        [str, tok] = strtok(ttt, ';');
        if (~isempty(str))
            nvar = nvar + 1; [sss,ooo] = strtok(str, '.'); 
            cat{nvar} = sss; var{nvar} = ooo(2:numel(ooo));
        else
            ok = 0;
        end
        tok = tok(2:numel(tok)); %%%the leftovers in the input str
        if (~isempty(tok))
           kk = 1;
           for (i = 1:numel(tok))
            if (tok(i) ~= ' ')
               kk = i; break
            end
           end
           ttt = tok(kk:numel(tok)); %%%get rid of the preceding spaces
           if (isempty(ttt)) ok = 0; end
        else
            ok = 0;
        end
    end
end





