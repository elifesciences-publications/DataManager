function DataManager_RettMice_OverlappingPlaceCellPairs
%%identifying cell pairs with overlapping place fields
hf = gcbf; dbtype = getappdata(hf, 'dbtype'); celllistvar = getappdata(hf, 'celllistvar'); celllisttitle = getappdata(hf, 'celllisttitle');
%%%%%%%added options for all types of databases
ok = 1;
if ~strcmp(dbtype, '.crrdb')
    msgbox('Not a crrdb. Quit'); ok = 0;
end
if ok
    hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
    if ~isfield(pinfo, 'field1') || ~isfield(pinfo, 'field2')
        msgbox('Field properties not imported. Quit'); ok = 0;
    end
end
if ok
    if ~isfield(pinfo.field1, 'PF1Devt1') || ~isfield(pinfo.field2, 'PF1Devt2')
        msgbox('Field trajectories (PF1Devt) not imported. Quit'); ok = 0;
    elseif ~isfield(pinfo.field1, 'PF1DBoundStart1') || ~isfield(pinfo.field2, 'PF1DBoundStart2')
        msgbox('Field start boundaries (PF1DBoundStart) not imported. Quit'); ok = 0;
    elseif ~isfield(pinfo.field1, 'PF1DBoundEnd1') || ~isfield(pinfo.field2, 'PF1DBoundEnd2')
        msgbox('Field end boundaries (PF1DBoundEnd) not imported. Quit'); ok = 0;    
    end
end
if ok
    inputIII=inputdlg({'Threshold for overlap (% of smaller field)'; 'New group name:'}, 'Pair Selection', 1, {'25'; 'FOverlap'}, 'on');
    if (~isempty(inputIII))
        [overlapthreshold, ok] = str2num(inputIII{1});
    else
        ok = 0;
    end
end
if ok
    fieldevt1 = pinfo.field1.PF1Devt1; fieldevt2 = pinfo.field2.PF1Devt2;
    Bstart1 = pinfo.field1.PF1DBoundStart1; Bstart2 = pinfo.field2.PF1DBoundStart2;
    Bend1 = pinfo.field1.PF1DBoundEnd1; Bend2 = pinfo.field2.PF1DBoundEnd2;
    %%%group selection
    hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection'); 
    grouplistpos= getappdata(hf, 'grouplistpos');
    %%%cell selection
    hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection');
    htext = getappdata(hspike, 'htext'); %every line (spike name) is a text object
    ntext = numel(htext); %number of text lines displayed: not all spikes being displayed
    cellind = []; grpind = find(groupselection == 1); 
    if numel(grpind) == 0
        msgbox('No group(s) selected. Quit'); ok = 0;
    end
end
if ok
    for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
    isoverlap = zeros(1, numel(cellind));
    for (k = 1:numel(cellind))
        i = cellind(k);
        isoverlap(k) = findifplacefieldoverlap(fieldevt1{i}, fieldevt2{i}, Bstart1{i}, Bstart2{i}, Bend1{i}, Bend2{i}, overlapthreshold);
    end
    indnow = cellind(find(isoverlap));
    %%%update groups
    ngroup = numel(data.grouplist.groupname);
    data.grouplist.groupname{ngroup+1} = inputIII{2}; data.grouplist.groupindex{ngroup+1} = indnow;
    ntype = ['Overlap fields > ' inputIII{1} '%('];
    ntype = strcat(ntype, data.grouplist.groupname{grpind(1)}, '_');
    for (ti = 2:numel(grpind))
           if (ti == numel(grpind))
               endstr = ')';
           else
               endstr = '_';
           end
           ntype = strcat(ntype, data.grouplist.groupname{grpind(ti)}, endstr); 
    end
    data.grouplist.grouptype{ngroup+1} = ntype; data.grouplist.groupcrit{ngroup+1} = [];
    %%%%%update the group plot
    delete(hgroup);
    hgroup = TextDisplayer(hf, grouplistpos, data.grouplist.groupname, 'GroupList', 'normalized');
    %save data
    setappdata(hf, 'hgroup', hgroup); setappdata(hf, 'data', data);
    %now re-set groupselectioin/spikeselection
    groupselection = getappdata(hgroup, 'selection'); spikeselection = 0*spikeselection; 
    groupsel = find(groupselection == 1); spikeselectindex = [];
    for (i = 1:numel(groupsel))
         spikeselectindex = union(spikeselectindex, data.grouplist.groupindex{groupsel(i)});
    end
    spikeselection(spikeselectindex) = ones(numel(spikeselectindex), 1);
    %update spike selection
    for (i = 1:ntext)
        tagtag = get(htext(i), 'Tag');
        [str, rem] = strtok(tagtag, '_');
        linenum = str2num(str); %current line number that selected
        if (spikeselection(linenum) == 0)
            set(htext(i), 'Color', [0 0 0]);
        else
            set(htext(i), 'Color', [1 0 0]);
        end
    end
    setappdata(hspike, 'selection', spikeselection); setappdata(hspike, 'htext', htext);
end

function isoverlap = findifplacefieldoverlap(fieldevt1, fieldevt2, Bstart1, Bstart2, Bend1, Bend2, thres)
isoverlap = 0; nfield1 = numel(Bstart1); nfield2 = numel(Bstart2);
for (i = 1:nfield1)
    S1 = Bstart1(i); E1 = Bend1(i);
    for (j = 1:nfield2)
        S2 = Bstart2(j); E2 = Bend2(j);
        if strcmp(fieldevt1{i}, fieldevt2{j})
           if (S2 < S1)
               if (E2 >= E1) 
                   isoverlap = 1; break; 
               elseif (E2<E1) && (E2>S1)
                   overlap = E2-S1;
                   if (overlap/(E1-S1) >= thres) || (overlap/(E2-S2) >= thres)
                       isoverlap = 1; break;
                   end
               end
           elseif (S2>=S1) && (S2<E1)
               if (E2<=E1)
                   isoverlap = 1; break;
               elseif (E2>E1)
                   overlap = E1-S2;
                   if (overlap/(E1-S1) >= thres) || (overlap/(E2-S2) >= thres)
                       isoverlap = 1; break;
                   end
               end
           end
        end
    end
end
    

