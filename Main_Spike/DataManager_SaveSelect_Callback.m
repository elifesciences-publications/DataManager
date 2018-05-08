function DataManager_SaveSelect_Callback
%%Save selected spikes into a new spike database
%%%%%%% Attention %%%%%%%%%%%%%
%%%%%Only cell arrays and non-empty numeric matrix assignments are allowed as entries such as:
%%%%%  pinfo.general.celltype{i}  or   pinfo.general.fieldnumber(i)

hf = gcbf; fname = get(hf, 'Name'); hgroup = getappdata(hf, 'hgroup');

groupselection = getappdata(hgroup, 'selection'); cellind = []; 
grpind = find(groupselection == 1); 

if (~isempty(grpind))
    mmind = strfind(fname, '__'); extfile = fname(mmind+2:numel(fname)); 
   [pp, nn, ee] = fileparts(extfile);
   if (strcmp(ee, '.behavdb'))
      pinfo = getappdata(hf, 'behav'); data = getappdata(hf, 'bhdata');
   elseif (strcmp(ee, '.eegdb'))
      pinfo = getappdata(hf, 'eeg'); data = getappdata(hf, 'eegdata');
   else %(strcmp(ee, '.spikedb'))
      pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
   end
   for (kk = 1:numel(grpind)) 
       if (~strncmpi(data.grouplist.groupname{grpind(kk)}, 'List0', 5))
          cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); 
       end
   end
   if (~isempty(cellind))
       disp('Save selected items into a new database');
       disp('-----> sorting data structure');
       pinfo = sortspike(pinfo, cellind);
       data = sortspike(data, cellind);
       disp('-----> writing data back to a file');
       [fname, pname] = uiputfile(fullfile(cd, '*.db'), 'Write the new spike database to:');
       writefilename = fullfile(pname, fname);
       if (strcmp(ee, '.behavdb'))
          behav = pinfo; bhdata = data; save(writefilename, 'behav', 'bhdata');
       elseif (strcmp(ee, '.eegdb'))
          eeg = pinfo; eegdata = data; save(writefilename, 'eeg', 'eegdata');
       else %(strcmp(ee, '.spikedb'))
          save(writefilename, 'pinfo', 'data');
       end
   else
       disp('-----> no items in the selected groups or only List0 selected. Nothing to save.');
   end
else
   disp('-----> no groups selected. Nothing to save.');
end

disp('**********************************');

function pinfo = sortspike(pinfo, selectindex)
%keep info if selected
%%%%%get all the fields in pinfo
fieldlist = fieldnames(pinfo);
nfield = numel(fieldlist);
for (iii = 1:nfield)
    if (strncmpi(fieldlist{iii}, 'grouplist', 6))
        pinfo.grouplist.groupname = {'List0'}; pinfo.grouplist.groupindex = {[1:numel(selectindex)]};
        pinfo.grouplist.grouptype = {'Manual'}; pinfo.grouplist.groupcrit = {[]};
    else
        subfieldlist = fieldnames(pinfo.(fieldlist{iii})); %%all the subfields
        for (jjj = 1:numel(subfieldlist))
            pinfo.(fieldlist{iii}).(subfieldlist{jjj}) = pinfo.(fieldlist{iii}).(subfieldlist{jjj})(selectindex);
        end
    end
end

