function DataManager_Save_Callback
hf = gcbf; fname = get(hf, 'Name');
mmind = strfind(fname, '__'); extfile = fname(mmind+2:numel(fname)); 
[pp, nn, ee] = fileparts(extfile);

if (strcmp(ee, '.behavdb'))
   behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
   save(extfile, 'behav', 'bhdata');
elseif (strcmp(ee, '.eegdb'))
   eeg = getappdata(hf, 'eeg'); eegdata = getappdata(hf, 'eegdata');
   save(extfile, 'eeg', 'eegdata');
else %(strcmp(ee, '.spikedb'))
   pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
   save(extfile, 'pinfo', 'data');
end

