function [pinfo, data] = DataManager_CrrTransferCellVariables_Callback(pinfo, data) 
%transfer variables from the parent spikedb

%get the parent spikedb first
ok = 1; cpathname = fullfile(cd, '*.spikedb'); 
[fname, pname] = uigetfile(cpathname, 'Select a parent spike database of the current crrdb:');
if (fname ~= 0)
    spikedbfile = fullfile(pname, fname);
    S = load(spikedbfile, '-mat'); %load the file
    spikepinfo = S.pinfo; spikedata = S.data; S = []; %load the plot structure
end

%get category and variable names
if ok
    subf = fieldnames(spikepinfo); catnow = []; varnow = []; 
    [sel, ok] = listdlg('ListString', subf, 'PromptString', 'Select one or more categories', 'SelectionMode', 'multiple');
    if ok
       for (i = 1:numel(sel))
            catnow{i} = subf{sel(i)}; 
            allvar = fieldnames(spikepinfo.(catnow{i}));
            strok = ['Select ', catnow{i}, ' variables to transfer:'];
            [sss, ok] = listdlg('ListString', allvar, 'PromptString', strok, 'SelectionMode', 'multiple');
            if ok
               for (j = 1:numel(sss)) var{i}{j} = allvar{sss(j)}; end
            end
       end
    end
end
%%%%%transfer varaibles
if ok
   %%%%first need to identify cell indices
   nspike = numel(spikepinfo.general.clname); npair = numel(pinfo.general.clname); 
   spikeid = cell(nspike,1); pairind1 = zeros(npair, 1); pairind2 = zeros(npair, 1);
   for (i = 1:nspike)
       spikeid{i} = strcat(spikepinfo.general.finaldir{i}, '_', spikepinfo.general.clname{i});
   end
   for (i = 1:npair)
       pdir = pinfo.general.finaldir{i}; pcl = pinfo.general.clname{i};
       ss = strfind(pdir, '__'); tt = strfind(pcl, '__');
       pairid1 = strcat(pdir(1:ss-1), '_', pcl(1:tt-1));
       pairid2 = strcat(pdir(ss+2:numel(pdir)), '_', pcl(tt+2:numel(pcl)));
       pairind1(i) = find(strcmp(spikeid, pairid1)); pairind2(i) = find(strcmp(spikeid, pairid2));
   end
   %%%%transfer variable values
   for (i = 1:numel(catnow))
       newcat1 = strcat(catnow{i},'1'); newcat2 = strcat(catnow{i},'2'); 
       for (j = 1:numel(var{i}))
           newvar1 = strcat(var{i}{j},'1'); newvar2 = strcat(var{i}{j},'2'); 
           pinfo.(newcat1).(newvar1) = spikepinfo.(catnow{i}).(var{i}{j})(pairind1);
           pinfo.(newcat2).(newvar2) = spikepinfo.(catnow{i}).(var{i}{j})(pairind2);
       end
   end
end
if ~ok disp(['-----> action canceled']); end
disp('**********************');







