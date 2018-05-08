function DataManager_EEGShowEEGdata_Callback
%%Plot the computed eegdata
hf = gcbf; eeg = getappdata(hf, 'eeg'); eegdata = getappdata(hf, 'eegdata');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
hfield = getappdata(hf, 'hfield'); plotparm = getappdata(hf, 'plotparm');
hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); 
dateind = find(spikeselection == 1); %%%%%dateind here is actually eegind (eeg IDs)
plotparm = getappdata(hf, 'plotparm');
grpind = find(groupselection == 1); grpname = eegdata.grouplist.groupname(grpind);
hfield = getappdata(hf, 'hfield');

tagmark = get(gcbo, 'Tag'); 
if (strcmp(tagmark, 'psd'))
    plotPSD(eeg, eegdata, dateind, 'PSD', plotparm);
elseif (strcmp(tagmark, 'powerspec'))
    plotPowerSpec(eeg, eegdata, dateind, 'Power', plotparm);
elseif (strcmp(tagmark, 'signaldist'))
    plotsignaldist(eeg, eegdata, dateind, 'Raw signal', plotparm);
elseif (strcmp(tagmark, 'crosscorrelation'))
    plotcrosscorrelation(eeg, eegdata, dateind, 'crosscorrelation', plotparm);
elseif (strcmp(tagmark, 'triggeraverage'))
    plottriggeredaverages(hfield, eeg, eegdata, dateind, 'triggeredaverage', plotparm, grpname);
% elseif (strcmp(tagmark, 'evtraj'))
%     plotevVelDir(eeg, eegdata, dateind, 'evTraj');
% elseif (strcmp(tagmark, 'sessoccup'))
%     plotsesoccup(eeg, eegdata, dateind, 'sessOccupancy');
% elseif (strcmp(tagmark, 'evoccup'))
%     plotevoccup(eeg, eegdata, dateind, 'evOccupancy');
elseif (strcmp(tagmark, 'showeegdata')) | (strcmp(tagmark, 'showcluster0')) %%%for all other options to plot data
     showalldata(eeg, eegdata, dateind, 'eegdata');
end
disp('************************');

function plotcrosscorrelation(eeg, eegdata, dateind, tmark, plotparm)
p = {'Auto/cross?'; 'Bin size (s)'; 'Max lag (s)'}; ok = 1;
d = {'cross'; '1'; '50'};
II = inputdlg(p, 'Correlation parameters', 3, d, 'on'); %%%resizable window
if ~isempty(II)
      crrmode = II{1}; binsize = str2num(II{2}); maxlag = str2num(II{3}); 
else
      ok = 0;
end
if ok
    if strncmpi(crrmode, 'cross', 2) && (numel(dateind)<2)
        disp('----> no more than two traces are selected'); ok = 0;
    else
        ind = 1:min([2 numel(dateind)]);
        dateind = dateind(ind); %%%only do it for the first 2 selections
        if strncmpi(crrmode, 'auto', 2)
            for (j = 1:numel(dateind))
                i = dateind(j);
                if ~strcmp(eeg.parm.band{i}, 'cluster0')
                    fname = eeg.general.eegfile{i};
                    [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
                else
                    dat = eegdata.cluster0.cnt{i}; timestamp = eegdata.cluster0.timepoint{i};
                    fs = 1/eeg.parm.cl0Binsize(i); fname = eeg.general.eegfile{i};
                end
                if (plotparm.evselect == 0)
                    sT = eeg.general.sessstartT{i}; eT = eeg.general.sessendT{i}; 
                    plotcrrnow(dat, dat, timestamp, fs, binsize, maxlag, tmark, sT, eT, 'Whole session', fname, fname)
                else
                    eventnames = eeg.general.eventname{i};
                    for (k = 1:numel(eventnames))
                        evtTimes = eegdata.event.eventtimes{i}{k};
                        plotcrrnow(dat, dat, timestamp, fs, binsize, maxlag, tmark, evtTimes.start, evtimes.ent, eventnames{k}, fname, fname);
                    end
                end
            end
        else
            fdir1 = eeg.general.finaldir{dateind(1)}; fdir2 = eeg.general.finaldir{dateind(2)};
            sess1 = eeg.general.sessname{dateind(1)}; sess2 = eeg.general.sessname{dateind(2)};
            if strcmp(fdir1, fdir2) && strcmp(sess1, sess2)
                for (j = 1:numel(dateind))
                    i = dateind(j); fname{j} = eeg.general.eegfile{i};
                    if ~strcmp(eeg.parm.band{i}, 'cluster0')
                       [timestamp{j}, dat{j}, gain, fs(j)] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
                    else
                       dat{j} = eegdata.cluster0.cnt{i}; fs(j) = 1/eeg.parm.cl0Binsize(i);  timestamp{j} = eegdata.cluster0.timepoint{i};
                    end
                end
                if (numel(dat{1})==numel(dat{2})) && (fs(1)==fs(2))
                    if (plotparm.evselect == 0)
                        sT = eeg.general.sessstartT{i}; eT = eeg.general.sessendT{i}; 
                        plotcrrnow(dat{1}, dat{2}, timestamp{1}, fs(1), binsize, maxlag, tmark, sT, eT, 'Whole session', fname{1}, fname{2});
                    else
                        eventnames = eeg.general.eventname{i};
                        for (k = 1:numel(eventnames))
                            evtTimes = eegdata.event.eventtimes{i}{k};
                            plotcrrnow(dat{1}, dat{2}, timestamp{1}, fs(1), binsize, maxlag, tmark, evtTimes.start, evtTimes.ent, eventnames{k}, fname{1}, fname{2});
                        end
                    end
                else
                    disp('----> number of data points or sampling frequency does not match'); ok = 0; 
                end
            else
                 disp('----> first two traces not in the same session'); ok = 0;
            end
        end
    end
end

function plottriggeredaverages(hfield, eeg, eegdata, dateind, tmark, plotparm, grpname)
catname = fieldnames(eeg); ntimevar = 0; triggertime = []; triggername = []; triggerstr = []; ok = 1;
for (i = 1:numel(catname))
    subfield = fieldnames(eeg.(catname{i}));
    fieldselection = getappdata(hfield(i), 'selection');
    for (j = 1:numel(subfield))
        if (fieldselection(j) == 1) %%if a subfield is selected, delete it
            varnow = eeg.(catname{i}).(subfield{j});
            if (~strcmp(subfield{j}, 'eventname')) || (~strcmp(catname{i}, 'general'))
                for (k = 1:numel(varnow))
                    if isnumeric(varnow{k})
                        ntimevar = ntimevar + 1; triggername{ntimevar} = subfield{j}; 
                        triggertime{ntimevar} = varnow; triggerstr{ntimevar} = []; break;
                    end
                end
            else %%%if selecet an event file
                unievent = [];
                for (k = 1:numel(varnow))
                     unievent = unique(union(unievent, varnow{k}));
                end
                [sel,ok] = listdlg('ListString', unievent, 'SelectionMode', 'single', 'PromptString', 'Select events to average');
                if (ok)
                    ntimevar = ntimevar + 1; triggername{ntimevar} = subfield{j}; 
                    triggertime{ntimevar} = unievent(sel); disp(triggertime{ntimevar});
                    lststr = {'Start time'; 'End time'; 'Reference time'};
                    [sel,ok] = listdlg('ListString', lststr, 'SelectionMode', 'single', 'PromptString', 'Select which time to average');
                    if ok
                       triggerstr{ntimevar} = lststr{sel};
                    end
                end
            end
        end
    end
end
if ok
  if ntimevar == 0
    disp('----> no time variables selected'); ok = 0;
  else
    p = {'Bin size (s)'; 'Max lag (s)'};
    d = {'0.01'; '0.1'};
    II = inputdlg(p, 'Average parameters', 3, d, 'on'); %%%resizable window
    if ~isempty(II)
        binsize = str2num(II{1}); maxlag = str2num(II{2});
    else
        ok = 0;
    end
  end
end
if ok
    input = inputdlg({'Event keyword';'Event type'; 'Average over session?'}, 'Filtering', 3, {'sws'; 'sws'; 'no'}); 
    if (~isempty(input))
        evkeyword = input{1}; evkeytype = input{2}; 
        if strncmpi(input{3}, 'yes', 1)
            isavg = 1;
        else
            isavg = 0;
        end
    else
        ok = 0;
    end
end
if ok
   spikedbfile = [];
   button = questdlg('Average a spike database?');
   if (strcmp(button, 'Yes'))
      cpathname = fullfile(cd, '*.spikedb'); 
      [fname, pname] = uigetfile(cpathname, 'Select a spike database file to open:');
      if (fname ~= 0)
          spikedbfile = fullfile(pname, fname); 
          S = load(spikedbfile, '-mat'); pinfo = S.pinfo; data = S.data; S = []; %load the plot structure
          allgroupname = data.grouplist.groupname; 
          [sel,ok] = listdlg('ListString', allgroupname, 'SelectionMode', 'single', 'PromptString', 'Select a group to average');
          if (ok)
             dataindthere = data.grouplist.groupindex{sel};
          end
      else
          ok = 0;
      end
   elseif strcmp(button, 'Cancel') || (isempty(button))
        ok = 0;
   end
end
if ok
for (tt = 1:ntimevar)
    avg = []; xbin = [];
    for (j = 1:numel(dateind))
         i = dateind(j);
         fdir = eeg.general.finaldir{i}; sess = eeg.general.sessname{i}; 
         evName = eeg.general.eventname{i}; evType = eeg.parm.eventtype{i}; evT = eegdata.event.eventtimes{i}; 
         if strcmp(triggername{tt}, 'eventname') && (~isempty(triggerstr{tt})) %if select an event file to trigger
            triggertimenow = matchEventtimeout(evT, evType, evName, triggertime{tt}, triggerstr{tt});
            disp(triggertimenow(1:5)');
         else %%if select a time variable
            triggertimenow = matchtimeout(eeg, eegdata, fdir, sess, triggertime{tt});
            disp(triggertimenow(1:5));
         end
         triggertimenow = filtertriggertime(triggertimenow,evkeyword, evkeytype, evName, evType, evT); 
         if isempty(triggertimenow)
             disp(['-------> warning: trigger times not found in ', fdir, ' -> ', sess, ': ', triggername{tt}]);
         else
              if ~isempty(spikedbfile)
                 [avgnow, xbinnow] = findotherdbavgs(pinfo, data, dataindthere, fdir, sess, triggertimenow, binsize, maxlag, isavg);
              else
                 fname = eeg.general.eegfile{i};
                 if ~strcmp(eeg.parm.band{i}, 'cluster0')
                    [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
                 else
                    dat = eegdata.cluster0.cnt{i}; fs = 1/eeg.parm.cl0Binsize(i);  timestamp = eegdata.cluster0.timepoint{i};
                 end
                 [avgnow, xbinnow] = findtriggeredaverage(dat, timestamp, fs, binsize, maxlag, triggertimenow, isavg);
              end
              if ~isempty(avgnow)
                 if isempty(avg)
                     avg = [avg; avgnow]; xbin = xbinnow;
                 else
                     if (numel(xbinnow) == numel(xbin))
                         avg = [avg; avgnow]; 
                     else
                         disp(['-------> warning: trace does not match with the previous one: ', fdir, ' -> ', sess]);
                     end
                 end
              end
         end
    end
    [mm,nn] = size(avg);
    if (mm>1)
       se = std(avg)/sqrt(mm); avg = mean(avg);
       plottrigavg(avg, xbin, se, triggername{tt}, strcat(grpname, '_', evkeyword, '_', evkeytype));
    else
       disp('-------> nothing to average');
    end
end
end

function [avgout, xbinnow] = findotherdbavgs(pinfo, data, dataind, fdir, sessname, triggertimenow, binsize, maxlag, isavg)
avgout = []; xbinnow = [];
ncell = numel(dataind); sel = strcmp(pinfo.general.finaldir(dataind), fdir); 
for (i = 1:ncell)
    if isempty(find(strcmp(pinfo.general.sessionname{dataind(i)}, sessname))) sel(i) = 0; end
end
nsel = numel(find(sel==1));
if nsel > 0
   spiketime = data.spike.spiketime(dataind); spiketime = spiketime(sel==1);
   for (i = 1:numel(spiketime)) spiketime{i} = spiketime{i}'; end
   spiketime = cell2mat(spiketime);
   maxlagbin = round(maxlag/binsize); xbinnow = (-maxlagbin:maxlagbin)*binsize;
   ntime = numel(triggertimenow); nbin = numel(xbinnow); avgnow = zeros(ntime, nbin);
   for (i = 1:ntime)
       bintimes = xbinnow + triggertimenow(i);
       spiketimenow = spiketime( (spiketime>=min(bintimes)) & (spiketime<=max(bintimes)) );
       cnt = hist(spiketimenow, bintimes); %%%center bin times
       avgnow(i,:) = cnt/nsel/binsize;
   end
   avgout = avgnow;
   if isavg
      if ntime > 1
         avgout = mean(avgnow);
      end
   end
end

function plotPowerSpec(eeg, eegdata, dateind, tmark, plotparm)
neeg = numel(dateind);
for (k = 1:neeg)
    i = dateind(k);
    x = eeg.spec.freq{i}; sess = eeg.general.eegfile{i};
    if (plotparm.evselect == 0)
        y{1} = eeg.spec.sessPower{i}; tinfix{1} = 'Whole session';
        y{2} = eeg.spec.sessNormPower{i}; tinfix{2} = 'Whole session normalized';
    else
        y{1} = eeg.spec.evtPower{i}; tinfix{1} = 'Selected events';
        y{2} = eeg.spec.evtNormPower{i}; tinfix{2} = 'Selected events normalized';
        %%%below is the old way to plot event PSD
%         for (iik = 1:numel(eeg.general.eventname{i}))
%             if (~isempty(eegdata.spec.evtpower{i})) 
%                 y{iik} = eegdata.spec.evtpower{i}{iik}; 
%             else
%                 y{iik} = [];
%             end
%             tinfix{iik} = eeg.general.eventname{i}{iik};
%         end
    end
    for (iik = 1:numel(y))
        if (plotparm.setlog >0) 
            y{iik} = log10(y{iik}); stry = strcat('Log10(', tmark, ')');
        else
            stry = tmark;
        end
        %if (~isempty(y{iik}))
           hf = figure('Name', strcat(sess, '---', tmark)); hax = axes('Parent', hf, 'NextPlot', 'add'); str = tinfix{iik};
           line(x, y{iik}, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]); xlabel('Frequency (Hz)'); ylabel(stry);
           text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
        %end
    end
end

function plotsignaldist(eeg, eegdata, dateind, tmark, plotparm)
neeg = numel(dateind);
binvector = [plotparm.range(1):plotparm.bin:plotparm.range(2)]; 
for (k = 1:neeg)
    i = dateind(k); 
    if ~strcmp(eeg.parm.band{i}, 'cluster0')
       %%%read the eegfile data
       [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, eeg.parm.timeunit(i), eeg.parm.buffersize(i)); % timestamp in second, dat in mV 
       dat = dat-mean(dat); sess = eeg.general.eegfile{i}; 
    else
       dat = eegdata.cluster0.cnt{i}; sess = eeg.general.eegfile{i}; 
       timestamp = eegdata.cluster0.timepoint{i};
    end
    if (plotparm.evselect == 0)
        y{1} = dat; tinfix{1} = 'Whole session';
    else
        for (iik = 1:numel(eeg.general.eventname{i}))
            tinfix{iik} = eeg.general.eventname{i}{iik}; iii = [];
            startT = eegdata.event.eventtimes{i}{iik}.start; endT = eegdata.event.eventtimes{i}{iik}.ent;
            for (tt = 1:numel(startT))
                ii = find( (timestamp>=startT(tt)) & (timestamp<=endT(tt)) ); iii = union(iii, ii);
            end
            y{iik} = dat(iii);
        end
    end
    for (iik = 1:numel(y))
        if (plotparm.setlog == 1) | (plotparm.setlog == 4) 
            y{iik} = log10(abs(y{iik})); stry = strcat('Log10(', tmark, ')');
        else
            stry = tmark;
        end
        hf = figure('Name', strcat(sess, '---', tmark)); hax = axes('Parent', hf, 'NextPlot', 'add');
        Y = histc(y{iik}, binvector); bar(hax, binvector, Y, 'EdgeColor', [0 0 1], 'FaceColor', [0 0 1]); 
        xlabel(stry); ylabel('count');
        text('Interpreter', 'none', 'Parent', hax, 'String', tinfix{iik}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96]);
        strnow = strcat('nn=', num2str(numel(y{iik})), ';mean=', num2str(mean(y{iik})), ';std=', num2str(std(y{iik})), 'mV');
        text('Interpreter', 'none', 'Parent', hax, 'String', strnow, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92]);
    end
end

function plotPSD(eeg, eegdata, dateind, tmark, plotparm)
smooth = 0; Tsig = 0; NTsig = 0; Fsig = 0; NFsig = 0; ok = 1; 
showtimepoint = 100;  %default number of winodws to plot (1 window ~= 1 second as defined in parm.specWinShift)
SS = questdlg(['Smoothing the spectrogram?']);
if strcmp(SS, 'Yes')
   smooth = 1; 
   p = {'Time sigma (bin)'; 'Freq sigma (bin)'}; 
   d = {'1'; '1'};
   II = inputdlg(p, 'Smoothing parameters', 4, d, 'on'); %%%resizable window
   if ~isempty(II)
      Tsig = str2num(II{1}); NTsig = 5; Fsig = str2num(II{2}); NFsig = 5;
   else
      ok = 0;
   end
elseif strcmp(SS, 'Cancel')
   ok = 0;
end


if ok 
   p = {'Min freq (Hz)'; 'Max freq (Hz):'}; 
   d = {'1'; '400'};
   II = inputdlg(p, 'Frequency range', 2, d, 'on'); %%%resizable window
   if ~isempty(II)
      Fmin = str2num(II{1}); Fmax = str2num(II{2}); 
   else
      ok = 0;
   end
elseif strcmp(SS, 'Cancel')
   ok = 0;
end

if ok
if (~isempty(dateind))
    [MCroot, MCname, DAname, DEname] = CurrentVersion;
    hmain = figure('Name', DAname, 'NumberTitle', 'off', 'NextPlot', 'add', 'MenuBar', 'figure', 'Unit', 'normalized',...
       'OuterPosition', [0.204 0.425 0.779 0.576], 'Position', [0.207 0.429 0.773 0.496]);
    hf = uimenu(hmain, 'Label','Add Plot');
    uimenu(hf,'Label','Add New Spike','Callback','DataAnimator_AddSpike_Callback');
    uimenu(hf,'Label','Add New Position','Callback','DataAnimator_AddPosition_Callback');
    uimenu(hf,'Label','Add New EEG','Callback','DataAnimator_AddEEG_Callback');
    uimenu(hf,'Label','Add New Powergram/ratio','Callback','DataAnimator_AddPower_Callback', 'Separator', 'on');
    uimenu(hf,'Label','Add New PSD','Callback','DataAnimator_AddPSD_Callback');
    uimenu(hf,'Label','Add New SleepClass','Callback','DataAnimator_AddSleepClass_Callback');
    ndatafile = 0; displaysetting = []; haxes = {}; data = {};
    setappdata(hmain, 'ndatafile', ndatafile); setappdata(hmain, 'displaysetting', displaysetting);
    setappdata(hmain, 'haxes', haxes); setappdata(hmain, 'data', data);
    ndatafile = 0;
    ndisplay = 0; 
end    

for (j = 1:numel(dateind))
    k = dateind(j); filename = eeg.general.eegfile{k}; freqY = eeg.spec.freq{k}; timewinX = eeg.spec.timewin{k};
    iii = find( (freqY>=Fmin) & (freqY<=Fmax) ); freqY = freqY(iii);
    winpower = eegdata.spec.winpower{k}(iii,:); %psd[fy timewin]
    if (plotparm.setlog>0) winpower = log10(abs(winpower)); end%%%this change the plot into log-scale
    [mm,nn] = size(winpower);
    for (im = 1:mm)
        for (in = 1:nn)
            if (winpower(im,in) == -Inf) winpower(im,in) = NaN; end
        end
    end
 
    ntimepoint = numel(timewinX); nfreqpoint = numel(freqY);
    showtimepoint = min([showtimepoint ntimepoint]);
if (ntimepoint ~= 0)
    ndatafile = ndatafile + 1;
    displaysetting{ndatafile}{1} = ndatafile;
    displaysetting{ndatafile}{2} = filename;
    displaysetting{ndatafile}{3} = 'd'; % data type
    displaysetting{ndatafile}{4} = [0 0]; %default plot, continuous mode
    displaysetting{ndatafile}{5} = [ntimepoint nfreqpoint 0 ndisplay showtimepoint 2 nfreqpoint+1]; % 
    displaysetting{ndatafile}{6} = [0.5 0.5 0.5]; %default backgrnd
    displaysetting{ndatafile}{7} = [0 0 1]; %default line color
    displaysetting{ndatafile}{8} = 0.5; %default dot line width
    displaysetting{ndatafile}{9} = 'none'; %default marker
    set(0, 'ShowHiddenHandles', 'on');
    H = get(hmain, 'Children');
    TypeH = get(H, 'Type');
    for (i = 1:size(TypeH))
        if (strcmp(TypeH{i}, 'uicontrol')) delete(H(i)); end
    end
    delete(findobj(hmain, 'Tag', 'colorbar'));
    set(0, 'ShowHiddenHandles', 'off');
    
    setappdata(hmain, 'ndatafile', ndatafile);
    setappdata(hmain, 'displaysetting', displaysetting);  

    %%re-allocacte spaces for all the plots
    [axposvector, uiposvector] = DataAnimator_PlotSpaceAuto(ndatafile,0.8);
    %%re-plot all axes and uicontrols
    DataAnimator_Replot(haxes, ndatafile-1, axposvector, uiposvector);
    
    %%add new pSD plot
    displaysetting = getappdata(hmain, 'displaysetting');
    data = getappdata(hmain, 'data');

%%read settings
idnum = ndatafile;
fname = displaysetting{idnum}{2};
ntimepoint = displaysetting{idnum}{5}(1); %total number of time points
nfreqpoint = displaysetting{idnum}{5}(2); %total number of frequency points
finitpos = displaysetting{idnum}{5}(3); %binary data start postion after the header
ndisplay = displaysetting{idnum}{5}(4); %start display number
showtimepoint = displaysetting{idnum}{5}(5); %display time period in time points
startfreq = displaysetting{idnum}{5}(6); %starting frequency index
endfreq = displaysetting{idnum}{5}(7); %end frequency index
maxdisplay = ceil(ntimepoint/showtimepoint);

data{idnum}(1,1) = 0; %first element of the data matrix
data{idnum}(1,2:ntimepoint+1) = timewinX; %read time stamps to first row
data{idnum}(2:nfreqpoint+1,1) = freqY; %read frequency points to first column
data{idnum}(2:nfreqpoint+1, 2:ntimepoint+1) = winpower; 

minfreq = min(data{idnum}(2:nfreqpoint+1,1)); 
maxfreq = max(data{idnum}(2:nfreqpoint+1,1)); 

freqindex = find((data{idnum}(:,1)>=minfreq)&(data{idnum}(:,1)<=maxfreq));
startfreq = max([min(freqindex) 2]); 
endfreq = max(freqindex);
displaysetting{idnum}{5}(6) = startfreq;
displaysetting{idnum}{5}(7) = endfreq;

A = data{idnum}(startfreq:endfreq, 2:ntimepoint+1); 
% %%%%smoothing the power spectrum
if smooth == 1
   disp('----------> smoothing the spectrogram ......');  
   A = TwoDSmooth_separate(A, Tsig, NTsig, Fsig, NFsig);
   %B = TwoDSmooth_fast(A, Tsig, Fsig, NTsig, NFsig);
   %occu = ones(size(A)); B = TwoDSmooth_new(A, occu, Tsig, Fsig, NTsig, NFsig); 
   data{idnum}(startfreq:endfreq, 2:ntimepoint+1) = A; 
end
minvalue = 0; maxvalue = 1;
if ~isempty(A)
   [mm,nn] = size(A); BB = reshape(A, mm*nn, 1); BB = BB(~isnan(BB));
   ma = mean(BB); ss = std(BB);
   minvalue = ma-5*ss; maxvalue =ma+5*ss; 
end
%minvalue = min(min(data{idnum}(startfreq:endfreq, 2:ntimepoint+1)));
%maxvalue = max(max(data{idnum}(startfreq:endfreq, 2:ntimepoint+1)));
if (isempty(minvalue)) | (isempty(maxvalue)) | (minvalue == maxvalue) | isnan(minvalue) | isnan(maxvalue)
    minvalue = 0; maxvalue = 1; disp('----------> psd data seem not meaningful!'); 
end
    
posvec = axposvector{idnum};
haxes{idnum} = axes('Parent', hmain, 'Units', 'normalized', 'Position', posvec, 'FontSize', 8, 'Visible', 'off', 'Color', [1 1 1],...
    'CLim', [minvalue maxvalue], 'CLimMode', 'manual', 'XTickMode', 'manual', 'XTickLabelMode', 'manual');
displaysetting{idnum}{10} = [minvalue maxvalue];
%displaysetting{idnum}{10} = [-3 6];
ndatapoint=size(data{idnum},1);
   disp('----------> display the spectrogram ......'); 
   starttimepoint = ndisplay*showtimepoint+2; % PSD time starts from second column
   endtimepoint = min(starttimepoint+showtimepoint-1, ntimepoint);
   hp = pcolor(haxes{idnum}, data{idnum}(1,starttimepoint:endtimepoint), data{idnum}(startfreq:endfreq,1), data{idnum}(startfreq:endfreq, starttimepoint:endtimepoint));
   xlabel('Time (s)'); ylabel('Frequency (Hz)');
   set(hp, 'EdgeColor', 'none', 'LineStyle', 'none');
   colormap(jet);  
   hc = colorbar('vert', 'peer', haxes{idnum}, 'CLim', [minvalue maxvalue], 'CLimMode', 'manual');
   set(hc, 'Tag', 'colorbar');
   
   tickdis = (data{idnum}(1,endtimepoint)-data{idnum}(1,starttimepoint))/5; %10 ticks
   labelnumber = data{idnum}(1,starttimepoint):tickdis:data{idnum}(1,endtimepoint); 
   for (i = 1:numel(labelnumber))
       labelstr{i} = num2str(labelnumber(i), 8);
   end
   set(haxes{idnum}, 'XTickMode', 'manual', 'XTick', labelnumber, 'XTickLabelMode', 'manual', 'XTickLabel', labelstr);
   set(gca, 'FontSize', 8);
   
   titletext = strcat(num2str(idnum), ' = ', fname,' -PSD');
   htitle=text('String', titletext, 'Interpreter', 'none', 'Parent', haxes{idnum}, 'Tag', 'title',...
              'Position', [0.01 0.95], 'Units', 'normalized', 'Color', [1 0 0], 'FontSize', 8); % add a title
   %displaysetting{idnum}{10} = [0 50]; %
   %displaysetting{idnum}{10} = get(haxes, 'CLim');
   set(haxes{idnum}, 'CLim', displaysetting{idnum}{10});
   set(haxes{idnum}, 'CLimMode', 'manual');

%%plot a control bar with the axes (%%%main figure (hmain))
barpos = [posvec(1) posvec(2)-0.20*posvec(4)/0.75 posvec(3) posvec(4)*0.06];
hbar = uicontrol('Style', 'slider', 'Parent', hmain, 'Units', 'normalized', 'Position', barpos,...
                 'Min', 0, 'Max', maxdisplay, 'SliderStep', [1/maxdisplay, 10/maxdisplay], 'Callback', 'DataAnimator_PSDBar_Callback',...
                 'Tag', strcat('slider_',num2str(idnum)));
set(hbar, 'Value', ndisplay); %%initial bar position
setappdata(hmain, 'displaysetting', displaysetting);
setappdata(hmain, 'data', data);
setappdata(haxes{idnum}, 'hcolorbar', hc); setappdata(haxes{idnum}, 'hpcolor', hp);

    %%add new control panel
    DataAnimator_PSDControl(hmain, haxes{ndatafile}, ndatafile, uiposvector{ndatafile});
    setappdata(hmain, 'haxes', haxes);
end
end
end

function showalldata(eeg, eegdata, dateind, tmark)
if (~isempty(dateind))
    [MCroot, MCname, DAname, DEname] = CurrentVersion;
    hmain = figure('Name', DAname, 'NumberTitle', 'off', 'NextPlot', 'add', 'MenuBar', 'figure', 'Unit', 'normalized',...
       'OuterPosition', [0.204 0.425 0.779 0.576], 'Position', [0.207 0.429 0.773 0.496]);
    hf = uimenu(hmain, 'Label','Add Plot');
    uimenu(hf,'Label','Add New Spike','Callback','DataAnimator_AddSpike_Callback');
    uimenu(hf,'Label','Add New Position','Callback','DataAnimator_AddPosition_Callback');
    uimenu(hf,'Label','Add New EEG','Callback','DataAnimator_AddEEG_Callback');
    uimenu(hf,'Label','Add New Powergram/ratio','Callback','DataAnimator_AddPower_Callback', 'Separator', 'on');
    uimenu(hf,'Label','Add New PSD','Callback','DataAnimator_AddPSD_Callback');
    uimenu(hf,'Label','Add New SleepClass','Callback','DataAnimator_AddSleepClass_Callback');
    ndatafile = 0; displaysetting = []; haxes = {}; data = {};
    setappdata(hmain, 'ndatafile', ndatafile); setappdata(hmain, 'displaysetting', displaysetting);
    setappdata(hmain, 'haxes', haxes); setappdata(hmain, 'data', data);
    nread = 0; nreadbuffer = 4;  %default number of buffers being read and plot (4 buffs ~= 1 second data)
    defaultbuffsize = 512; %%%this is the default buffersize -EEG data will use actual buffersize
    defaultfreq = 2034; %%%this is the default sampling frequency - EEG data will use actual vale
end
for (j = 1:numel(dateind))
    k = dateind(j);
    if ~strcmp(eeg.parm.band{k}, 'cluster0') %%%if EEG traces
      filename = eeg.general.eegfile{k}; if exist(filename, 'file')~=2 filename = strrep(filename, 'I:', 'J:'); end
      displaysetting{j}{2} = filename; 
      gain = eeg.general.eeggain{k}; %%gain already in mV
      freq = eeg.general.freq{k}; timeunit = eeg.parm.timeunit(k); buffersize = eeg.parm.buffersize(k);
      fid = fopen(filename);  %first open a file for read
      while 1               %find where header ends
          tline = fgets(fid);
          if (strncmpi(tline, '-ADMaxValue', 10)) 
             [str, ent] = strtok(tline, ' '); ent = ent(2:numel(ent));  maxADvalue = str2num(ent); %maximum AD value in the data
          end
          if (strncmpi(tline, '%%ENDHEADER', 8)), break, end;
      end
      finitpos = ftell(fid);       %header end position, initial data reading position
      ok = fclose(fid);
      %get total number of buffers in the file
      maxread = EEG_ReadTotalBuffer(filename, nreadbuffer, finitpos, buffersize);
      %%%%assign parameters to display setting file and variables
      displaysetting{j}{1} = j;  % id number
      displaysetting{j}{3} = 'e'; % data type
      displaysetting{j}{4} = [0 0]; %default plot, continuous mode
      displaysetting{j}{5} = [nread nreadbuffer finitpos buffersize gain maxread freq maxADvalue];
      displaysetting{j}{6} = [0.5 0.5 0.5]; %default backgrnd
      displaysetting{j}{7} = [0 0 1]; %default line color
      displaysetting{j}{8} = 0.5; %default dot line width
      displaysetting{j}{9} = 'none'; %default marker
    
    else %%%if cluster0 power data
      filename = eeg.general.eegfile{k};
      displaysetting{j}{2} = filename;
      freq = NaN; finitpos = NaN;
      wstime = eeg.parm.cl0Binsize(k); %%%window size time
      shifttime = wstime; %%window shift time
      ntimestamp = numel(eegdata.cluster0.timepoint{k});
      %%%%assign parameters to display setting file and variables
      ndisplay = 0; %default start display: first display
      showtimepoint = (nreadbuffer*defaultbuffsize/defaultfreq)/wstime; %computed to align with EEG display time axes
      displaysetting{j}{1} = j;  % id number
      displaysetting{j}{3} = 'w'; % data type
      displaysetting{j}{4} = [0 0]; %default plot, continuous mode
      displaysetting{j}{5} = [wstime shifttime finitpos ndisplay showtimepoint ntimestamp]; % 
      displaysetting{j}{6} = [0.5 0.5 0.5]; %default backgrnd
      displaysetting{j}{7} = [0 0 1]; %default line color
      displaysetting{j}{8} = 0.5; %default dot line width
      displaysetting{j}{9} = 'none'; %default marker
    end     
        
      %%%%%get the band-specific data
      startT{j} = []; endT{j} = []; peakT{j} = []; startA{j} = []; endA{j} = []; peakA{j} = [];
      switch eeg.parm.band{k}
      case 'ripple'
           if isfield(eeg, 'ripple')
               if isfield(eeg.ripple, 'sessStartT')
                  startT{j} = eeg.ripple.sessStartT{k}; endT{j} = eeg.ripple.sessEndT{k}; peakT{j} = eeg.ripple.sessPeakT{k};
                  nrip = numel(eeg.ripple.sessStartT{k}); peakA{j} = eeg.ripple.sessAmp{k};
               elseif isfield(eeg.ripple, 'StartT')
                  startT{j} = eeg.ripple.StartT{k}; endT{j} = eeg.ripple.EndT{k}; peakT{j} = eeg.ripple.PeakT{k};
                  nrip = numel(eeg.ripple.StartT{k}); peakA{j} = eeg.ripple.Amp{k};
               end
               startA{j} = -1*eeg.ripple.startThreshold{k}*ones(1, nrip); endA{j} = -1*eeg.ripple.startThreshold{k}*ones(1, nrip); 
               %startA{j} = eeg.ripple.startThreshold{k}*ones(1, nrip); endA{j} = eeg.ripple.startThreshold{k}*ones(1, nrip); 
           end
      case 'theta'
           if (~isempty(strfind(filename, 'smooth'))) && (isfield(eeg, 'theta'))
              startT{j} = eeg.theta.maxTime{k}; endT{j} = eeg.theta.minTime{k}; 
              startA{j} = eeg.theta.maxAmp{k}; endA{j} = eeg.theta.minAmp{k}; 
           end
      case 'spindle'
           if (isfield(eeg, 'spindle'))
              startT{j} = eeg.spindle.maxTime{k}; endT{j} = eeg.spindle.minTime{k}; 
              startA{j} = eeg.spindle.maxAmp{k}; endA{j} = eeg.spindle.minAmp{k}; 
           end
      case 'hvs'
           if (isfield(eeg, 'hvs'))
              startT{j} = eeg.hvs.maxTime{k}; endT{j} = eeg.hvs.minTime{k}; 
              startA{j} = eeg.hvs.maxAmp{k}; endA{j} = eeg.hvs.minAmp{k}; 
           end
      case 'cluster0'
           if (isfield(eeg, 'cluster0'))
              ncl0 = numel(eeg.cluster0.sessStartT{k}); 
              startT{j} = eeg.cluster0.sessStartT{k}; endT{j} = eeg.cluster0.sessEndT{k}; peakT{j} = eeg.cluster0.sessPeakT{k}; 
              startA{j} = eeg.cluster0.startThreshold{k}*ones(1, ncl0); endA{j} = eeg.cluster0.startThreshold{k}*ones(1, ncl0); 
              peakA{j} = eeg.cluster0.sessAmp{k};
           end
      end
end

set(0, 'ShowHiddenHandles', 'on');
H = get(hmain, 'Children');
TypeH = get(H, 'Type');
for (i = 1:size(TypeH))
     if (strcmp(TypeH{i}, 'uicontrol')) delete(H(i)); end
end
delete(findobj(hmain, 'Tag', 'colorbar'));
set(0, 'ShowHiddenHandles', 'off');

ndatafile = numel(dateind);
setappdata(hmain, 'ndatafile', ndatafile);
setappdata(hmain, 'displaysetting', displaysetting);  
%%re-allocacte spaces for all the plots
[axposvector, uiposvector] = DataAnimator_PlotSpaceAuto(ndatafile,0.8);

for (j = 1:numel(dateind))
    k = dateind(j);
    if ~strcmp(eeg.parm.band{k}, 'cluster0') %%%if EEG traces add new EEG plot/control
       [haxes{j}, hbar] = DataAnimator_EEGDisplay(j, hmain, axposvector{j});
       DataAnimator_EEGControl(hmain, haxes{j}, hbar, j, uiposvector{j});
    else %%%if EEG traces add new EEG plot/control
        size(eegdata.cluster0.timepoint{k})
        size(eegdata.cluster0.cnt{k})
       haxes{j} = DataAnimator_PowergramDisplay(j, hmain, axposvector{j},...
           eegdata.cluster0.timepoint{k}', eegdata.cluster0.cnt{k}');
       DataAnimator_PowergramControl(hmain, haxes{j}, j, uiposvector{j}); 
    end
    setappdata(haxes{j}, 'startT', startT{j}); setappdata(haxes{j}, 'endT', endT{j}); setappdata(haxes{j}, 'peakT', peakT{j}); 
    setappdata(haxes{j}, 'startA', startA{j}); setappdata(haxes{j}, 'endA', endA{j}); setappdata(haxes{j}, 'peakA', peakA{j});
    %%%%%%%%%%plot events
    TT = get(haxes{j}, 'XLim'); minT = TT(1); maxT = TT(2);  VV = get(haxes{j}, 'YLim'); minV = VV(1); maxV = VV(2);
    indstart = find( (startT{j}>=minT) & (startT{j}<=maxT) );
    indend = find( (endT{j}>=minT) & (endT{j}<=maxT) );
    indpeak = find( (peakT{j}>=minT) & (peakT{j}<=maxT) );
    if ~strcmp(eeg.parm.band{k}, 'cluster0') %%%if EEG traces add new EEG plot/control
        for (k = 1:numel(indstart))
            line(startT{j}(indstart(k)), startA{j}(indstart(k)), 'Parent', haxes{j}, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20,...
                  'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]);
        end
        for (k = 1:numel(indend))
            line(endT{j}(indend(k)), endA{j}(indend(k)), 'Parent', haxes{j}, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20,...
                  'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
        end    
    else
        for (k = 1:numel(indstart))
            line([startT{j}(indstart(k)) startT{j}(indstart(k))], [0.8*minV 0.8*maxV], 'Parent', haxes{j}, 'LineWidth', 1, 'Color', [0 1 0]);
        end
        for (k = 1:numel(indend))
            line([endT{j}(indend(k)) endT{j}(indend(k))], [0.8*minV 0.8*maxV], 'Parent', haxes{j}, 'LineWidth', 1, 'Color', [1 0 0]);
        end  
    end
    for (k = 1:numel(indpeak))
         line(peakT{j}(indpeak(k)), peakA{j}(indpeak(k)), 'Parent', haxes{j}, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30,...
                  'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
    end
end
setappdata(hmain, 'haxes', haxes);

function plotcrrnow(dat1, dat2, timestamp, fs, binsize, maxlag, tmark, st, et, eventname, f1name, f2name)
binnow = round(binsize*fs); %%%%need to re-sample the data, otherwise the computation is too slow
ind = 1:binnow:numel(dat1); dat1 = dat1(ind); dat2 = dat2(ind); timestamp = timestamp(ind);
maxlagbin = round(maxlag/binsize); lags = (-maxlagbin:maxlagbin)*binsize;
STind = ones(size(st)); ETind = ones(size(et));
for (i = 1:numel(st))
    iii = find( (timestamp>=st(i)) & (timestamp<=et(i)) );
    if ~isempty(iii)
       STind(i) = min(iii); ETind(i) = max(iii);
    end
end
ccc = Utilities_FindXCorrCoef_self(dat1', dat2', maxlagbin, STind, ETind);
hg = figure('Name', tmark); 
hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
xlabel ('time lag (s)'); ylabel('Correlation'); 
line(lags, ccc, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]);
text('Interpreter', 'none', 'Parent', hax, 'String', f1name, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.94]);
text('Interpreter', 'none', 'Parent', hax, 'String', f2name, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.91]);
text('Interpreter', 'none', 'Parent', hax, 'String', eventname, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.88]);

function plottrigavg(avg, xbin, se, tmark, grpname)
hg = figure('Name', tmark); 
hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
xlabel ('time lag (s)'); ylabel('Triggered average'); 
line(xbin, avg, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]);
line(xbin, avg-se, 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
line(xbin, avg+se, 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
for (i = 1:numel(grpname))
     text('Parent', hax, 'Interpreter', 'none', 'String', grpname{i}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.96-(i-1)*0.04]);
end

function triggertimenow = matchEventtimeout(evT, evType, evName, eventname, timepoint)
triggertimenow = []; 
for (i = 1:numel(evName))
    if strcmp(eventname, evName{i})
        if strcmp(timepoint, 'Start time')
           triggertimenow = [triggertimenow; evT{i}.start];
        elseif strcmp(timepoint, 'End time')
           triggertimenow = [triggertimenow; evT(i).ent];
        elseif strcmp(timepoint, 'Reference time')
           triggertimenow = [triggertimenow; evT(i).ref];
        end
    end
end
function triggertimenow = matchtimeout(eeg, eegdata, fdir, sess, triggertime)
triggertimenow = [];
iii = find(strcmp(eeg.general.finaldir, fdir) & strcmp(eeg.general.sessname, sess) ); %%%this is the session index
nm = numel(iii);
if nm>0
    triggertime = triggertime(iii); sel = zeros(1, nm);
    for (i = 1:nm)
        if (isnumeric(triggertime{i}) && (~isempty(triggertime{i}))) sel(i) = 1; end
    end
    ij = find(sel==1); %disp(iii(ij))
    if numel(ij) == 1
        triggertimenow = triggertime{ij};
    else
        disp(['-------> warning: multiple matches in ', fdir, ' -> ', sess]);
    end
end
function triggertime = filtertriggertime(triggertime,evkeyword, evkeytype, evName, evType, evT) 
evsel = ones(1,numel(evName));
if ~isempty(evkeyword)
   for (i = 1:numel(evName))
       if isempty(strfind(lower(evName{i}), lower(evkeyword))) evsel(i) = 0; end 
   end
end
if ~isempty(evkeytype)
   for (i = 1:numel(evName))
       if isempty(strfind(lower(evType{i}), lower(evkeytype))) evsel(i) = 0; end 
   end
end
evpos = find(evsel == 1); evT = evT(evpos);
startT = []; entT = [];
for (i = 1:numel(evT))
    startT = [startT evT{i}.start]; entT = [entT evT{i}.ent];
end
if ~isempty(startT)
    iii = [];
    for (i = 1:numel(startT))
        iiok = find( (triggertime>=startT(i)) & (triggertime<=entT(i)) );
        iii = union(iii, iiok);
    end
    triggertime = triggertime(iii);
end

function [avgnow, xbinnow] = findtriggeredaverage(dat, timestamp, fs, binsize, maxlag, triggertimenow, isavg)
binnow = round(binsize*fs); %%%%need to re-sample the data, otherwise the computation is too slow
ind = 1:binnow:numel(dat); dat = dat(ind); timestamp = timestamp(ind);
maxlagbin = round(maxlag/binsize); xbinnow = (-maxlagbin:maxlagbin)*binsize; xbin = -maxlagbin:maxlagbin;
iii = find( (triggertimenow>timestamp(1)+maxlag) & (triggertimenow<timestamp(numel(timestamp))-maxlag) );
triggertimenow = triggertimenow(iii);
ntime = numel(triggertimenow); nbin = numel(xbinnow); avgnow = zeros(ntime, nbin);
for (i = 1:ntime)
    [~,iinow] = min(abs(timestamp-triggertimenow(i)));
    iii = xbin + iinow;
    avgnow(i,:) = dat(iii);
end
if isavg
   if ntime > 1
      avgnow = mean(avgnow);
   end
end
% for (i = 1:ntime)
%     iii = find( (timestamp>=triggertimenow(i)-maxlag) & (timestamp<=triggertimenow(i)+maxlag) );
%     numel(iii)
%     if (numel(iii)==nbin+1)
%         iii = iii(1:nbin);
%     elseif (numel(iii)==nbin-1)
%         iii = [iii; max(iii)];
%     end
%     size(avgnow)
%     size(dat(iii))
%     avgnow(i,:) = dat(iii);
% end





