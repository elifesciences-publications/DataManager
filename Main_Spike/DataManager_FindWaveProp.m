function [pinfo,data] = DataManager_FindWaveProp(pinfo,data, cellind, vv)
%%find waveform properties for each spike: channel 1,2,3,4 amplitudes, maximum height, maximum width, all average values
%%all amplitude values in AD units (need gain to switch to mV), time values in data point (need sampling frequency to switch to second)
if (~isempty(cellind))
   nspike = numel(pinfo.general.animalname);  if (~isfield(pinfo, 'waveform')) pinfo.waveform = []; end
   if (~isfield(pinfo.waveform, 'champ')) pinfo.waveform.champ = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'maxwd')) pinfo.waveform.maxwd = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'maxht'))    pinfo.waveform.maxht = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'posamp'))    pinfo.waveform.posamp = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'negamp'))    pinfo.waveform.negamp = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'npampratio'))    pinfo.waveform.npampratio = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'hfwd'))    pinfo.waveform.hfwd = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'hlwd'))    pinfo.waveform.hlwd = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'hlht'))    pinfo.waveform.hlht = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'maxposamp'))    pinfo.waveform.maxposamp = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'maxnegamp'))    pinfo.waveform.maxnegamp = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'maxnpampratio'))    pinfo.waveform.maxnpampratio = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'maxhfwd'))    pinfo.waveform.maxhfwd = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'maxhlwd'))    pinfo.waveform.maxhlwd = cell(1,nspike); end
   if (~isfield(pinfo.waveform, 'maxhlht'))    pinfo.waveform.maxhlht = cell(1,nspike); end
end

for (jjjk = 1:numel(cellind))
    i = cellind(jjjk);
    disp(strcat('-----> wave prob ---', pinfo.general.parmfile{i}));
    gain = pinfo.general.gain{i}; fs = pinfo.parm.fs(i); spikebasepoint = pinfo.parm.spikebasepoint(i);
    if (~isempty(pinfo.general.wavefile{i}))
        wavefile = pinfo.general.wavefile{i};
        [spiketime, wv] = ReadSpikeWaveSPW(wavefile);
        %%%%%spiketime in second, wv = [32, 4, nspike]
        nspike = numel(spiketime);
        if (nspike >0)
           hfwd = zeros(4,1); hlwd = zeros(4,1); hlht = zeros(4,1); posamp = zeros(4,1); negamp = zeros(4,1); ampratio = zeros(4,1); 
           %calculate average waveform
           wvavg = sum(wv,3)/nspike;
           for (j = 1:4)
               [hfwd(j), hlwd(j), hlht(j), posamp(j), negamp(j), ampratio(j)] = findparam(wvavg(:,j)*gain(j), fs, spikebasepoint);
           end
           if (vv == 1) plotspikewaveavg(wavefile, wvavg, gain); end
           wvavg = []; wv = []; spiketime = [];
           pinfo.waveform.posamp{i} = posamp;
           pinfo.waveform.negamp{i} = negamp;
           pinfo.waveform.npampratio{i} = ampratio;
           pinfo.waveform.hfwd{i} = hfwd;
           pinfo.waveform.hlwd{i} = hlwd;
           pinfo.waveform.hlht{i} = hlht;
           [pinfo.waveform.maxposamp{i}, ikk] = max(posamp);
           pinfo.waveform.maxnegamp{i} = negamp(ikk);
           pinfo.waveform.maxnpampratio{i} = ampratio(ikk);
           pinfo.waveform.maxhfwd{i} = hfwd(ikk);
           pinfo.waveform.maxhlwd{i} = hlwd(ikk);
           pinfo.waveform.maxhlht{i} = hlht(ikk);
        end
    end
    if (~isempty(pinfo.general.parmfile{i}))
        %get spike parm file (pinfo.filename.parmfile{i})
        parmfile = pinfo.general.parmfile{i};
        S = load(parmfile, '-ascii'); %this loading method suppose not sensitive to different parameter columns
        if (~isempty(S))
           pinfo.waveform.champ{i} = [mean(S(:,2))*gain(1) mean(S(:,3))*gain(2) mean(S(:,4))*gain(3) mean(S(:,5))*gain(4)];
           pinfo.waveform.maxwd{i} = mean(S(:,6))/fs;
           pinfo.waveform.maxht{i} = mean(S(:,7))*max(gain);
           S = [];
        end
    end
end

function [hfwd, hlwd, hlht, posamp, negamp, ampratio] = findparam(wvavg, fs, spikebasepoint)
if (spikebasepoint >= 1)
    baseline = mean(wvavg(1:spikebasepoint)); %take the first xx points as baseline
else
    baseline = 0;
end
wvavg = wvavg - baseline;
[posamp, postime] = max(wvavg); [negamp, negtime] = min(wvavg); 
if (posamp ~= 0) 
    ampratio = abs(negamp/posamp); 
else
    ampratio = NaN;
end
hlht = posamp - negamp;
hlwd = (negtime - postime)/fs; 

%%%pos peak half width in s
hfwd = NaN; hftime1 = []; hftime2 = [];
if (negtime > postime)
    for (i = 1:postime)
        if (wvavg(i) >= posamp/2)
            hftime1 = i; break
        end 
    end
    for (i = postime:negtime)
        if (wvavg(i) <= posamp/2)
            hftime2 = i; break
        end
   end
   if ((~isempty(hftime1)) && (~isempty(hftime2)))
      hfwd = (hftime2 - hftime1)/fs;
   end
end

function plotspikewaveavg(wavefile, wvavg, gain)
[npoint,nch] = size(wvavg);
tlim = [-2 npoint+2]/32.556; maxV = max([max(gain(1)*wvavg(:,1)) max(gain(2)*wvavg(:,2)) max(gain(3)*wvavg(:,3)) max(gain(4)*wvavg(:,4))]);
minV = min([min(gain(1)*wvavg(:,1)) min(gain(2)*wvavg(:,2)) min(gain(3)*wvavg(:,3)) min(gain(4)*wvavg(:,4))]);
vlim = [minV-5 maxV+5];
hf = figure('Name', strcat(wavefile, '_AverageWaveforms'), 'NumberTitle', 'off');
ha1 = axes('Parent', hf, 'NextPlot', 'add', 'Position', [0.1 0.1 0.15 0.9],...
    'XLim', tlim, 'YLim', vlim); xlabel('ms'); ylabel('uV');
ha2 = axes('Parent', hf, 'NextPlot', 'add', 'Position', [0.35 0.1 0.15 0.9],...
    'XLim', tlim, 'YLim', vlim);
ha3 = axes('Parent', hf, 'NextPlot', 'add', 'Position', [0.6 0.1 0.15 0.9],...
    'XLim', tlim, 'YLim', vlim);
ha4 = axes('Parent', hf, 'NextPlot', 'add', 'Position', [0.85 0.1 0.15 0.9],...
    'XLim', tlim, 'YLim', vlim);
plot ([1:npoint]/32.556, gain(1)*wvavg(:,1), 'Parent', ha1, 'LineWidth', 2, 'Color', [1 0 0]);
xlabel('ms'); ylabel('uV');
plot ([1:npoint]/32.556, gain(2)*wvavg(:,2), 'Parent', ha2, 'LineWidth', 2, 'Color', [1 0 0]);
plot ([1:npoint]/32.556, gain(3)*wvavg(:,3), 'Parent', ha3, 'LineWidth', 2, 'Color', [1 0 0]);
plot ([1:npoint]/32.556, gain(4)*wvavg(:,4), 'Parent', ha4, 'LineWidth', 2, 'Color', [1 0 0]);


