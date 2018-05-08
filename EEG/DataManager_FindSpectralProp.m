function [eeg, eegdata] = DataManager_FindSpectralProp(eeg, eegdata, eegind, vv)
%%%do spectral analysis on broadband EEG files: (1) overall spectrum at different behavioral stages
%%%                                             (2) window-by-window spectrum
%%% require the following parameters
%%%    eeg.parm.specMinFreq/specMaxFreq; eeg.parm.specWinSize/specWinShift; eeg.parm.timeunit/buffersize;

%%% variables to assign
neeg = numel(eegind); ok = 1; alleeg = numel(eeg.general.finaldir);
if (neeg > 10)
    SS = questdlg(['Too many (>10) files to compute, Continue?']);
    if (~strcmp(SS, 'Yes')) ok = 0; end
end

if ok
   if (~isfield(eeg, 'spec')) eeg.spec = []; end
   if (~isfield(eeg.spec, 'freq')) eeg.spec.freq = cell(1, alleeg); end
   if (~isfield(eeg.spec, 'sessPower')) eeg.spec.sessPower = cell(1, alleeg); end
   if (~isfield(eeg.spec, 'sessNormPower')) eeg.spec.sessNormPower = cell(1, alleeg); end
   if (~isfield(eeg.spec, 'evtPower')) eeg.spec.evtPower = cell(1, neeg); end
   if (~isfield(eeg.spec, 'evtNormPower')) eeg.spec.evtNormPower = cell(1, neeg); end
   if (~isfield(eegdata, 'spec')) eegdata.spec = []; end
   if (~isfield(eegdata.spec, 'timewin')) eegdata.spec.timewin = cell(1, alleeg); end
   if (~isfield(eegdata.spec, 'winpower')) eegdata.spec.winpower = cell(1, alleeg); end
   if (~isfield(eegdata.spec, 'evtpower')) eegdata.spec.evtpower = cell(1, alleeg); end
   
   for (iiik = 1:neeg)
    i = eegind(iiik);
    %if strcmp(eeg.parm.band{i}, 'broad') |  %%only do this for broadband EEG traces - no no
    disp(['--------> spectral analysis: ', eeg.general.eegfile{i}]);
    %%%%get the EEg data
    if (~isempty(strfind(eeg.general.eegfile{i}, '.eeg'))) %if an EEG file 
       buffsize = eeg.parm.buffersize(i); timeunit = eeg.parm.timeunit(i);
       [timestamp, dat, gain, fs] = ReadEEGFile(eeg.general.eegfile{i}, timeunit, buffsize); % timestamp in second, dat in mV 
    elseif (~isempty(strfind(eeg.general.eegfile{i}, '.spm'))) %%if a binned cluster file
        timestamp = eegdata.cluster0.timepoint{i}; dat = eegdata.cluster0.cnt{i};
        fs = 1/eeg.parm.cl0Binsize(i);
    end
    dat = dat-mean(dat);
    
    normstart = eeg.parm.specNormMinFreq(i); 
    normend = eeg.parm.specNormMaxFreq(i); %%%frequency range for power normalization
    minfy = eeg.parm.specMinFreq(i); maxfy = eeg.parm.specMaxFreq(i); maxfy = min([maxfy fs/2]);
    fstep = eeg.parm.specFreqStep(i);
    F = minfy:fstep:maxfy; %%default frequency range, may be changed 
    evkeyword = eeg.parm.specEvtKeyword{i}; evkeytype = eeg.parm.specEvtType{i};
    
    %%%%power spectrum density window by window for a session
    ws = eeg.parm.specWinSize(i);
    if ws > 0
       wshift = eeg.parm.specWinShift(i);  
       nws = ceil((timestamp(numel(timestamp)) - timestamp(1))/wshift);
       psd = zeros(numel(F), nws);
       %psdnow = cell(1, nws);
       disp('-----------> computing spectrogram for the session...');
       for (tti = 1:nws)
            sT = timestamp(1) + tti*wshift - ws/2; eT = timestamp(1) + tti*wshift + ws/2;
            psd(:,tti) = findmtpower(dat, timestamp, F, fs, sT, eT, 0); %%%new algorithm - multitaper using Slepian sequences
       end
       eeg.spec.timewin{i} = timestamp(1)+(1:nws)'*wshift; 
       eegdata.spec.winpower{i} = abs(psd); %psd[fy, ntime] for each window of a session
    end
    eeg.spec.freq{i} = F; 
   
    %%%%average power for the session
    disp('-----------> computing PSD for the session...');
    eeg.spec.sessPower{i} = findmtpower(dat, timestamp, F, fs, timestamp(1), timestamp(numel(timestamp)), 0);
    eeg.spec.sessNormPower{i} = findnormpower(eeg.spec.sessPower{i}, F, normstart, normend);
    %eeg.spec.sessNormPower{i} = findnormpower(eeg.spec.sessPower{i}, freq, normstart, normend);
    %eeg.spec.sessNormPower{i} = findmtpower(dat, timestamp, F, fs, timestamp(1), timestamp(numel(timestamp)), 1);
    
      %%%%This is the old one that compute PSD for run events
%     %%%run event power 
%     disp('-----------> computing PSD within run events...');
%     evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i};
%     runstarttime = []; runendtime = [];
%     for (j = 1:numel(evTime))
%          if strcmp(evType{j}, 'run')
%              runstarttime = [runstarttime; evTime{j}.start]; runendtime = [runendtime; evTime{j}.ent];
%          end
%     end
%     %eegdata.spec.runevtpower{i} = findmtpower(dat, timestamp, F, fs, runstarttime, runendtime); %%%average power for all run event
%     eeg.spec.runevPower{i} = findmtpower(dat, timestamp, F, fs, runstarttime, runendtime, 0); %%%average power for all run event
%     eeg.spec.runevNormPower{i} = findnormpower(eeg.spec.runevPower{i}, F, normstart, normend);
%     %eeg.spec.runevNormPower{i} = findmtpower(dat, timestamp, F, fs, runstarttime, runendtime, 1); %%%average power for all run event
    if (~isempty(evkeyword))||(~isempty(evkeytype));
       evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i}; evNames = eeg.general.eventname{i};
       evstarttime = []; evendtime = [];
       for (j = 1:numel(evTime))
           evsel = 1;
           if ~isempty(evkeyword)
              if isempty(strfind(lower(evNames{j}), lower(evkeyword)))  evsel = 0; end
           end
           if ~isempty(evkeytype)
               if ~strcmp(evType{j}, evkeytype) evsel = 0; end
           end
           if evsel
              evstarttime = [evstarttime; evTime{j}.start]; evendtime = [evendtime; evTime{j}.ent];
           end
       end
       if (~isempty(evstarttime))
           disp('-----------> computing PSD within selected events...');
           %eegdata.spec.runevtpower{i} = findmtpower(dat, timestamp, F, fs, runstarttime, runendtime); %%%average power for all run event
           %eeg.spec.runevNormPower{i} = findmtpower(dat, timestamp, F, fs, runstarttime, runendtime, 1); %%%average power for all run event
           eeg.spec.evtPower{i} = findmtpower(dat, timestamp, F, fs, evstarttime, evendtime, 0); %%%average power for all run event
           eeg.spec.evtNormPower{i} = findnormpower(eeg.spec.evtPower{i}, F, normstart, normend);
           %eeg.spec.evtNormPower{i} = findnormpower(eeg.spec.evtPower{i}, freq, normstart, normend);
       else
           disp('-----------> no events match with the specified keyword and eventtype');
       end
    end

    
%         %%This is even the older one that computes psd for every event
%     evTime = eegdata.event.eventtimes{i}; evType = eeg.parm.eventtype{i}; evNames = eeg.general.eventname{i};
%     for (j = 1:numel(evTime))
%         disp(['-----------> computing PSD for event file:', evNames{j}]);
%         eegdata.spec.evtpower{i}{j} = findmtpower(dat, timestamp, F, fs, evTime{j}.start, evTime{j}.ent); %%%average power for each event
%     end
   end
end

function power = findmtpower(dat, timestamp, F, fs, startT, endT, normflag) %multi-taper power spectral density estimate - using Chronux package
power = zeros(size(F)); pp = zeros(size(F)); 
%power = []; freq = []; pp = [];
nev = 0; winleng = fs*10; %shift = fs*4; %do it for every 4x data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% chopped to smaller segments because: (1) long data length takes too long to compute Slepian sequecnes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                      (2) like Welch method tend to reduce variance; but
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                      (3) reduce resolution (widen main lobe)
for (k = 1:numel(startT))
     iii = find( (timestamp>=startT(k)) & (timestamp<endT(k)) ); np = numel(iii);
     nwin = ceil(np/winleng); datnow = dat(iii) - mean(dat(iii));
     for (i = 1:nwin)
         datok = datnow( 1+(i-1)*winleng : min([numel(datnow) i*winleng]) ); datok = datok - mean(datok);
         if (numel(datok) >= 0.5*fs) %a hard parameter imposed here for the minimum number of data points ~500ms EEg data
             nev = nev + 1; nnw = 4; 
              nft = 2^nextpow2(numel(datok)); 
              [ppnow,freq] = pmtm(datok, nnw, nft, round(fs)); %%%use 2*nnw-1 tapers 
             %pp = pp + pmtm(datok, nnw, F, round(fs)); %%%use 2*nnw-1 tapers ------ This makes the computation too slow, results very similar
              pp = pp + sumuppower(freq, ppnow, F, normflag);
         end
     end
end
if (nev >0) power = pp/nev; end


% function power = findpower(dat, timestamp, F, fs, startT, endT) %regular fft power spectral density estimate
% power = zeros(size(F)); pp = zeros(size(F)); nev = 0; leng = fs*10; shift = fs*5;
% for (k = 1:numel(startT))
%      iii = find( (timestamp>=startT(k)) & (timestamp<endT(k)) ); np = numel(iii); 
%      nwin = ceil(np/shift); datnow = dat(iii) - mean(dat(iii));
%      %disp(['------> timestamp (first, last): ', num2str(timestamp(1)), '  ', num2str(timestamp(numel(timestamp)))]);
%      %disp(['------> start/end time: ', num2str(startT(k)), '  ', num2str(endT(k))]);
%      for (i = 1:nwin)
%          %disp('-----> get in now');
%          datok = datnow( max([1 i*shift-round(leng/2)]) : min([numel(datnow) i*shift+round(leng/2)]) ); datok = datok - mean(datok);
%          if (numel(datok) >= 0.5*fs) %a hard parameter imposed here for the minimum number of data points ~500ms EEg data
%              nev = nev + 1; nft = 2^nextpow2(numel(datok));
%              freq = fs/2*linspace(0,1,nft/2+1); fftnow = fft(datok, nft)/np; ppnow = 2*abs(fftnow(1:nft/2+1));  %a regular fft
%              pp = pp + sumuppower(freq, ppnow, F);
%          end
%      end
%      
% %      if (np >= 100) %a hard parameter imposed here for the minimum number of data points ~50ms EEg data
% %          datnow = dat(iii) - mean(dat(iii)); nev = nev + 1; nft = min([4096 2^nextpow2(np)]);
% %          freq = fs/2*linspace(0,1,nft/2+1); fftnow = fft(datnow, nft)/np; ppnow = 2*abs(fftnow(1:nft/2+1));  %a regular fft
% %          pp = pp + sumuppower(freq, ppnow, F);
% %      end
% end
% if (nev >0) power = pp/nev; end

function power = sumuppower(freq, pp, F, normflag)
power = zeros(size(F)); %%%need to find the corresponding interval between freq and FF and do the normalization if necessary
fgap = freq(2)-freq(1); Fgap = F(2)-F(1);
if normflag
   pp = pp/(fgap*sum(pp)); %%%normalized by the total power
end
if ~isempty(F)
   [tt, kk] = min( abs(freq-F(1)) ); 
   if (~isempty(kk))
       power(1) = pp(kk);  
   end
   for (i = 2:numel(F))
       kk = find( (freq>F(i-1)) & (freq<=F(i)) ); 
       if (~isempty(kk))
           power(i) = mean(pp(kk)); %*numel(kk)*fgap/Fgap; %%%make sure total power conserved for the same frequency bin --- IT IS CONSERVED.
       end
   end
end

function NormPower = findnormpower(Power, F, normstart, normend)
%NormPower = Power;
iii = find( (F>=normstart) & (F<=normend) );
NormPower = Power/sum(Power(iii))/(F(2)-F(1)); %%%The sumxinterval of (iii-bins) = 1




