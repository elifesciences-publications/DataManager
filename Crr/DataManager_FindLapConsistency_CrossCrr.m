function DataManager_FindLapConsistency_CrossCrr
hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); tagmark = get(gcbo, 'Tag');
hgroup = getappdata(hf, 'hgroup');hfield = getappdata(hf, 'hfield');
plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
[writefilename, okk] = getoutputfile(hf, ow);
if okk
   %get selected cellind
   groupselection = getappdata(hgroup, 'selection'); cellind = []; grpind = find(groupselection == 1); 
   for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end   
end
if okk
   if ~isfield(pinfo, 'crr')
       disp(['--------> not a crrdb; lap crr not computed; aborted.']); okk = 0;
   end
end
if okk 
    if(~isempty(cellind))
        [fname, pname] = uigetfile(fullfile(cd, '*.spikedb'), 'Select a matching spike database:');
        if (numel(fname)>1)
            S = load(fullfile(pname, fname), '-mat'); spikepinfo = S.pinfo;
        else
            okk = 0; disp('--------------> no matching spikedb selected; aborted');
        end 
        if okk
           [pinfo,data] = DoItNow(pinfo,data,spikepinfo,cellind);
            %%%%%%%%save and plot the new database
            save(writefilename, 'pinfo', 'data');
            if (ow)
               iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
               fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
               newname = fname(1:ftt-1); set(hf, 'Name', newname); hmain = hf; 
            else
               hmain = DataManager_DataManager_Callback; plotparm.linkspike = 1; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
            end
            setappdata(hmain,'plotparm', plotparm);
            DataManager_PlotSpikeDatabase(hmain, pinfo, data, 'CellPair', 'clname', '.crrdb'); 
            set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
            setappdata(hmain, 'pinfo', pinfo); setappdata(hmain, 'data', data); pinfo = []; data = []; 
        end
    else
        disp('--------------> no groups selected or groups do not contain any cells');
    end
end
disp('**********************');

function [pinfo,data] = DoItNow(pinfo,data,spikepinfo,cellind)
%%compute lap-by-lap pair-wise cross-correlations: cellind here is the pairind
%%fields assigned here:
if (~isempty(cellind))
  nspike = numel(pinfo.general.animalname); %%%here is the total number of pairs
  if (~isfield(pinfo, 'lapcrr')) pinfo.lapcrr = []; end
  if (~isfield(pinfo.lapcrr, 'evtName')) pinfo.lapcrr.evtName = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evtMeanRate1')) pinfo.lapcrr.evtMeanRate1 = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evtMeanRate2')) pinfo.lapcrr.evtMeanRate2 = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapcrrMeanR')) pinfo.lapcrr.evLapcrrMeanR = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMeanIntcrr')) pinfo.lapcrr.evLapMeanIntcrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean1stPcrr')) pinfo.lapcrr.evLapMean1stPcrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean1stPZsc')) pinfo.lapcrr.evLapMean1stPZsc = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean1stPtime')) pinfo.lapcrr.evLapMean1stPtime = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean2ndPcrr')) pinfo.lapcrr.evLapMean2ndPcrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean2ndPZsc')) pinfo.lapcrr.evLapMean2ndPZsc = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean2ndPtime')) pinfo.lapcrr.evLapMean2ndPtime = cell(1, nspike); end 
  
  if (~isfield(pinfo.lapcrr, 'evLapSlidecrrMeanR')) pinfo.lapcrr.evLapSlidecrrMeanR = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMeanIntSlidecrr')) pinfo.lapcrr.evLapMeanIntSlidecrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean1stPSlidecrr')) pinfo.lapcrr.evLapMean1stPSlidecrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean1stPSlideZsc')) pinfo.lapcrr.evLapMean1stPSlideZsc = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean1stPSlidetime')) pinfo.lapcrr.evLapMean1stPSlidetime = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean2ndPSlidecrr')) pinfo.lapcrr.evLapMean2ndPSlidecrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean2ndPSlideZsc')) pinfo.lapcrr.evLapMean2ndPSlideZsc = cell(1, nspike); end 
  if (~isfield(pinfo.lapcrr, 'evLapMean2ndPSlidetime')) pinfo.lapcrr.evLapMean2ndPSlidetime = cell(1, nspike); end 
end
sss = [num2str(numel(cellind)) ' pairs found. Compute lap-by-lap crr?']; 
aa = questdlg(sss);
if (strcmp(aa, 'Yes'))
for (jjjjk = 1:numel(cellind))
    i = cellind(jjjjk); crrtypenow = pinfo.general.crrtype{i};
%if (strcmp(pinfo.general.crrtype{i}, 'cross'))
    disp(['-----------> compute lap crr (', num2str(jjjjk), ' out of ', num2str(numel(cellind)), ') -- ', pinfo.general.clname{i}]);
    %%%%%%get all parameters
    crrmode = pinfo.parm.crrmode{i}; binsize = pinfo.parm.timebin(i); maxlag = pinfo.parm.maxlag(i); 
    timeunit = pinfo.parm.timeunit(i); pairind = data.crr.cellind{i};
    intbin = pinfo.parm.intbin{i}; crrtype = pinfo.general.crrtype{i};
    smoothmode = pinfo.parm.smoothmode{i}; smoothbin = pinfo.parm.smoothbin(i);
    P1searchwin = pinfo.parm.P1stSearchWin{i}; P2searchwin = pinfo.parm.P2ndSearchWin{i}; 
    searchmode = pinfo.parm.searchPmode{i}; %%%%%%either 'peack', 'trough', or 'both'
    setbacktime = pinfo.parm.setbacktime(i);
    %%%%rectify spikes all into row vectors %%%%using default column vectors does not work well, since matlab treat scalar as row vectors
    spike1 = data.spike.spiketime{pairind(1)}*timeunit; 
    [A, B] = size(spike1); if (B == 1) spike1 = spike1'; end
    spike2 = data.spike.spiketime{pairind(2)}*timeunit;
    [A, B] = size(spike2); if (B == 1) spike2 = spike2'; end
    %%%%set up time bins
    np = ceil(maxlag/binsize); lagbin = (-np:1:np); %%%%lagbin is a row vector
    P1win = P1searchwin/binsize; P1winind = find( (lagbin>=P1win(1)) & (lagbin<=P1win(2)) );
    %%%%event crr 
    evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i};  evType = pinfo.parm.eventtype{i}; nev = 0;
    for (j = 1:numel(evTime))
        if strcmp(evType{j}, 'run')
           nev = nev + 1; pinfo.lapcrr.evtName{i}{nev} = evName{j}; 
           pinfo.lapcrr.evtMeanRate1{i}{nev} = spikepinfo.firing.evtmeanrate{pairind(1)}(j);
           pinfo.lapcrr.evtMeanRate2{i}{nev} = spikepinfo.firing.evtmeanrate{pairind(2)}(j);
           sT = evTime{j}.start-setbacktime; eT = evTime{j}.ent+setbacktime; nlap = numel(sT);
           crr = zeros(nlap, numel(lagbin));  
           for (k = 1:nlap)
               crrnow = computecrrhere(spike1, spike2, sT(k), eT(k), lagbin, crrmode, binsize);
               if (strncmpi(smoothmode, 'yes', 1) & strcmp(crrtypenow, 'cross')) crrnow = smoothcrr(crrnow, smoothbin); end
               crr(k,:) = crrnow;
           end
           rr = []; [~,nlag] = size(crr);
           for (ki = 1:nlap)
                for (kj = ki+1:nlap)
                    ccc = computerrcorr(crr(ki,P1winind), crr(kj,P1winind)); %%%both row vectors - only do this for the central window
                    rr = [rr ccc];
                end
           end
           pinfo.lapcrr.evLapcrrMeanR{i}(nev) = mean(rr(~isnan(rr)));
           crrnow = NaN*zeros(1, nlap);
           for (ki = 1:nlag)
                valnow = crr(:,ki); crrnow(ki) = mean(valnow(~isnan(valnow)));
           end
           [intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime] = determinepeaks(crrnow, lagbin, intbin, binsize, P1searchwin, P2searchwin, searchmode);
           pinfo.lapcrr.evLapMeanIntcrr{i}(nev) = intcrr;
           pinfo.lapcrr.evLapMean1stPcrr{i}(nev) = Ppeakcrr; pinfo.lapcrr.evLapMean1stPZsc{i}(nev) = PZsc; pinfo.lapcrr.evLapMean1stPtime{i}(nev) = Ppeaktime; 
           pinfo.lapcrr.evLapMean2ndPcrr{i}(nev) = Speakcrr; pinfo.lapcrr.evLapMean2ndPZsc{i}(nev) = SZsc; pinfo.lapcrr.evLapMean2ndPtime{i}(nev) = Speaktime; 
           %%%%%compute a sliding shuffled version
           crr = zeros(nlap, numel(lagbin));  
           for (k = 1:nlap)
               crrnow = computeslidingcrrhere(spike1, spike2, sT(k), eT(k), lagbin, crrmode, binsize);
               if (strncmpi(smoothmode, 'yes', 1)) crrnow = smoothcrr(crrnow, smoothbin); end
               crr(k,:) = crrnow;
           end
           rr = [];
           for (ki = 1:nlap)
                for (kj = ki+1:nlap)
                    ccc = computerrcorr(crr(ki,P1winind), crr(kj,P1winind)); %%%both row vectors - only do this for the central window
                    rr = [rr ccc];
                end
           end
           pinfo.lapcrr.evLapSlidecrrMeanR{i}(nev) = mean(rr(~isnan(rr)));
           crrnow = NaN*zeros(1, nlap);
           for (ki = 1:nlag)
                valnow = crr(:,ki); crrnow(ki) = mean(valnow(~isnan(valnow)));
           end
           [intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime] = determinepeaks(crrnow, lagbin, intbin, binsize, P1searchwin, P2searchwin, searchmode);
           pinfo.lapcrr.evLapMeanIntSlidecrr{i}(nev) = intcrr;
           pinfo.lapcrr.evLapMean1stPSlidecrr{i}(nev) = Ppeakcrr; pinfo.lapcrr.evLapMean1stPSlideZsc{i}(nev) = PZsc; pinfo.lapcrr.evLapMean1stPSlidetime{i}(nev) = Ppeaktime; 
           pinfo.lapcrr.evLapMean2ndPSlidecrr{i}(nev) = Speakcrr; pinfo.lapcrr.evLapMean2ndPSlideZsc{i}(nev) = SZsc; pinfo.lapcrr.evLapMean2ndPSlidetime{i}(nev) = Speaktime;
        end
    end
end
end


function [intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime] = determinepeaks(crr, lagbin, intbin, binsize, P1searchwin, P2searchwin, searchmode)
Ppeakcrr = NaN; Ppeaktime = NaN; Speakcrr = NaN; Speaktime = NaN; PZsc = NaN; SZsc = NaN; zcrr = zscore(crr);
%%%%for intcrr
intbin = intbin/binsize; intwin = find( (lagbin>=intbin(1)) & (lagbin<=intbin(2)) );
centercrr = crr(intersect(intwin, 1:numel(crr)));
intcrr = mean(centercrr(~isnan(centercrr))); 
%%%%find primary/secondary peaks
devalue = NaN; crr = Position_Dezero(crr, devalue); %%%interpolate NaN values
[minindex, maxindex] = FindLocal(crr);

%%%%noe analyze search windows and search mode
if (strncmpi(searchmode, 'peak', 2))
   allind = intersect(maxindex+1, 1:numel(crr)); %%%somehow the identified indices are shiftd by -1 
elseif (strncmpi(searchmode, 'trough', 2))
   allind = intersect(minindex+1, 1:numel(crr)); %%%somehow the identified indices are shiftd by -1 
else
   allind = intersect(union(minindex+1, maxindex+1), 1:numel(crr)); %%%somehow the identified indices are shiftd by -1
end
P1searchwin = P1searchwin/binsize; P1win = find( (lagbin>=P1searchwin(1)) & (lagbin<=P1searchwin(2)) );
P2searchwin = P2searchwin/binsize; P2win = find( (lagbin>=P2searchwin(1)) & (lagbin<=P2searchwin(2)) );
if (~isempty(allind))
   %%%%%%%%primary 1st peak
   all1stind = intersect(allind, P1win); iii = [];
   if (~isempty(all1stind))
      if (strncmpi(searchmode, 'peak', 2))
         [~, iii] = max(crr(all1stind));
      elseif (strncmpi(searchmode, 'trough', 2))
         [~, iii] = min(crr(all1stind));
      else
         [~, iii] = max(abs(crr(all1stind)));
      end
      Ppeakcrr = crr(all1stind(iii)); Ppeaktime = binsize*lagbin(all1stind(iii)); PZsc = zcrr(all1stind(iii));
   end
   %%%secondary peak
   all2ndind = intersect(setdiff(allind, all1stind(iii)), P2win);
   if (~isempty(all2ndind))
      if (strncmpi(searchmode, 'peak', 2))
         [~, iii] = max(crr(all2ndind));
      elseif (strncmpi(searchmode, 'trough', 2))
         [~, iii] = min(crr(all2ndind));
      else
         [~, iii] = max(abs(crr(all2ndind)));
      end
      Speakcrr = crr(all2ndind(iii)); Speaktime = binsize*lagbin(all2ndind(iii)); SZsc = zcrr(all1stind(iii));
   end
end

function newcrr = smoothcrr(crr, smoothbin)
newcrr = NaN*ones(size(crr)); %zeros(size(crr)); %
halfbin = floor(smoothbin/2);
for (i = 1:numel(crr))
    kkk = [i-halfbin:i+halfbin]; kkk = kkk( (kkk>=1) & (kkk<=numel(crr)) );
    crrnow = crr(kkk); newcrr(i) = mean(crrnow(~isnan(crrnow)));
end

function [crr, norm] = computeslidingcrrhere(spike1, spike2, sT, eT, lagbin, crrmode, binsize)
crr = NaN*ones(size(lagbin)); %zeros(size(lagbin)); % 
norm = 0;
if (~isempty(spike1)) && (~isempty(spike2))
    iii = find(eT-sT>=binsize); sT = sT(iii); eT = eT(iii); %%%filter out short events
    %%%sliding spike2 here
    spike1 = slideshuffle(spike1, sT, eT);
    spike2 = slideshuffle(spike2, sT, eT); %%%%here in this particular function  numel(sT) = numel(eT) = 1: just 1 lap
    if (strncmpi(crrmode, 'rate', 2))
        crr = computeratecrr(spike1, spike2, sT, eT, lagbin, binsize);
    elseif (strncmpi(crrmode, 'count', 2))
        crr = computecountcrr(spike1, spike2, sT, eT, lagbin, binsize);
    elseif (strncmpi(crrmode, 'normcount', 2))
        [crr, norm] = computenormcountcrr(spike1, spike2, sT, eT, lagbin, binsize);
    end
end

function datanow = slideshuffle(spikedata, tmin, tmax)
spikedata = spikedata( (spikedata>=tmin) & (spikedata<=tmax) );
tintv = tmax-tmin; datanow = spikedata-tmin; 
datanow = datanow + rand*(tmax-tmin); datanow = mod(datanow, tintv) + tmin;

function [crr, norm] = computecrrhere(spike1, spike2, sT, eT, lagbin, crrmode, binsize)
crr = NaN*ones(size(lagbin)); %zeros(size(lagbin)); 
norm = 0;
if (~isempty(spike1)) && (~isempty(spike2))
    iii = find(eT-sT>=binsize); sT = sT(iii); eT = eT(iii); %%%filter out short events
    if (strncmpi(crrmode, 'rate', 2))
        crr = computeratecrr(spike1, spike2, sT, eT, lagbin, binsize);
    elseif (strncmpi(crrmode, 'count', 2))
        crr = computecountcrr(spike1, spike2, sT, eT, lagbin, binsize);
    elseif (strncmpi(crrmode, 'normcount', 2))
        [crr, norm] = computenormcountcrr(spike1, spike2, sT, eT, lagbin, binsize);
    end
end

function ccc = computeratecrr(t1, t2, sT, eT, lagbin, binsize) %%%t1, t2 are row vectors
nev = numel(sT); timebin = cell(1,nev); nc = numel(lagbin); ccc = NaN*ones(size(lagbin)); %ccc = zeros(size(lagbin));  
cnt1 = cell(1,nev); cnt2 = cell(1,nev);
for (i = 1:nev)
     nbin = ceil( (eT(i)-sT(i))/binsize );
     timebin{i} = sT(i) + ( (1:nbin) - 1 )*binsize;
end
%tttt = cputime;
for (i = 1:nev)
    t1now = t1( (t1>=sT(i)) & (t1<=eT(i)) ); t2now = t2( (t2>=sT(i)) & (t2<=eT(i)) ); %%%filter here for saving hisc time?
    nbinnow = numel(timebin{i}) - 1; 
    cnt1{i} = zeros(1, nbinnow); cnt2{i} = zeros(nbinnow, 1);
    if (~isempty(t1now))
       cnt = histc(t1now, timebin{i}); cnt1{i} = cnt(1:nbinnow); %output row vector (histc output same row or column vector as t1now)
    end
    if (~isempty(t2now))
       cnt = histc(t2now, timebin{i}); cnt2{i} = cnt(1:nbinnow)'; %output column vector
    end
end
%tttt1 = cputime;
%disp(strcat('*****multiple ep binning time:*****', num2str(tttt1-tttt)));
%move 2nd t2 firing rate time series accroding to time lag: if lagbin <0, move t2 to the right against t1
for (i= 1:nc)
    %tttt = cputime;
    nbin = 0; sum1 = 0; sum2 = 0; sumsq1 = 0; sumsq2 = 0; sumsq12 = 0; %catecated count vector cross all episodes for corrcoeff
    for (j = 1:nev)
        cn1 = cnt1{j}; cn2 = cnt2{j}; %%%%using matrix instead of cell saves a lot of time
        nbinnow = numel(timebin{j}) - 1; 
        startbin1 = max([1 1-lagbin(i)]);
        endbin1 = min([nbinnow nbinnow-lagbin(i)]);
        startbin2 = max([1 1+lagbin(i)]);
        endbin2 = min([nbinnow nbinnow+lagbin(i)]);
        if ( (startbin1<endbin1) & (startbin2<endbin2) ) %if both series have space for the lag in this episode
            nbin = nbin + endbin1 - startbin1 + 1;
            sum1 = sum1 + sum(cn1(startbin1:endbin1)); sum2 = sum2 + sum(cn2(startbin2:endbin2));
            sumsq1 = sumsq1 + sum(cn1(startbin1:endbin1).^2); sumsq2 = sumsq2 + sum(cn2(startbin2:endbin2).^2);
            sumsq12 = sumsq12 + cn1(startbin1:endbin1) * cn2(startbin2:endbin2); %%%row vector times column vector
        end
    end
    %tttt1 = cputime;
    %disp(strcat('*****crr single point computing time:*****', num2str(tttt1-tttt)));
    if (nbin ~= 0)
       co1 = sumsq1-sum1^2/nbin; co2 = sumsq2-sum2^2/nbin;
       if ( (co1 == 0) | (co2 == 0) )
           ccc(i) = NaN; %0;
       else
           ccc(i) = (sumsq12 - sum1*sum2/nbin) / sqrt(co1*co2);
       end
    end
end

function count = computecountcrr(spike1, spike2, sT, eT, lagbin, bin)
count = zeros(size(lagbin));   %initial assignment for spike coun
np = lagbin(numel(lagbin)); %%%number of left and right lag points
t1 = []; t2 = [];
for (i = 1:numel(sT))
    t1 = [t1 spike1( (spike1>=sT(i)) & (spike1<=eT(i)) )];
    t2 = [t2 spike2( (spike2>=sT(i)) & (spike2<=eT(i)) )]; %%%do not filter the second train to avoid the edge effect?
                                                            %%%better do it to reflect the true event windows;
                                                            %%%   be aware of the edge effect for short event files
end
for (i = 1:numel(t1))
    % for each spike in first train, count spike number in second train that within a time window defined by bin
%     for (k = -np:1:np)
%         minn = t1(i) + k*bin - bin/2;  %defining the window
%         maxx = t1(i) + k*bin + bin/2;
%         count(k+np+1) = count(k+np+1) + numel(find( (t2>=minn) & (t2<maxx)));
%     end
    k = -np:1:np; wind = t1(i) + k*bin - bin/2;
    count(k+np+1) = count(k+np+1) + histc(t2, wind);
end

function [count, mm] = computenormcountcrr(spike1, spike2, sT, eT, lagbin, bin)
%%%%normalized count: assuming spike1 and spike2 are independent poisson, at each time lag bin, the mean expected spike count is
%%%%       mm = numel(spike1)*numel(spike2)*bin/totallengthofspiketrains; 
%%%% and the std (for Poisson) is ss = sqrt(mm)
%%%% normalized as (count-mm)/ss
%%%% lagbin = (-np:1:np); 
count = zeros(size(lagbin)); mm = 0;   %initial assignment for spike coun = row vector
np = lagbin(numel(lagbin)); %%%number of left and right lag points
t1 = []; t2 = []; tleng = 0;
for (i = 1:numel(sT))
    tleng = tleng + (eT(i)-sT(i));
    t1 = [t1 spike1( (spike1>=sT(i)) & (spike1<=eT(i)) )];
    t2 = [t2 spike2( (spike2>=sT(i)) & (spike2<=eT(i)) )]; %%%%need to filter to compute the norm factor
end
%tttt1 = cputime;
if (~isempty(t2))
for (i = 1:numel(t1)) % for each spike in first train, count spike number in second train that within a time window defined by bin
%     for (k = -np:1:np)
%         minn = t1(i) + k*bin - bin/2;  %defining the window
%         maxx = t1(i) + k*bin + bin/2;
%         count(k+np+1) = count(k+np+1) + numel(find( (t2>=minn) & (t2<maxx)));
%     end
      k = -np:1:np; wind = t1(i) + k*bin - bin/2;
      count(k+np+1) = count(k+np+1) + histc(t2, wind);
end
%tttt2 = cputime;
%disp(strcat('*****crr single point counting time:*****', num2str(tttt2-tttt1)));
%%%%%normalization
mm = numel(t1)*numel(t2)*bin/tleng;
if (mm>0) count = (count-mm)/sqrt(mm); end
end


function [minindex, maxindex] = FindLocal(dat)
%%NOW IT WORKS REALLY WELL!!!!
%%the problem with the derivative polarity change search is that EEG traces
%%are not that smooth, it gets more rugged on peaks and troughs
%%solved by smoothing EEG traces
deriv = diff(dat); %difference of EEG points
deriva = deriv(1:numel(deriv)-1);
derivb = deriv(2:numel(deriv)); %%deriva and derivb shift by one point
pola = deriva .* derivb; %chech for polarity change (product <0)
mindex = find(pola <= 0); %mindex are indices in pola, deriva and derivb, mindex+1 is the real singular point index in dat
%disp(strcat('-----> totally:', num2str(numel(mindex)), ' found'));
if (numel(mindex) < 2)
    minindex = []; maxindex = []; %if only o or 1 singular point: error
else
    %for each polarity change, check if minima or maxima
    mvalue = ClassPol(mindex, deriva);
    %now classify singular values to minima and maxima
    maxindex = find(mvalue == 1); %return real maximum indices in dat
    minindex = find(mvalue == -1); %return real minimum indices in dat
end
nm = numel(minindex); nM = numel(maxindex);
if (abs(nm-nM)>1) disp('-----------> Warning: minima and maxima do not match!'); end
nn = min([nm nM]); minindex = minindex(1:nn); maxindex = maxindex(1:nn);

function mvalue = ClassPol(mindex, deriva)
%have to classify in detail to get an accurate peaks and troughs
%mvalue = 1 if maxima, =-1 if minima, otherwise =0
%strategy: deal with point one by one, slope1~=0, then move to next until another
%slope2~=0, mvalue and position decide by slope1*slope2
mvalue = zeros(1, numel(deriva)); %initial assignment =0
pointnow = 1;
pointend = numel(mindex);
while (pointnow < pointend)
    slope1 = deriva(mindex(pointnow));
    zeropoint = 0; %how many zeros on peaks or troughs
    if (slope1 > 0) %if up maximum candidate
        while (mindex(pointnow)+1+zeropoint < numel(deriva))
            if (deriva(mindex(pointnow)+1+zeropoint) == 0) %if next point flat (zero) go on getting next point
                zeropoint = zeropoint + 1;
            else
                break
            end
        end
        if (mindex(pointnow)+1+zeropoint < numel(deriva)) %if still valid index
            if (deriva(mindex(pointnow)+1+zeropoint) < 0) %if next none-zero point is down then a maximum
               mvalue(mindex(pointnow)+floor((zeropoint+1)/2)) = 1; %take middle point in the zero points
            end
        end
    elseif (slope1 < 0) %if down minimum candidate
        while (mindex(pointnow)+1+zeropoint < numel(deriva))
            if (deriva(mindex(pointnow)+1+zeropoint) == 0) %if next point flat (zero) go on getting next point
                zeropoint = zeropoint + 1;
            else
                break
            end
        end
        if (mindex(pointnow)+1+zeropoint < numel(deriva)) %if still valid index
            if (deriva(mindex(pointnow)+1+zeropoint) > 0) %if next none-zero point is down then a maximum
                mvalue(mindex(pointnow)+floor((zeropoint+1)/2)) = -1; %take middle point in the zero points
            end
        end
    end
    pointnow = pointnow + zeropoint + 1;
end

function [writefilename, okk] = getoutputfile(hf, ow)
okk = 1; writefilename = [];
if (ow == 0)
   [fname, pname] = uiputfile(fullfile(cd, '*.crrdb'), 'Write the new crr database to:');
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

function ccc = computerrcorr(rr1, rr2) %%%both row vectors
ind = find( (~isnan(rr1)) & (~isnan(rr2)) ); rr1 = rr1(ind); rr2 = rr2(ind);
sum1 = sum(rr1); sum2 = sum(rr2); nbin = numel(rr1);
sumsq1 = sum(rr1.^2); sumsq2 = sum(rr2.^2);
sumsq12 = rr1 * rr2'; %%%row vector times column vector
co1 = sumsq1-sum1^2/nbin; co2 = sumsq2-sum2^2/nbin;
if ( (co1 == 0) | (co2 == 0) )
      ccc = NaN;
else
      ccc = (sumsq12 - sum1*sum2/nbin) / sqrt(co1*co2);
end
