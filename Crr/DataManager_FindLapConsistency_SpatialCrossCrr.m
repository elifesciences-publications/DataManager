function DataManager_FindLapConsistency_SpatialCrossCrr
hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data'); tagmark = get(gcbo, 'Tag');
hgroup = getappdata(hf, 'hgroup'); hfield = getappdata(hf, 'hfield');
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
    if isempty(cellind)
       disp('--------------> no groups selected or groups do not contain any cells');
    else
       if (plotparm.linkbehav == 0);
          disp(['--------> no behav data linked']);
       else
          behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
          [fname, pname] = uigetfile(fullfile(cd, '*.spikedb'), 'Select a matching spike database:');
          if (numel(fname)>1)
             S = load(fullfile(pname, fname), '-mat'); spikepinfo = S.pinfo;
          else
             okk = 0; disp('--------------> no matching spikedb selected; aborted');
          end
          if okk
            [pinfo,data] = DoItNow(pinfo,data,behav,bhdata,spikepinfo,cellind);
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
       end
    end
end
disp('**********************');

function [pinfo,data] = DoItNow(pinfo,data,behav, bhdata, spikepinfo,cellind)
%%compute lap-by-lap pair-wise cross-correlations: cellind here is the pairind
%%fields assigned here:
if (~isempty(cellind))
  nspike = numel(pinfo.general.animalname); %%%here is the total number of pairs
  if (~isfield(pinfo.parm, 'spatMaxlag')) pinfo.parm.spatMaxlag = 250*ones(1, nspike); end 
  if (~isfield(pinfo, 'lapspatialcrr')) pinfo.lapspatialcrr = []; end
  if (~isfield(pinfo.lapspatialcrr, 'evtName')) pinfo.lapspatialcrr.evtName = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evtMeanRate1')) pinfo.lapspatialcrr.evtMeanRate1 = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evtMeanRate2')) pinfo.lapspatialcrr.evtMeanRate2 = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapcrrMeanR')) pinfo.lapspatialcrr.evLapcrrMeanR = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMeanIntcrr')) pinfo.lapspatialcrr.evLapMeanIntcrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean1stPcrr')) pinfo.lapspatialcrr.evLapMean1stPcrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean1stPZsc')) pinfo.lapspatialcrr.evLapMean1stPZsc = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean1stPloc')) pinfo.lapspatialcrr.evLapMean1stPloc = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean2ndPcrr')) pinfo.lapspatialcrr.evLapMean2ndPcrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean2ndPZsc')) pinfo.lapspatialcrr.evLapMean2ndPZsc = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean2ndPloc')) pinfo.lapspatialcrr.evLapMean2ndPloc = cell(1, nspike); end 
  %if (~isfield(pinfo.lapspatialcrr, 'evLapMean3rdPcrr')) pinfo.lapspatialcrr.evLapMean3rdPcrr = cell(1, nspike); end 
  %if (~isfield(pinfo.lapspatialcrr, 'evLapMean3rdPZsc')) pinfo.lapspatialcrr.evLapMean3rdPZsc = cell(1, nspike); end 
  %if (~isfield(pinfo.lapspatialcrr, 'evLapMean3rdPloc')) pinfo.lapspatialcrr.evLapMean3rdPloc = cell(1, nspike); end 
  
  if (~isfield(pinfo.lapspatialcrr, 'evLapshufcrrMeanR')) pinfo.lapspatialcrr.evLapshufcrrMeanR = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMeanIntshufcrr')) pinfo.lapspatialcrr.evLapMeanIntshufcrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean1stPshufcrr')) pinfo.lapspatialcrr.evLapMean1stPshufcrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean1stPshufZsc')) pinfo.lapspatialcrr.evLapMean1stPshufZsc = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean1stPshufloc')) pinfo.lapspatialcrr.evLapMean1stPshufloc = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean2ndPshufcrr')) pinfo.lapspatialcrr.evLapMean2ndPshufcrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean2ndPshufZsc')) pinfo.lapspatialcrr.evLapMean2ndPshufZsc = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean2ndPshufloc')) pinfo.lapspatialcrr.evLapMean2ndPshufloc = cell(1, nspike); end
  %if (~isfield(pinfo.lapspatialcrr, 'evLapMean3rdPshufcrr')) pinfo.lapspatialcrr.evLapMean3rdPshufcrr = cell(1, nspike); end 
  %if (~isfield(pinfo.lapspatialcrr, 'evLapMean3rdPshufZsc')) pinfo.lapspatialcrr.evLapMean3rdPshufZsc = cell(1, nspike); end 
  %if (~isfield(pinfo.lapspatialcrr, 'evLapMean3rdPshufloc')) pinfo.lapspatialcrr.evLapMean3rdPshufloc = cell(1, nspike); end
  
  if (~isfield(pinfo.lapspatialcrr, 'evLapSlidecrrMeanR')) pinfo.lapspatialcrr.evLapSlidecrrMeanR = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMeanIntSlidecrr')) pinfo.lapspatialcrr.evLapMeanIntSlidecrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean1stPSlidecrr')) pinfo.lapspatialcrr.evLapMean1stPSlidecrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean1stPSlideZsc')) pinfo.lapspatialcrr.evLapMean1stPSlideZsc = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean1stPSlideloc')) pinfo.lapspatialcrr.evLapMean1stPSlideloc = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean2ndPSlidecrr')) pinfo.lapspatialcrr.evLapMean2ndPSlidecrr = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean2ndPSlideZsc')) pinfo.lapspatialcrr.evLapMean2ndPSlideZsc = cell(1, nspike); end 
  if (~isfield(pinfo.lapspatialcrr, 'evLapMean2ndPSlideloc')) pinfo.lapspatialcrr.evLapMean2ndPSlideloc = cell(1, nspike); end
  %if (~isfield(pinfo.lapspatialcrr, 'evLapMean3rdPSlidecrr')) pinfo.lapspatialcrr.evLapMean3rdPSlidecrr = cell(1, nspike); end 
  %if (~isfield(pinfo.lapspatialcrr, 'evLapMean3rdPSlideZsc')) pinfo.lapspatialcrr.evLapMean3rdPSlideZsc = cell(1, nspike); end 
  %if (~isfield(pinfo.lapspatialcrr, 'evLapMean3rdPSlideloc')) pinfo.lapspatialcrr.evLapMean3rdPSlideloc = cell(1, nspike); end
end
sss = [num2str(numel(cellind)) ' pairs found. Compute lap-by-lap spatial crr?']; 
aa = questdlg(sss);
if (strcmp(aa, 'Yes'))
for (jjjjk = 1:numel(cellind))
    i = cellind(jjjjk);
%if (strcmp(pinfo.general.crrtype{i}, 'cross'))
    disp(['-----------> compute lap spatial crr (', num2str(jjjjk), ' out of ', num2str(numel(cellind)), ') -- ', pinfo.general.clname{i}]);
    %%%%%%get all parameters
    pairind = data.crr.cellind{i}; timeunit = pinfo.parm.timeunit(i);
    %%%finaldirnow = strtok(pinfo.general.finaldir{i}, '_'); %does not work
    aaa = strfind(pinfo.general.finaldir{i}, '__'); 
    finaldirnow = pinfo.general.finaldir{i}(1:aaa(1)-1);
    smParm.d1sigma = spikepinfo.parm.fSmooth1DSigma(pairind(1)); 
    smParm.d1Nsig = spikepinfo.parm.fSmooth1DNSigma(pairind(1)); smParm.d1sm = spikepinfo.parm.fSmooth1D{pairind(1)};
    %smoothmode = pinfo.parm.smoothmode{i}; smoothbin = pinfo.parm.smoothbin(i);
    
    %crrmode = pinfo.parm.crrmode{i}; binsize = pinfo.parm.timebin(i); maxlag = pinfo.parm.maxlag(i); 
    %intbin = pinfo.parm.intbin{i}; crrtype = pinfo.general.crrtype{i};
    
    %P1searchwin = pinfo.parm.P1stSearchWin{i}; P2searchwin = pinfo.parm.P2ndSearchWin{i}; 
    %searchmode = pinfo.parm.searchPmode{i}; %%%%%%either 'peack', 'trough', or 'both'
    %%%%rectify spikes all into row vectors %%%%using default column vectors does not work well, since matlab treat scalar as row vectors
    spike1 = data.spike.spiketime{pairind(1)}*timeunit; 
    [A, B] = size(spike1); if (B == 1) spike1 = spike1'; end
    spike2 = data.spike.spiketime{pairind(2)}*timeunit;
    [A, B] = size(spike2); if (B == 1) spike2 = spike2'; end
    %%%%set up time bins
    %np = ceil(maxlag/binsize); lagbin = (-np:1:np); %%%%lagbin is a row vector
    %P1win = P1searchwin/binsize; P1winind = find( (lagbin>=P1win(1)) & (lagbin<=P1win(2)) );
    %%%%event crr 
    evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i}; nev = 0;
    for (j = 1:numel(evTime))
        if strcmp(evType{j}, 'run')
            %%%%locate event position data
            evSess = identifysession(evTime{j}, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
            posid = []; evid = [];
            if (~isempty(evSess))
                posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess) );
            end
            if numel(posid) == 1
                if (isfield(behav.general, 'eventname'))
                    evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
                else
                    evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
                end
            end
            if (numel(posid)~=1)||(numel(evid)~=1)
                disp(['------------------> Warning: lap rate map consistency not computed for this event: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
            else
                nev = nev + 1; pinfo.lapspatialcrr.evtName{i}{nev} = evName{j}; 
                pinfo.lapspatialcrr.evtMeanRate1{i}{nev} = spikepinfo.firing.evtmeanrate{pairind(1)}(j);
                pinfo.lapspatialcrr.evtMeanRate2{i}{nev} = spikepinfo.firing.evtmeanrate{pairind(2)}(j);
                spatbin = behav.parm.s1dbinsize(posid); spatmaxlag = pinfo.parm.spatMaxlag; nlap = numel(evTime{j}.start);
                xbin = bhdata.event.Xbin{posid}{evid}; D1rate = cell(1, nlap); D2rate = cell(1, nlap); 
                for (jnow = 1:nlap)
                     occutimenow{jnow} = bhdata.event.LapOccuptime{posid}{evid}{jnow}; 
                     lappostime{1} = bhdata.event.LapAllPostimestamp{posid}{evid}{jnow}; lapx{1} = bhdata.event.LapAllX{posid}{evid}{jnow};
                     evok.start = evTime{j}.start(jnow); evok.ent = evTime{j}.ent(jnow);
                     [ratenow, ~] = getlinearmap(spike1, evok, lappostime, lapx, xbin, occutimenow{jnow}, smParm); D1rate{jnow} = ratenow';
                     [ratenow, ~] = getlinearmap(spike2, evok, lappostime, lapx, xbin, occutimenow{jnow}, smParm); D2rate{jnow} = ratenow';
                end
                [crr, T] = computecrrhere(D1rate, D2rate, spatmaxlag, spatbin);
                [rr, intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime, Tpeakcrr, TZsc, Tpeaktime] = findcrrparameters(crr, T, spatbin); %, smoothmode, smoothbin);
                pinfo.lapspatialcrr.evLapcrrMeanR{i}(nev) = rr;
                pinfo.lapspatialcrr.evLapMeanIntcrr{i}(nev) = intcrr;
                pinfo.lapspatialcrr.evLapMean1stPcrr{i}(nev) = Ppeakcrr; pinfo.lapspatialcrr.evLapMean1stPZsc{i}(nev) = PZsc; 
                pinfo.lapspatialcrr.evLapMean1stPloc{i}(nev) = Ppeaktime; 
                pinfo.lapspatialcrr.evLapMean2ndPcrr{i}(nev) = Speakcrr; pinfo.lapspatialcrr.evLapMean2ndPZsc{i}(nev) = SZsc; 
                pinfo.lapspatialcrr.evLapMean2ndPloc{i}(nev) = Speaktime;
                %pinfo.lapspatialcrr.evLapMean3rdPcrr{i}(nev) = Tpeakcrr; pinfo.lapspatialcrr.evLapMean3rdPZsc{i}(nev) = TZsc; 
                %pinfo.lapspatialcrr.evLapMean3rdPloc{i}(nev) = Tpeaktime;
                %%%%%%lap shuffles
                [crr, T] = computecrrhere(D1rate, D2rate(randperm(nlap)), spatmaxlag, spatbin);
                [rr, intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime, Tpeakcrr, TZsc, Tpeaktime] = findcrrparameters(crr, T, spatbin); %, smoothmode, smoothbin);
                pinfo.lapspatialcrr.evLapshufcrrMeanR{i}(nev) = rr;
                pinfo.lapspatialcrr.evLapMeanIntshufcrr{i}(nev) = intcrr;
                pinfo.lapspatialcrr.evLapMean1stPshufcrr{i}(nev) = Ppeakcrr; pinfo.lapspatialcrr.evLapMean1stPshufZsc{i}(nev) = PZsc; 
                pinfo.lapspatialcrr.evLapMean1stPshufloc{i}(nev) = Ppeaktime; 
                pinfo.lapspatialcrr.evLapMean2ndPshufcrr{i}(nev) = Speakcrr; pinfo.lapspatialcrr.evLapMean2ndPshufZsc{i}(nev) = SZsc; 
                pinfo.lapspatialcrr.evLapMean2ndPshufloc{i}(nev) = Speaktime;
                %pinfo.lapspatialcrr.evLapMean3rdPshufcrr{i}(nev) = Tpeakcrr; pinfo.lapspatialcrr.evLapMean3rdPshufZsc{i}(nev) = TZsc; 
                %pinfo.lapspatialcrr.evLapMean3rdPshufloc{i}(nev) = Tpeaktime;
                %%%%%%sliding shuffles
                [crr, T] = computecrrhere(circleslideratecurve(D1rate), circleslideratecurve(D2rate), spatmaxlag, spatbin);
                [rr, intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime, Tpeakcrr, TZsc, Tpeaktime] = findcrrparameters(crr, T, spatbin); %, smoothmode, smoothbin);
                pinfo.lapspatialcrr.evLapSlidecrrMeanR{i}(nev) = rr;
                pinfo.lapspatialcrr.evLapMeanIntSlidecrr{i}(nev) = intcrr;
                pinfo.lapspatialcrr.evLapMean1stPSlidecrr{i}(nev) = Ppeakcrr; pinfo.lapspatialcrr.evLapMean1stPSlideZsc{i}(nev) = PZsc; 
                pinfo.lapspatialcrr.evLapMean1stPSlideloc{i}(nev) = Ppeaktime; 
                pinfo.lapspatialcrr.evLapMean2ndPSlidecrr{i}(nev) = Speakcrr; pinfo.lapspatialcrr.evLapMean2ndPSlideZsc{i}(nev) = SZsc; 
                pinfo.lapspatialcrr.evLapMean2ndPSlideloc{i}(nev) = Speaktime;
                %pinfo.lapspatialcrr.evLapMean3rdPSlidecrr{i}(nev) = Tpeakcrr; pinfo.lapspatialcrr.evLapMean3rdPSlideZsc{i}(nev) = TZsc; 
                %pinfo.lapspatialcrr.evLapMean3rdPSlideloc{i}(nev) = Tpeaktime;
            end
        end
    end
end
end
%end

function [rr, intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime, Tpeakcrr, TZsc, Tpeaktime] = findcrrparameters(crr, T, binsize) %, smoothmode, smoothbin)
ncr = numel(crr);
simi = [];
for (i = 1:ncr)
    for (j = i+1:ncr)
        r1 = crr{i}; r2 = crr{j}; iii = find( (~isnan(r1)) & (~isnan(r2)) ); 
        if (~isempty(iii))
           r1 = r1(iii); r2 = r2(iii); cr = computerrcorr(r1, r2); simi= [simi cr]; 
        end
    end
end
rr = NaN; nn = numel(simi); if (nn > 0) rr = mean(simi(~isnan(simi))); end
%%%%%%%%%%%%%%also plot the averaged trace
if (ncr == 1)
    mmcc = crr{1};
else
   crnow = zeros(ncr, numel(T)); mmcc = NaN*ones(1, numel(T)); %zeros(1, numel(T));
   for (i = 1:ncr)
        crnow(i,:) = crr{i};
   end
   for (j = 1:numel(T))
       rrr = crnow(:,j); rrr = rrr(~isnan(rrr));
       if (~isempty(rrr)) mmcc(j) = mean(rrr); end
   end
end
%hf = figure; plot(T, mmcc);
intcrr = NaN; Ppeakcrr = NaN; Ppeaktime = NaN; Speakcrr = NaN; Speaktime = NaN; PZsc = NaN; SZsc = NaN; Tpeakcrr=NaN; TZsc=NaN; Tpeaktime = NaN;
mmlag = max(T); bintime = 5*binsize; searchrange = mmlag; %searchbin = round(searchrange/binsize); 
ind = find( (T>=-searchrange) & (T <= searchrange) ); %mmccnow = mmcc(ind); %TTnow = T(ind);
%[mV, iit] = max(mmccnow); mT = TTnow(iit); 
%for (i = 1:ncr)
    [intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime, Tpeakcrr, TZsc, Tpeaktime] = ...
        determinepeaks(mmcc, T/binsize, [-bintime bintime], binsize, [-searchrange searchrange], [-searchrange searchrange], 'peak');
%end
%disp('************');
%disp([rr intcrr Ppeakcrr PZsc Ppeaktime Speakcrr SZsc Speaktime]);
% allV = [Ppeakcrr Speakcrr mV]; allT = [Ppeaktime Speaktime mT];
% ii = find(allV ~= mV); allV = allV(ii); allT = allT(ii); [sV, ii] = max(allV); sT = allT(ii);
% Ppeakcrr = mV; Ppeaktime = mT; Speakcrr = sV; Speaktime = sT; 
% Zmm = zscore(mmccnow); PZsc = max(Zmm); 
% ii = find(mmccnow==sV); if ~isempty(ii) SZsc =Zmm(ii(1)); end 
% disp([[rr intcrr Ppeakcrr PZsc Ppeaktime Speakcrr SZsc Speaktime]]);
%disp('************');

function [intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime, Tpeakcrr, TZsc, Tpeaktime] = determinepeaks(crr, lagbin, intbin, binsize, P1searchwin, P2searchwin, searchmode)
Ppeakcrr = NaN; Ppeaktime = NaN; Speakcrr = NaN; Speaktime = NaN; PZsc = NaN; SZsc = NaN; 
Tpeakcrr=NaN; TZsc=NaN; Tpeaktime = NaN; zcrr = zscore(crr);
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
         [mmx, iii] = max(crr(all1stind));
      elseif (strncmpi(searchmode, 'trough', 2))
         [mmx, iii] = min(crr(all1stind));
      else
         [mmx, iii] = max(abs(crr(all1stind)));
      end
      Ppeakcrr = crr(all1stind(iii)); Ppeaktime = binsize*lagbin(all1stind(iii)); PZsc = zcrr(all1stind(iii));
      exind1 = find(crr(allind) == mmx);
   end
   %%%secondary peak
   all2ndind = intersect(setdiff(allind, allind(exind1)), P2win);
   if (~isempty(all2ndind))
      if (strncmpi(searchmode, 'peak', 2))
         [smx, jjj] = max(crr(all2ndind));
      elseif (strncmpi(searchmode, 'trough', 2))
         [smx, jjj] = min(crr(all2ndind));
      else
         [smx, jjj] = max(abs(crr(all2ndind)));
      end
      Speakcrr = crr(all2ndind(jjj)); Speaktime = binsize*lagbin(all2ndind(jjj)); SZsc = zcrr(all2ndind(jjj));
      exind2 = find(crr(allind) == smx);
   end
%    %%%Tertiary peak
%    allgone = union(allind(exind1), allind(exind2));
%    all3rdind = intersect(setdiff(allind, allgone), P2win);
%    if (~isempty(all3rdind))
%       if (strncmpi(searchmode, 'peak', 2))
%          [~, iii] = max(crr(all3rdind));
%       elseif (strncmpi(searchmode, 'trough', 2))
%          [~, iii] = min(crr(all3rdind));
%       else
%          [~, iii] = max(abs(crr(all3rdind)));
%       end
%       Tpeakcrr = crr(all3rdind(iii)); Tpeaktime = binsize*lagbin(all3rdind(iii)); TZsc = zcrr(all3rdind(iii));
%    end
end

function newcrr = smoothcrr(crr, smoothbin)
newcrr = zeros(size(crr)); %NaN*ones(size(crr)); 
halfbin = floor(smoothbin/2);
for (i = 1:numel(crr))
    kkk = [i-halfbin:i+halfbin]; kkk = kkk( (kkk>=1) & (kkk<=numel(crr)) );
    crrnow = crr(kkk); newcrr(i) = mean(crrnow(~isnan(crrnow)));
end

function [minindex, maxindex] = FindLocal(dat)
%%NOW IT WORKS REALLY WELL!!!! -- Be careful: data need to be continuous and non-NaN
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
      ccc = NaN; %0;
else
      ccc = (sumsq12 - sum1*sum2/nbin) / sqrt(co1*co2);
end

function [crr, T] = computecrrhere(spike1, spike2, maxlag, binsize)
np = ceil(maxlag/binsize); lagbin = (-np:1:np); %%%%lagbin is a row vector
T = lagbin*binsize; crr = cell(numel(spike1),1);
%[m1, n1] = size(spike1); if (m1>1) spike1 = spike1'; end
%[m1, n1] = size(spike2); if (m1>1) spike2 = spike2'; end %%%rates are now all row vectors
for (i = 1:numel(spike1))
    crr{i} = computelapcrr(spike1{i}, spike2{i}, lagbin);
end
function ccc = computelapcrr(cn1, cn2, lagbin) 
nc = numel(lagbin); nbinnow = numel(cn1); ccc = NaN*ones(size(lagbin)); %ccc = zeros*ones(size(lagbin)); 
for (i= 1:nc)
    %%%%using matrix instead of cell saves a lot of time
    startbin1 = max([1 1-lagbin(i)]);
    endbin1 = min([nbinnow nbinnow-lagbin(i)]);
    startbin2 = max([1 1+lagbin(i)]);
    endbin2 = min([nbinnow nbinnow+lagbin(i)]);
    if ( (startbin1<endbin1) & (startbin2<endbin2) ) %if both series have space for the lag in this episode
        rr1 = cn1(startbin1:endbin1); rr2 = cn2(startbin2:endbin2);
        iii = find( (~isnan(rr1)) & (~isnan(rr2)) ); nbin = numel(iii);
        if (~isempty(iii))
            rr1 = rr1(iii); rr2 = rr2(iii);
            sum1 = sum(rr1); sum2 = sum(rr2);
            sumsq1 = sum(rr1.^2); sumsq2 = sum(rr2.^2);
            sumsq12 = rr1 * rr2'; %%%row vector times column vector
            co1 = sumsq1-sum1^2/nbin; co2 = sumsq2-sum2^2/nbin;
            if ( (co1 == 0) | (co2 == 0) )
               ccc(i) = NaN; %0;
            else
               ccc(i) = (sumsq12 - sum1*sum2/nbin) / sqrt(co1*co2);
            end
        end
    end
end

function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end   

function spikerate = circleslideratecurve(spikerate)
nev = numel(spikerate);
for (i = 1:nev)
     nbin = numel(spikerate{i}); IX = mod((1:nbin)+ ceil(rand*nbin), nbin) + 1; spikerate{i} = spikerate{i}(IX);
end

function [D1rate, nact] = getlinearmap(spiketime, evTime, lappostime, lapx, xbin, occutime, smParm)
nx = numel(xbin); D1rate = zeros(nx,1); nlap = numel(lappostime); count = zeros(nx,1); spikex = []; nact = NaN;
if (nx > 1)
    binsize = xbin(2) - xbin(1); sigma = round(smParm.d1sigma/binsize); %%sigma now in number of bins
%%%%%count spikes in xbin lap bu lap
for (i = 1:nlap)
    counti = zeros(nx,1);
    spikenow = sort( spiketime( (spiketime>=evTime.start(i)) & (spiketime<=evTime.ent(i)) ) );
    spikex = NaN*ones(size(spikenow)); [laptimenow, iii] = sort(lappostime{i}); allxnow = lapx{i}(iii);
    lastpoint = 1;  %this only for saving time
    for (j = 1:numel(spikenow))
         for (k = lastpoint:numel(laptimenow)-1) %find corresponding time in position data
             if (laptimenow(k) <= spikenow(j)) && (laptimenow(k+1) > spikenow(j)) 
                 spikex(j) = allxnow(k); lastpoint = k; break; 
             end
         end
    end
    for (j = 1:nx-1) counti(j) = numel(find((spikex>=xbin(j)) & (spikex<xbin(j+1)))); end
    counti(nx) = numel(find(spikex>=xbin(nx)));
    count = count + counti;
end
%%%%%compute firing rate
for (i = 1:nx)
     if (occutime(i) ~= 0) D1rate(i) = count(i)/occutime(i); end
end
nact = numel(find(count>0))/numel(find(occutime>0));
%%now do a 1d smoothing
if (strcmpi(smParm.d1sm, 'yes'))
    %D1rate = OneDSmooth_new(D1rate, occutime, sigma, smParm.d1Nsig);
    XX = [-smParm.d1Nsig*sigma:1:smParm.d1Nsig*sigma]; 
    wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
    D1rate = conv(D1rate, wind); D1rate = D1rate((nw-1)/2+1:numel(D1rate)-(nw-1)/2);
end
end