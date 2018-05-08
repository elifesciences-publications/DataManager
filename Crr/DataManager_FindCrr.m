function [pinfo, data] = DataManager_FindCrr(pinfo, data, cellind, vv)
%%compute pair-wise auto-/cross-correlations: cellind here is the pairind
%%fields assigned here:
if (~isempty(cellind))
  nspike = numel(pinfo.general.animalname); %%%here is the total number of pairs
  if (~isfield(pinfo, 'crr')) pinfo.crr = []; end
  if (~isfield(pinfo.crr, 'sessIntcrr')) pinfo.crr.sessIntcrr = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'sessCrrNorm')) pinfo.crr.sessCrrNorm = cell(1, nspike); end   
  if (~isfield(pinfo.crr, 'sess1stPcrr')) pinfo.crr.sess1stPcrr = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'sess1stPZsc')) pinfo.crr.sess1stPZsc = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'sess1stPtime')) pinfo.crr.sess1stPtime = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'sess2ndPcrr')) pinfo.crr.sess2ndPcrr = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'sess2ndPZsc')) pinfo.crr.sess2ndPZsc = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'sess2ndPtime')) pinfo.crr.sess2ndPtime = cell(1, nspike); end 
  %if (~isfield(pinfo.crr, 'sess3rdPcrr')) pinfo.crr.sess3rdPcrr = cell(1, nspike); end 
  %if (~isfield(pinfo.crr, 'sess3rdPZsc')) pinfo.crr.sess3rdPZsc = cell(1, nspike); end 
  %if (~isfield(pinfo.crr, 'sess3rdPtime')) pinfo.crr.sess3rdPtime = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'evIntcrr')) pinfo.crr.evIntcrr = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'evCrrNorm')) pinfo.crr.evCrrNorm = cell(1, nspike); end  
  if (~isfield(pinfo.crr, 'ev1stPcrr')) pinfo.crr.ev1stPcrr = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'ev1stPZsc')) pinfo.crr.ev1stPZsc = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'ev1stPtime')) pinfo.crr.ev1stPtime = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'ev2ndPcrr')) pinfo.crr.ev2ndPcrr = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'ev2ndPZsc')) pinfo.crr.ev2ndPZsc = cell(1, nspike); end 
  if (~isfield(pinfo.crr, 'ev2ndPtime')) pinfo.crr.ev2ndPtime = cell(1, nspike); end 
  %if (~isfield(pinfo.crr, 'ev3rdPcrr')) pinfo.crr.ev3rdPcrr = cell(1, nspike); end 
  %if (~isfield(pinfo.crr, 'ev3rdPZsc')) pinfo.crr.ev3rdPZsc = cell(1, nspike); end 
  %if (~isfield(pinfo.crr, 'ev3rdPtime')) pinfo.crr.ev3rdPtime = cell(1, nspike); end 
  if (~isfield(data, 'crr')) data.crr = []; end
  if (~isfield(data.crr, 'sesscrr')) data.crr.sesscrr = cell(1, nspike); end
  if (~isfield(data.crr, 'evcrr')) data.crr.evcrr = cell(1, nspike); end
  if (~isfield(data.crr, 'timebin')) data.crr.timebin = cell(1, nspike); end
end
sss = [num2str(numel(cellind)) ' pairs found. Compute crr?']; 
aa = questdlg(sss);
if (strcmp(aa, 'Yes'))
for (jjjjk = 1:numel(cellind))
    i = cellind(jjjjk);
    disp(['-----------> compute crr (', num2str(jjjjk), ' out of ', num2str(numel(cellind)), ') -- ', pinfo.general.clname{i}]);
    %%%%%%get all parameters
    crrmode = pinfo.parm.crrmode{i}; binsize = pinfo.parm.timebin(i); maxlag = pinfo.parm.maxlag(i); setbacktime = pinfo.parm.setbacktime(i);
    timeunit = pinfo.parm.timeunit(i); pairind = data.crr.cellind{i};
    intbin = pinfo.parm.intbin{i}; crrtype = pinfo.general.crrtype{i};
    smoothmode = pinfo.parm.smoothmode{i}; smoothbin = pinfo.parm.smoothbin(i);
    P1searchwin = pinfo.parm.P1stSearchWin{i}; P2searchwin = pinfo.parm.P2ndSearchWin{i}; 
    searchmode = pinfo.parm.searchPmode{i}; %%%%%%either 'peack', 'trough', or 'both'
    sesskeyword = pinfo.parm.sesskeyword{i}; %  = 'track';
    sesskeytype = pinfo.parm.sesskeytype{i}; % = 'linear';
    evtkeyword = pinfo.parm.evtkeyword{i}; %  = 'sleep';
    evtkeytype = pinfo.parm.evtkeytype{i}; % = 'sws';
    %%%%rectify spikes all into row vectors %%%%using default column vectors does not work well, since matlab treat scalar as row vectors
    spike1 = data.spike.spiketime{pairind(1)}*timeunit; 
    [A, B] = size(spike1); if (B == 1) spike1 = spike1'; end
    spike2 = data.spike.spiketime{pairind(2)}*timeunit;
    [A, B] = size(spike2); if (B == 1) spike2 = spike2'; end
    %%%%set up time bins
    np = ceil(maxlag/binsize); lagbin = (-np:1:np); %%%%lagbin is a row vector
    timebinnow = lagbin*binsize; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The following is for within ripple cross-correlaton normalizaton to
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ripple baseline [-0.4 -018] and
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [0.18 0.4]; theta baseline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [0.24 0.36]s
    normbin = find( ((timebinnow>-0.36001)&(timebinnow < -0.23999)) | ((timebinnow>0.239999)&(timebinnow < 0.36001)) );
    if isempty(normbin) && strcmp(crrmode, 'normbase')
        disp(['---------------------> warning: no normalization baseline']);
    end
    data.crr.timebin{i} = timebinnow;
    %%%%%session crr   
    sessions = pinfo.general.sessionname{i}; nsess = numel(sessions);
    sesstype = pinfo.parm.sessType{i};
    for (tt = 1:nsess)
      if isselected(sessions{tt}, sesstype{tt}, sesskeyword, sesskeytype)
          %disp(sessions{tt});
        sT = pinfo.general.sessionstartT{i}(tt); eT = pinfo.general.sessionendT{i}(tt); 
        [crr, norm] = computecrrhere(spike1, spike2, sT, eT, lagbin, crrmode, binsize, normbin); %%%spikes here are column vectors
        if (strncmpi(smoothmode, 'yes', 1))
            crr = smoothcrr(crr, smoothbin);
%             if strncmpi(crrtype, 'auto', 2) 
%                 if strncmpi(crrmode, 'rate', 2)
%                    crr(lagbin == 0) = 1;
%                 %elseif strncmpi(crrmode, 'count', 2) %%%%do not forget to change the event crr as well
%                 %   crr(lagbin == 0) = numel(spike1);
%                 %elseif strncmpi(crrmode, 'normcount', 2)
%                 %   if (norm>0) crr(lagbin == 0) = (numel(spike1)-norm)/sqrt(norm); end
%                 end
%             end
        end
        data.crr.sesscrr{i}{tt} = crr; 
        [intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime, Tpeakcrr, TZsc, Tpeaktime] = determinepeaks(crr, lagbin, intbin, binsize, P1searchwin, P2searchwin, searchmode);
        pinfo.crr.sessIntcrr{i}{tt} = intcrr; pinfo.crr.sessCrrNorm{i}{tt} = norm;
        pinfo.crr.sess1stPcrr{i}{tt} = Ppeakcrr; pinfo.crr.sess1stPZsc{i}{tt} = PZsc; pinfo.crr.sess1stPtime{i}{tt} = Ppeaktime; 
        pinfo.crr.sess2ndPcrr{i}{tt} = Speakcrr; pinfo.crr.sess2ndPZsc{i}{tt} = SZsc; pinfo.crr.sess2ndPtime{i}{tt} = Speaktime; 
        %pinfo.crr.sess3rdPcrr{i}{tt} = Tpeakcrr; pinfo.crr.sess3rdPZsc{i}{tt} = TZsc; pinfo.crr.sess3rdPtime{i}{tt} = Tpeaktime; 
      else
        data.crr.sesscrr{i}{tt} = []; 
        pinfo.crr.sessIntcrr{i}{tt} = []; pinfo.crr.sessCrrNorm{i}{tt} = [];
        pinfo.crr.sess1stPcrr{i}{tt} = []; pinfo.crr.sess1stPZsc{i}{tt} = []; pinfo.crr.sess1stPtime{i}{tt} = []; 
        pinfo.crr.sess2ndPcrr{i}{tt} = []; pinfo.crr.sess2ndPZsc{i}{tt} = []; pinfo.crr.sess2ndPtime{i}{tt} = [];  
      end
    end
    %%%%event crr 
    evTime = data.events.eventtimes{pairind(1)}; evType = pinfo.parm.eventtype{i}; evName = pinfo.general.eventname{i};
    for (j = 1:numel(evTime))
      if isselected(evName{j}, evType{j}, evtkeyword, evtkeytype)  
          %disp(evName{j});
        sT = evTime{j}.start-setbacktime; eT = evTime{j}.ent+setbacktime;
        [crr, norm] = computecrrhere(spike1, spike2, sT, eT, lagbin, crrmode, binsize, normbin);
        if (strncmpi(smoothmode, 'yes', 1))
            crr = smoothcrr(crr, smoothbin);
%             if (strncmpi(crrtype, 'auto', 2')) && (strncmpi(crrmode, 'rate', 2))
%                 crr(lagbin == 0) = 1;
%             end
        end
        data.crr.evcrr{i}{j} = crr; 
        [intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime, Tpeakcrr, TZsc, Tpeaktime] = determinepeaks(crr, lagbin, intbin, binsize, P1searchwin, P2searchwin, searchmode);
        pinfo.crr.evIntcrr{i}{j} = intcrr; pinfo.crr.evCrrNorm{i}{j} = norm;
        pinfo.crr.ev1stPcrr{i}{j} = Ppeakcrr; pinfo.crr.ev1stPZsc{i}{j} = PZsc; pinfo.crr.ev1stPtime{i}{j} = Ppeaktime; 
        pinfo.crr.ev2ndPcrr{i}{j} = Speakcrr; pinfo.crr.ev2ndPZsc{i}{j} = SZsc; pinfo.crr.ev2ndPtime{i}{j} = Speaktime; 
        %pinfo.crr.ev3rdPcrr{i}{j} = Tpeakcrr; pinfo.crr.ev3rdPZsc{i}{j} = TZsc; pinfo.crr.ev3rdPtime{i}{j} = Tpeaktime;
      else
        data.crr.evcrr{i}{j} = []; 
        pinfo.crr.evIntcrr{i}{j} = []; pinfo.crr.evCrrNorm{i}{j} = [];
        pinfo.crr.ev1stPcrr{i}{j} = []; pinfo.crr.ev1stPZsc{i}{j} = []; pinfo.crr.ev1stPtime{i}{j} = []; 
        pinfo.crr.ev2ndPcrr{i}{j} = []; pinfo.crr.ev2ndPZsc{i}{j} = []; pinfo.crr.ev2ndPtime{i}{j} = []; 
      end
    end
end
end

function [intcrr, Ppeakcrr, PZsc, Ppeaktime, Speakcrr, SZsc, Speaktime, Tpeakcrr, TZsc, Tpeaktime] = determinepeaks(crr, lagbin, intbin, binsize, P1searchwin, P2searchwin, searchmode)
Ppeakcrr = NaN; Ppeaktime = NaN; Speakcrr = NaN; Speaktime = NaN; PZsc = NaN; SZsc = NaN; 
Tpeakcrr = NaN; TZsc = NaN; Tpeaktime = NaN; zcrr = zscore(crr);
%%%%for intcrr
intbin = intbin/binsize; intwin = find( (lagbin>=intbin(1)) & (lagbin<=intbin(2)) );
centercrr = crr(intersect(intwin, 1:numel(crr)));
intcrr = mean(centercrr(~isnan(centercrr))); 
%%%%find primary/secondary peaks
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
   all1stind = intersect(allind, P1win); iii = []; exind1 = [];
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

function [crr, norm] = computecrrhere(spike1, spike2, sT, eT, lagbin, crrmode, binsize, normbin)
crr = zeros(size(lagbin)); %NaN*ones(size(lagbin)); 
norm = 0;
if (~isempty(spike1)) && (~isempty(spike2))
    iii = find(eT-sT>=binsize); sT = sT(iii); eT = eT(iii); %%%filter out short events
    if (strncmpi(crrmode, 'rate', 2))
        crr = computeratecrr(spike1, spike2, sT, eT, lagbin, binsize);
    elseif (strncmpi(crrmode, 'count', 2))
        crr = computecountcrr(spike1, spike2, sT, eT, lagbin, binsize);
    elseif (strncmpi(crrmode, 'normcount', 2))
        [crr, norm] = computenormcountcrr(spike1, spike2, sT, eT, lagbin, binsize);
    elseif (strncmpi(crrmode, 'basenorm', 2))
        [crr, norm] = computebasenormcountcrr(spike1, spike2, sT, eT, lagbin, binsize, normbin);    
    end
end

function ccc = computeratecrr(t1, t2, sT, eT, lagbin, binsize) %%%t1, t2 are row vectors
nev = numel(sT); timebin = cell(1,nev); nc = numel(lagbin); ccc = zeros(size(lagbin)); %NaN*ones(size(lagbin)); 
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
           ccc(i) = 0;
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
%%%%         mm = numel(spike1)*numel(spike2)*bin/totallengthofspiketrains;
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
%disp(num2str([numel(t1) numel(t2) tleng mm]));
if (mm>0) count = (count-mm)/sqrt(mm); end
end

function [count, mm] = computebasenormcountcrr(spike1, spike2, sT, eT, lagbin, bin, normbin)
%%%%normalized count: assuming spike1 and spike2 are independent poisson, at each time lag bin, the mean expected spike count is
%%%%         mm = numel(spike1)*numel(spike2)*bin/totallengthofspiketrains;
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
if (~isempty(t2))
for (i = 1:numel(t1)) % for each spike in first train, count spike number in second train that within a time window defined by bin
      k = -np:1:np; wind = t1(i) + k*bin - bin/2;
      count(k+np+1) = count(k+np+1) + histc(t2, wind);
end
mm = numel(t1)*numel(t2)*bin/tleng;
if (mm>0) 
    count = count - mm;
    %%%select [-0.8 -0.6] and [0.6 0.8] as baseline for now just for within
    %%%ripple cross-correlation
    if ~isempty(normbin)
       count = count - mean(count(normbin)); % median(count(normbin)); 
    end
    count = count/sqrt(mm);
end
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

function ok = isselected(sessname, sesstype, sesskeyword, sesskeytype)
ok = 1;
if isempty(sesskeyword) && isempty(sesskeytype) ok = 0; end
if ~isempty(sesskeyword)
   if isempty(strfind(lower(sessname), lower(sesskeyword))) ok = 0; end
end
if ~isempty(sesskeytype)
   if ~strcmpi(sesstype, sesskeytype) ok = 0; end
end


