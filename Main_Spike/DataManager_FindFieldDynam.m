function [pinfo,data] = DataManager_FindFieldDynam(pinfo,data,behav,bhdata,cellind,vv)
%Compute the dynamics of spatial firing profile:
%     1. For open field, break a session into segments, compute rate maps for each segment, then compute the spatial correlation
%     2. For linear track, compute rate maps for each lap, then compute the spatial correlation
%     3. For linear track, also compute the lap-by-lap shifts of all place fields' peak/mean rate, peak/mean position, and field size 
%                          this is done only for the most prominent field (using PF1DInMeanRate) per cell

%require the following parameters (not all listed)
%   behav.parm.s1Dbinsize, behav.parm.s2Dbinsize, behav.parm.sessPmarker,
%   pinfo.parm.fd1DBaseEvt, pinfo.parm.fd1DBaseLap, pinfo.parm.fd2DBaseSess, pinfo.parm.fd2DBaseSeg; %%%pinfo.parm.MinLapSpeed, pinfo.parm.MinFieldSpeed
%   pinfo.parm.sessType, pinfo.parm.eventtype

%require the following variables (not all listed)
%   pinfo.general.sessionname, pinfo.general.eventname, pinfo.field.sessPairs; pinfo.field.run1DEvtPairs; pinfo.field.run1DTrajPairs

%requre the following data (not all listed)
%   data.spike.spiketime; data.event.eventtimes
%   bhdata.pos.postimestamp; bhdata.pos.XX/YY
%   bhdata.sess.gridXYbin; bhdata.sess.gridOccuptime; bhdata.sess.gridSegOccuptime
%   
%   bhdata.event.Xbin, bhdata.event.Xjoint, bhdata.event.Occuptime,
%   bhdata.event.LapAllPostimestamp, bhdata.event.LapAllX

%variable to assign
if (~isempty(cellind))
   nspike = numel(pinfo.general.parmfile);
   if (~isfield(pinfo, 'fielddynam')) pinfo.fielddynam = []; end
   %%%session field dynam
   if (~isfield(pinfo.fielddynam, 'sessRefSegNum')) pinfo.fielddynam.sessRefSegNum = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'sessName')) pinfo.fielddynam.sessName = cell(1, nspike); end 
   if (~isfield(pinfo.fielddynam, 'sessSegNum')) pinfo.fielddynam.sessSegNum = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'sessWithinSegNum')) pinfo.fielddynam.sessWithinSegNum = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'sessSegSparsity')) pinfo.fielddynam.sessSegSparsity = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'sessSegSpatialInfo')) pinfo.fielddynam.sessSegSpatialInfo = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'sessSegSpatialInfotime')) pinfo.fielddynam.sessSegSpatialInfotime = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'sessSegActarea')) pinfo.fielddynam.sessSegActarea = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'sessSegCrr')) pinfo.fielddynam.sessSegCrr = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'sessSegPeakrate')) pinfo.fielddynam.sessSegPeakrate = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'sessSegMeanrate')) pinfo.fielddynam.sessSegMeanrate = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'sessSegPeakshift')) pinfo.fielddynam.sessSegPeakshift = cell(1, nspike); end
   %%%evt field dynam
   if (~isfield(pinfo.fielddynam, 'runRefLapNum')) pinfo.fielddynam.runRefLapNum = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'run1DEvtname')) pinfo.fielddynam.run1DEvtname = cell(1, nspike); end 
   if (~isfield(pinfo.fielddynam, 'runLapNum')) pinfo.fielddynam.runLapNum = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runWithinLapNum')) pinfo.fielddynam.runWithinLapNum = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapSparsity')) pinfo.fielddynam.runLapSparsity = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapSpatialInfo')) pinfo.fielddynam.runLapSpatialInfo = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapSpatialInfotime')) pinfo.fielddynam.runLapSpatialInfotime = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapActarea')) pinfo.fielddynam.runLapActarea = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapCrr')) pinfo.fielddynam.runLapCrr = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapComX')) pinfo.fielddynam.runLapComX = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapPeakX')) pinfo.fielddynam.runLapPeakX = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapMeanRate')) pinfo.fielddynam.runLapMeanRate = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapPeakRate')) pinfo.fielddynam.runLapPeakRate = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapFieldSize')) pinfo.fielddynam.runLapFieldSize = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapSkew')) pinfo.fielddynam.runLapSkew = cell(1, nspike); end
   if (~isfield(pinfo.fielddynam, 'runLapKurt')) pinfo.fielddynam.runLapKurt = cell(1, nspike); end
   %%%data variables to assign
   if (~isfield(data, 'fielddynam')) data.fielddynam = []; end
   if (~isfield(data.fielddynam, 'sess2DsegRateMaps')) data.field.sess2DsegRateMaps = cell(1, nspike); end
   if (~isfield(data.fielddynam, 'evt1DlapRateMaps')) data.field.evt1DlapRateMaps = cell(1, nspike); end
end

for (jjjk = 1:numel(cellind))
    i = cellind(jjjk); 
    %%%%load computing parameters
    fd2DBaseSess = pinfo.parm.fd2DBaseSess{i};  fd2DBaseSeg = pinfo.parm.fd2DBaseSeg{i}; 
    fd1DBaseSess = pinfo.parm.fd1DBaseSess{i}; fd1DBaseLap = pinfo.parm.fd1DBaseLap{i}; 
    smParm.d2sigma = pinfo.parm.fSmooth2DSigma(i); smParm.d2Nsig = pinfo.parm.fSmooth2DNSigma(i); smParm.d2sm = pinfo.parm.fSmooth2D{i};
    smParm.d1sigma = pinfo.parm.fSmooth1DSigma(i); smParm.d1Nsig = pinfo.parm.fSmooth1DNSigma(i); smParm.d1sm = pinfo.parm.fSmooth1D{i};
    threshold1Drate = pinfo.parm.f1DThresholdRate(i); max1Dgap = pinfo.parm.f1DMaxGap(i); min1Dpeakrate = pinfo.parm.f1DMinPeakRate(i); 
    timeunit = pinfo.parm.timeunit(i); base1DrateN = pinfo.parm.f1DBaseRateGridN(i); 
    threshold2Drate = pinfo.parm.f2DThresholdRate(i); max2Dgap = pinfo.parm.f2DMaxGap(i); min2Dpeakrate = pinfo.parm.f2DMinPeakRate(i); 
    base2DrateN = pinfo.parm.f2DBaseRateGridN(i); 
    disp(strcat('-----> compute field dynamics ---', pinfo.general.parmfile{i}));
    spiketime = timeunit*data.spike.spiketime{i};
    %%%%%compute segment rate maps for each (open field) session: assuming all open sessions use the same maze at the same location
    fprate = pinfo.field.PF2DInPeakrate{i};
    if (~isempty(fprate))
    sessions = pinfo.general.sessionname{i}; nsess = numel(sessions); sesstype = pinfo.parm.sessType{i}; 
    XYbin = []; XYsess = []; D2ratenow = []; gridoccup = []; n2Dactnow = [];
    segSessName = []; sessNum = []; segNum = []; withinsegNum = []; nses = 0; nseg = 0;%%%this the final output maps
    for (tt = 1:nsess)
        if (strcmp(sesstype{tt}, 'open'))
            nses = nses + 1; XYbin{nses} = []; XYsess{nses} = []; 
            stime = pinfo.general.sessionstartT{i}(tt); etime = pinfo.general.sessionendT{i}(tt);
            %%%locate session position data
            finaldirnow = pinfo.general.finaldir{i}; sessnow = sessions{tt}; 
            posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, sessnow) );
            if (numel(posid) ~= 1)
               disp(['------------------> 2D field dynamics not computed for this session: no or more than 1 positon files match the session: ', finaldirnow, '___', sessnow]);
            else
               postime = bhdata.pos.postimestamp{posid}*behav.parm.timeunit(posid);
               XYbin{nses} = bhdata.sess.gridXYbin{posid}; XYsess{nses} = sessnow; XYsessNum(nses) = nses;
               segtime = behav.parm.sessSegTime(posid); nsegnow = floor((etime-stime)/segtime);
               %%%find which color to use
               cid = find( strcmp(behav.general.posMarker{posid}, behav.parm.sessPmarker{posid}) );        
               XX = bhdata.pos.XX{posid}{cid}; YY = bhdata.pos.YY{posid}{cid}; 
               for (j = 1:nsegnow)
                    nseg = nseg + 1;
                    stnow = stime + (j-1)*segtime; etnow = stnow + segtime; if (etnow>etime) etnow = etime; end
                    ep.start = stnow; ep.ent = etnow; [spiketimenow, epid] = SpikeEventFilter(spiketime, ep);
                    gridoccup{nseg} = bhdata.sess.gridSegOccuptime{posid}{j}; 
                    [D2map, n2Dact] = get2Dratemaps(spiketimenow, XYbin{nses}, gridoccup{nseg}, postime, XX, YY, smParm);
                    D2ratenow{nseg} = D2map; n2Dactnow(nseg) = n2Dact;
                    segSessName{nseg} = sessnow; sessNum(nseg) = nses; segNum(nseg) = nseg; withinsegNum(nseg) = j;
               end
            end
        end
    end
    data.fielddynam.sess2DRateMaps{i} = D2ratenow; 
    pinfo.fielddynam.sessName{i} = segSessName; pinfo.fielddynam.sessSegNum{i} = segNum; pinfo.fielddynam.sessWithinSegNum{i} = withinsegNum;
    
    %%%%%%%%%%%%%%%%%%%compute seg-by-seg: spatial info, sparsity, active area, correlation (stability)
    [D2refmap, XYrefbin, refsegNum] = resolve2Drefmap(segSessName, segNum, withinsegNum, D2ratenow, XYsess, XYbin, fd2DBaseSess, fd2DBaseSeg);
    spinfo = NaN*ones(1,nseg); spinfotime = NaN*ones(1,nseg); sparsity = NaN*ones(1,nseg); crr = NaN*ones(1,nseg);
    for (j = 1:nseg)
            [spinfo(j), spinfotime(j), sparsity(j)] = findratemapprop(D2ratenow{j}, gridoccup{j});
            if (~isempty(XYbin{sessNum(j)})) && (~isempty(XYrefbin))
               [crr(j), pp0] = correlateratemaps(D2ratenow{j}, D2refmap, XYbin{sessNum(j)}, XYrefbin, [0 0], 0);
            end
    end
    pinfo.fielddynam.sessSegSparsity{i} = sparsity; pinfo.fielddynam.sessSegSpatialInfo{i} = spinfo;
    pinfo.fielddynam.sessSegSpatialInfotime{i} = spinfotime;
    pinfo.fielddynam.sessSegActarea{i} = n2Dactnow;
    pinfo.fielddynam.sessSegCrr{i} = crr; pinfo.fielddynam.sessRefSegNum{i} = refsegNum;
    %%%%%%%%%%%%%%%%%%compute seg-by-seg 2D field properties: mean/peak rates and position shift from the ref
    %%%%%%%%%%%%%%%%%%%%%reference field properties
    refpeakrate = NaN; refmeanrate = NaN; refX = NaN; refY = NaN;
    [ny, nx] = size(D2refmap);
    if (nx > 1) && (ny > 1)
        baseratenow = computebaserate(D2refmap, base2DrateN);
        pfnow = find2Dfieldprop(D2refmap, XYrefbin, threshold2Drate, max2Dgap, min2Dpeakrate, baseratenow);
        prate = [];
        for (tj = 1:numel(pfnow)) prate(tj) = pfnow(tj).inpeakrate; end
        if (~isempty(prate))
            [refpeakrate, iii] = max(prate); refmeanrate = pfnow(iii).inmeanrate;
            refX = pfnow(iii).peakX; refY = pfnow(iii).peakY;
        end
    end
    prate = NaN*ones(1,nseg); mrate = NaN*ones(1,nseg); pshift = NaN*ones(1,nseg);
    for (j = 1:nseg)
        [ny, nx] = size(D2ratenow{j});
        if (nx > 1) && (ny > 1)
            baseratenow = computebaserate(D2ratenow{j}, base2DrateN);
            pf = find2Dfieldprop(D2ratenow{j}, XYbin{sessNum(j)}, threshold2Drate, max2Dgap, min2Dpeakrate, baseratenow);
            prateok = [];
            for (tj = 1:numel(pf)) prateok(tj) = pf(tj).inpeakrate; end
            if (~isempty(prateok))
               [prate(j), iii] = max(prateok); prate(j) = prate(j)-refpeakrate; mrate(j) = pf(iii).inmeanrate-refmeanrate;
               pshift(j) = sqrt( (refX-pf(iii).peakX)^2 + (refY-pf(iii).peakY)^2 );
            end
        end
    end
    pinfo.fielddynam.sessSegPeakrate{i} = prate; pinfo.fielddynam.sessSegMeanrate{i} = mrate;
    pinfo.fielddynam.sessSegPeakshift{i} = pshift;
    else
        disp('----------------> 2D dynamics not computed: no place fields found');
    end
    %%%%%%here all 2D map processings are done
        
    %%%%%%%%compute lap-by-lap 1D rate maps on the trajectory the cell's most prominent field (maximum peak rate) resides
    fprate = pinfo.field.PF1DInPeakrate{i}; 
    if (~isempty(fprate)) %%%if fields exist on linear track
       [~,iii] = max(fprate); evNamenow = pinfo.field.PF1Devt{i}{iii}; [~,evTrajnow] = strtok(evNamenow, '_');
       evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i}; 
       finaldirnow = pinfo.general.finaldir{i};
       D1ratenow = []; xbin = []; xev = []; occutime = []; n1Dact = []; 
       lapevName = []; evNum = []; lapNum = []; withinlapNum = []; nev = 0; nlap = 0;
       evD1rate = []; evOccutime = []; %%%%this is the ev-averaged ratemap and occutime, used to re-define avgratemap field bounds across evts 
       for (j = 1:numel(evTime)) 
         if (strcmp(evType{j}, 'run')) && (~isempty(strfind(evName{j},evTrajnow))) %% && (strncmpi(evName{j}, 'run', 3)) %%%if this is a run event on the selected traj 
             disp(['----------------------------> event name ', evName{j}]);
            %%%%locate event position data
            nev = nev + 1;
            evSess{j} = identifysession(evTime{j}, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
            posid = []; evid = [];
            if (~isempty(evSess{j}))
                posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess{j}) );
            end
            if numel(posid == 1)
                if (isfield(behav.general, 'eventname'))
                    evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
                else
                    evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
                end
            end
            if (numel(posid)~=1)||(numel(evid)~=1)
               disp(['------------------> 1D field dynamics not computed for this event: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
            else
               evD1rate{nev} = data.field.evt1DRateMaps{i}{j}; evOccutime{nev} = bhdata.event.Occuptime{posid}{evid};
               xbin{nev} = bhdata.event.Xbin{posid}{evid}; xev{nev} = evName{j}; 
               nlapnow = numel(evTime{j}.start);
               for (tj = 1:nlapnow)
                   nlap = nlap + 1;
                   occutime{nlap} = bhdata.event.LapOccuptime{posid}{evid}{tj}; 
                   lappostime{1} = bhdata.event.LapAllPostimestamp{posid}{evid}{tj}; lapx{1} = bhdata.event.LapAllX{posid}{evid}{tj};
                   evok.start = evTime{j}.start(tj); evok.ent = evTime{j}.ent(tj);
                   [D1ratenow{nlap}, n1Dact(nlap)] = getlinearmap(spiketime, evok, lappostime, lapx, xbin{nev}, occutime{nlap}, smParm);
                   lapevName{nlap} = evName{j}; lapNum(nlap) = nlap; withinlapNum(nlap) = tj; evNum(nlap) = nev;
               end
            end
        end
      end
      data.fielddynam.evt1DRateMaps{i} = D1ratenow; pinfo.fielddynam.run1DEvtname{i} = lapevName; pinfo.fielddynam.runLapNum{i} = lapNum;
      pinfo.fielddynam.runWithinLapNum{i} = withinlapNum;
    
    %%%%%%%%%%%%%%%%%%%compute lap-by-lap: spatial info, sparsity, active area, correlation (stability)
    [D1refmap, xrefbin, reflapNum] = resolve1Drefmap(lapevName, lapNum, withinlapNum, D1ratenow, xev, xbin, fd1DBaseSess, fd1DBaseLap);
    spinfo = NaN*ones(1,nlap); spinfotime = NaN*ones(1,nlap); sparsity = NaN*ones(1,nlap); crr = NaN*ones(1,nlap);
    for (j = 1:nlap)
            [spinfo(j), spinfotime(j), sparsity(j)] = findratemapprop(D1ratenow{j}, occutime{j});
            if (~isempty(xbin{evNum(j)})) && (~isempty(xrefbin))
               [crr(j), pp0] = correlateratemaps(D1ratenow{j}, D1refmap, xbin{evNum(j)}, xrefbin, 0, 0);
            end
    end
    pinfo.fielddynam.runLapSparsity{i} = sparsity; pinfo.fielddynam.runLapSpatialInfo{i} = spinfo;
    pinfo.fielddynam.runLapSpatialInfotime{i} = spinfotime; pinfo.fielddynam.runRefLapNum{i} = reflapNum;
    pinfo.fielddynam.runLapActarea{i} = n1Dact; pinfo.fielddynam.runLapCrr{i} = crr;
    
%     %%%%%%%%%%%%%%The following needs to be finalized. Disable now!!!!%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%compute lap-by-lap 1D field properties: mean/peak rates and position shift from the ref
%     %%%%%%a concern here: fboundS, fboundE -these are defined on a traj during 1 sess, dynam are computed for the traj on all sess 
%     %%%%%%trajectories in different sessions need to be aligned and field boundaries need to be redefined
%     [D1ratenow, D1refmap, D1avgmap, xxbin] = alignlaps(D1ratenow, D1refmap, evD1rate, evOccutime, evNum, xbin, xrefbin); 
%     %%%%compute 1D place fields
%     nx = numel(xxbin);
%     if (nx > 1)
%         baserate = computebaserate(D1avgmap, base1DrateN);
%         pf = find1Dfieldprop(D1avgmap, xxbin, threshold1Drate, max1Dgap, min1Dpeakrate, baserate);
%         if (numel(pf)>=1)
%            ppeakrate = [];
%            for (j = 1:numel(pf)) ppeakrate(j) = pf(j).inpeakrate; end
%            [~, ii] = max(ppeakrate); fstart = pf(ii).boundstart; fend = pf(ii).boundend;
%            %[mmrateref, pprateref, fsizeref, comXref, plocref, skewref, kurtref] = find1DFFprop(D1refmap, xxbin, fstart, fend);
%            mmrate = NaN*ones(1,nlap); pprate = NaN*ones(1,nlap); kurt = NaN*ones(1, nlap);
%            fsize = NaN*ones(1,nlap); comX = NaN*ones(1,nlap); ploc = NaN*ones(1, nlap); skew = NaN*ones(1, nlap);
%            for (j = 1:nlap)
%                [mmrate(j), pprate(j), fsize(j), comX(j), ploc(j), skew(j), kurt(j)] = find1DFFprop(D1ratenow{j}, xxbin, fstart, fend);
%            end
% %            pinfo.fielddynam.runLapMeanRate{i} = mmrate-mmrateref; pinfo.fielddynam.runLapPeakRate{i} = pprate-pprateref;
% %            pinfo.fielddynam.runLapFieldSize{i} = fsize-fsizeref; pinfo.fielddynam.runLapComX{i} = comX-comXref;
% %            pinfo.fielddynam.runLapPeakX{i} = ploc-plocref; 
% %            pinfo.fielddynam.runLapSkew{i} = skew-skewref; pinfo.fielddynam.runLapKurt{i} = kurt-kurtref;
%            pinfo.fielddynam.runLapMeanRate{i} = mmrate-mean(mmrate(~isnan(mmrate))); pinfo.fielddynam.runLapPeakRate{i} = pprate-mean(pprate(~isnan(pprate)));
%            pinfo.fielddynam.runLapFieldSize{i} = fsize-mean(fsize(~isnan(fsize))); pinfo.fielddynam.runLapComX{i} = comX-mean(comX(~isnan(comX)));
%            pinfo.fielddynam.runLapPeakX{i} = ploc-mean(ploc(~isnan(ploc))); 
%            pinfo.fielddynam.runLapSkew{i} = skew-mean(skew(~isnan(skew))); pinfo.fielddynam.runLapKurt{i} = kurt-mean(kurt(~isnan(kurt)));
%         else
%             disp('----------------> dynamics of 1D field properties not computed: no average 1D field found');
%         end
%     else
%         disp('----------------> dynamics of 1D field properties not computed: no common space found');
%     end
    else
        disp('----------------> 1D dynamics not computed: no place fields found');
    end
end
%%%%%%%%%%%%%%%END OF THE MAIN ROUTINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pf = find2Dfieldprop(rateingrid, xybin, zerothreshold, maxgap, minpeakrate, baselinerate) %%%%do 1D first, 2D is harder to delinearate
pf = [];  %%%%%ATTENTION: rateingrid(y,x)
%%algorithm: similar to the 1D case, subtract from baseline firing rate (for ctx), then examine maximum peaks one by one.
%%           For each peak, accumulate the pixels surround the peak
%%input: rateingrid - spike firing rate in each grid on lineared track grid
%%       zerothreshold - a threshold (in percentage of average firing rate)
%%                    below which the area is not in its place field
%%       maxgap - maximum gap (rate below zero threshold) in a number of pixels within a place field
%%output: numfield - number of place fields
numfield = 0; binsize = xybin{1}(2)-xybin{1}(1); maxgap = round(maxgap/binsize); %%maxgap now in number of bins
gridnum = size(rateingrid); %%%%this is the grid
%%%%subtract from the baseline rate
rateingrid = rateingrid - baselinerate; iii = find(rateingrid<0); if (~isempty(iii)) rateingrid(iii) = zeros(1, numel(iii)); end
%%%%find the peak: rate in peakratenow, loc = [peakind1 peakind2]
[peakR, peakind1] = max(rateingrid); [peakratenow, peakind2] = max(peakR); peakind1 = peakind1(peakind2); 
while (peakratenow >= minpeakrate) %%%%find all the non-zero rate pixels surrounding the peak
    numfield = numfield + 1;
    ind = []; ind(1,1) = peakind1; ind(1,2) = peakind2; np = 1; %%this store all the pixels in the field
    for (i = 1:gridnum(1))
        for (j = 1:gridnum(2))
            if (rateingrid(i,j)>=peakratenow*zerothreshold) %%if the grid rate is high (rates in grids that belong to an already identified field are asigned to 0, see below
                np = np + 1; ind(np,:) = [i j];
            end
        end
    end
    %%%%%now need to check the connectivity (xgap+ygap <= 2*maxgap) in ind
    isneibor = zeros(np,1); isneibor(1) = 1; nzero = find(isneibor == 0); %%%%decide if a grid is connected to the peak ind(1,1)
    for (i = 1:numel(nzero))
        gridtocheck = nzero(i); connectedgrid = ind(isneibor == 1,:);   
        if (min(abs(connectedgrid(:,1)-ind(gridtocheck,1)) + abs(connectedgrid(:,2)-ind(gridtocheck,2))) <= 2*maxgap)
                                                                        %%%this is the key gap check for connectivity
            isneibor(gridtocheck) = 1;
        end
    end
    nzeronow = find(isneibor == 0); 
    while (~isempty(nzeronow)) && (numel(nzeronow) < numel(nzero)) %%%if still room to check
        nzero = nzeronow;
        for (i = 1:numel(nzero))
            gridtocheck = nzero(i); connectedgrid = ind(isneibor == 1,:);   
            if (min(abs(connectedgrid(:,1)-ind(gridtocheck,1)) + abs(connectedgrid(:,2)-ind(gridtocheck,2))) <= 2*maxgap)
               isneibor(gridtocheck) = 1;
            end
        end
        nzeronow = find(isneibor == 0);
    end
    ind = ind(isneibor==1,:); np = size(ind,1); ratenow = zeros(np,1); %%%%this stores all the grids in the field
    for (i = 1:np)
        ratenow(i) = rateingrid(ind(i,1), ind(i,2)); rateingrid(ind(i,1), ind(i,2)) = 0; %%zero rates for the grids in the current field
    end
    %%%%%%characteristics of the current field
    pf(numfield).peakX = xybin{1}(peakind2); pf(numfield).peakY = xybin{2}(peakind1); pf(numfield).inpeakrate = peakratenow + baselinerate;
    pf(numfield).area = np*binsize*binsize; pf(numfield).size = binsize*binsize*sum(ratenow); 
    pf(numfield).inmeanrate = baselinerate + mean(ratenow);
    %%%%%%leave leftovers for computinng the next field
    [peakR, peakind1] = max(rateingrid); [peakratenow, peakind2] = max(peakR); peakind1 = peakind1(peakind2); 
end

function pf = find1Dfieldprop(rateingrid, xbin, zerothreshold, maxgap, minpeakrate, baselinerate)
pf = [];
%%algorithm: ubtract from baseline firing rate (for ctx), then examine maximum peak layer by layer 
%%break spike firing area into a number of place fields
%%work only with lineared track
%%input: rateingrid - spike firing rate in each grid on lineared track grid
%%       zerothreshold - a threshold (in percentage of average firing rate)
%%                    below which the area is not in its place field
%%       maxgap - maximum gap (rate below zero threshold) in a number of pixels within a place field
%%output: numfield - number of place fields
numfield = 0; binsize = xbin(2)-xbin(1); maxgap = round(maxgap/binsize); %%maxgap now in number of bins
gridnum = numel(rateingrid); gridleftindex = 1:gridnum; %%%this is the leftover indices in grid for new field identification
rateingrid = rateingrid - baselinerate; %iii = find(rateingrid<0); if (~isempty(iii)) rateingrid(iii) = zeros(1, numel(iii)); end
[peakratenow, peakindexnow] = max(rateingrid);
while (peakratenow >= minpeakrate) 
    numfield = numfield + 1;
    startindex = []; endindex = []; threstartindex = []; threendindex = []; %%%field locations in rateingrid
    %%%%%determine the left bound of the current field
    zeropoint = 0;
    if (peakindexnow == 1)
        startindex = 1; threstartindex = 1;
    else
        for (k = peakindexnow-1:-1:1)
            if (rateingrid(k) < peakratenow * zerothreshold) | (isnan(rateingrid(k))) 
               zeropoint = zeropoint + 1; 
            else
               zeropoint = 0;
            end
            if (zeropoint > maxgap)  
               startindex = k + floor(zeropoint/2); threstartindex = k + zeropoint; break;
            elseif ((~ismember(k, gridleftindex))|(k==1)) & (zeropoint >= floor(maxgap/2))
               startindex = k + zeropoint-floor(maxgap/2); threstartindex = k + zeropoint; break;
            elseif ((~ismember(k, gridleftindex))|(k==1)) & (zeropoint < floor(maxgap/2))
               startindex = k; threstartindex = k + zeropoint; break;   
            end
        end
    end
    %%%%determine the right bound of the current field
    zeropoint = 0;
    if (peakindexnow == gridnum)
        endindex = gridnum; threendindex = gridnum;
    else
        for (k = peakindexnow+1:gridnum)
            if (rateingrid(k) < peakratenow * zerothreshold) | (isnan(rateingrid(k))) 
               zeropoint = zeropoint + 1; 
            else
               zeropoint = 0;
            end
            if (zeropoint >= maxgap)  
               endindex = k - floor(zeropoint/2); threendindex = k - zeropoint; break;
            elseif ((~ismember(k, gridleftindex))|(k==gridnum)) & (zeropoint >= floor(maxgap/2))
               endindex = k - zeropoint+floor(maxgap/2); threendindex = k - zeropoint; break;
            elseif ((~ismember(k, gridleftindex))|(k==gridnum)) & (zeropoint < floor(maxgap/2))
               endindex = k; threendindex = k - zeropoint; break;   
            end
        end
    end
    %%%%%%characteristics of the current field
    pf(numfield).peakX = xbin(peakindexnow); pf(numfield).startX = xbin(threstartindex); 
    pf(numfield).endX = xbin(threendindex);
    pf(numfield).boundstart = xbin(startindex); pf(numfield).boundend = xbin(endindex);
    pf(numfield).leng = (threendindex-threstartindex+1)*binsize; pf(numfield).size = binsize*sum(rateingrid(threstartindex:threendindex));
    pf(numfield).inmeanrate = baselinerate + mean(rateingrid(threstartindex:threendindex));
    pf(numfield).inpeakrate = peakratenow + baselinerate;
           %%%the center location index and skewness: ATTENTION: rateingrid is a (unnormalized) density function, but not samples
    locmean = NaN; skew = NaN; kurt = NaN;
    if (endindex-startindex > 2)
        sumrate = sum(rateingrid(threstartindex:threendindex));
        if (sumrate >0)
            allloc = xbin(threstartindex:threendindex);  %%field location points
            densfun = rateingrid(threstartindex:threendindex)/sumrate;  %now normalized propability density function
            locmean = allloc * densfun; %sum(allloc .* densfun); 
            locvar = sqrt( ((allloc-locmean).^2) * densfun ); %sqrt( sum( ((allloc-locmean).^2) .* densfun ) );
            if (locvar > 0)
               skew = ((allloc-locmean).^3) * densfun  / (locvar^3);
               X = (allloc-locmean)/locvar; kurt = (X.^4)*densfun - 3;
               %skew = sum( ((allloc-locmean).^3) .* densfun ) / (locvar^3);
               %X = (allloc-locmean)/locvar; kurt = sum((X.^4).*densfun) - 3;
            end
        end
    end
    pf(numfield).comX = locmean; pf(numfield).skew = skew; pf(numfield).kurtosis = kurt;
    %%%%%%leave leftovers for computinng the next field
    gridleftindex = setdiff(gridleftindex, startindex:endindex); %%%leftover grid indices
    ratenow = rateingrid(gridleftindex); [peakratenow, peaknewindex] = max(ratenow); peakindexnow = gridleftindex(peaknewindex);
end

function [mrate, prate, fsize, comX, ploc, skew, kurt] = find1DFFprop(D1map, xxbin, fstart, fend)
mrate = NaN; prate = NaN; fsize=NaN; comX=NaN; ploc=NaN; skew=NaN; kurt = NaN; indstart = 1; indend = 1;
for (i = 2:numel(xxbin))
    if ((xxbin(i)>fstart) && (xxbin(i-1)<=fstart)) 
       indstart = i; break
    end
end
for (i = 2:numel(xxbin))
    if ((xxbin(i)>fend) && (xxbin(i-1)<=fend)) 
       indend = i; break
    end
end
if (indend > indstart)
   ratenow = D1map(indstart:indend); xnow = xxbin(indstart:indend);
   iii = find(ratenow > -10); 
   if (numel(iii) > 2)
       ratenow = ratenow(iii); xnow = xnow (iii); mrate = mean(ratenow); [prate,ii] = max(ratenow); ploc = xnow(ii); 
       fsize = sum(ratenow)*min(diff(xxbin));
       densfun = ratenow/sum(ratenow); comX = xnow * densfun; locvar = sqrt( ((xnow-comX).^2) * densfun );
       if (locvar > 0)
           X = (xnow-comX)/locvar; skew = (X.^3) * densfun; kurt = (X.^4)*densfun - 3;
       end
   end
end

function [crr,pp] = correlateratemaps(rate1, rate2, XY1, XY2, RC, theta)
%%%%here RC and theta are the rotation center and angle: RC[x y](coordinates in pixels), 
%%%% theta (angles in degree): 0, 90, -90, 180, 1 (x mirror), -1 (y mirror)  
%%%% assume XYgrid1 and XYgrid2 (could be 1D or 2D) have the same bin sizes
crr = NaN; pp = NaN;
%%%%first transform the rotated rate maps for 2D maps; for 1D maps: only do the x mirror (only choice 1 is valid... wait... look at directionaly algorithm first)
%%%%for 2D maps: transform XY grid coordinates, then choose common grids and their rates
%%%%%%%%%%determine 1D or 2D rate maps
if (~iscell(XY1)) %%if a 1D rate map: add the cell dimension and RC input replaced by traj2 center
    nd = 1; binsize = XY2(2) - XY2(1);
else
    nd = 2; binsize = XY2{1}(2) - XY2{1}(1);
end
%%%%%first, this is the first map rates and coordinates (x,y), & the second map without rotation/mirror, shifting, etc.
if (nd == 1)
    r1 = rate1; x1 = XY1; r2 = rate2; x2 = XY2; RC = (min(XY2) + max(XY2))/2; 
elseif (nd == 2)
    nx1 = numel(XY1{1}); ny1 = numel(XY1{2});
    r1 = reshape(rate1, nx1*ny1, 1); x1 = zeros(nx1*ny1,1); y1 = zeros(nx1*ny1,1); 
    for (i = 1:nx1)
        for (j = 1:ny1)
            x1((i-1)*ny1+j) = XY1{1}(i); y1((i-1)*ny1+j) = XY1{2}(j);
        end
    end
    nx2 = numel(XY2{1}); ny2 = numel(XY2{2});
    r2 = reshape(rate2, nx2*ny2, 1); x2 = zeros(nx2*ny2,1); y2 = zeros(nx2*ny2,1); 
    for (i = 1:nx2)
        for (j = 1:ny2)
            x2((i-1)*ny2+j) = XY2{1}(i); y2((i-1)*ny2+j) = XY2{2}(j);
        end
    end
end
%%%%%%second: work on the second coordinates (rotated, mirrored, shifted, etc)
%%%%rotation or mirror transformation: shift cartesian origin to RC, then rotate in polar coordinates, then shift cartesian orgin back
%%%%%%%%%%%%shift origin
x2 = x2 - RC(1); if (nd == 2) y2 = y2-RC(2); end
%%%%%%%%%%%%rotate or mirror
if (abs(theta) > 2) %if rotate 2D maps
    [ang, rr] = cart2pol(x2, y2); ang = ang + theta*2*pi/360; [x2, y2] = pol2cart(ang, rr); %coordinates for all rate2 points
elseif (theta == 1) %if mirror X
    x2 = -x2;
elseif (theta == -1) %if mirror y: only available for 2D maps
    if (nd ==2) y2 = -y2; end
end
%%%%%%%%%%%%shift origin back
if (nd ==1) x2 = x2 + RC; end
if (nd == 2) x2 = x2+RC(1); y2 = y2+RC(2); end
%%%%match up (r1, x1, y1) with (r2, x2, y2): coordinate (x y) need to be within binsize 
iii = find(r1>-10); r1 = r1(iii); x1 = x1(iii); if (nd ==2) y1 = y1(iii); end
jjj = find(r2>-10); r2 = r2(jjj); x2 = x2(jjj); if (nd ==2) y2 = y2(jjj); end
nmatch = 0; rr1 = []; rr2 = []; %%%%have to work point-by-point
if (nd ==1)
    for (i = 1:numel(x1))
        [mmm, iii] = min( abs(x2-x1(i)) );
        if (mmm <= binsize)
            nmatch = nmatch + 1; rr1(nmatch) = r1(i); rr2(nmatch) = r2(iii);
        end
    end
elseif (nd == 2)
    for (i = 1:numel(x1))
        [mmm, iii] = min( abs(x2-x1(i)) + abs(y2-y1(i)) );
        if (mmm <= 2*binsize)
            nmatch = nmatch + 1; rr1(nmatch) = r1(i); rr2(nmatch) = r2(iii);
        end
    end
end
if (nmatch>=10) %%%%if overlap more than 10 bins
    m1 = mean(rr1); m2 = mean(rr2); s1 = std(rr1); s2 = std(rr2);
    if (s1~=0)&&(s2~=0)
        [R,P] = corrcoef(rr1',rr2'); crr = R(1,2); pp = P(1,2);
        %rr1 = (rr1-m1)/s1; rr2 = (rr2-m2)/s2; crr = rr1*rr2'/nmatch; 
    end
end

%%%%%%%%%%%%%%The following needs to be finalized. Disable now!!!!%%%%%%%%%%%
function [D1rate, D1refmap, D1avgmap, xxbin] = alignlaps(D1ratenow, D1refmap, evD1rate, evOccutime, evNum, xbin, xrefbin)
nev = numel(xbin); nlap = numel(evNum);
D1rate = cell(1, nlap);
%%%%%D1ratenow{nlap}; evNum(nlap); xbin{nev}; xrefbin(nbin)
% %%%%%align at the location centers
% ind = cell(1,nev); nind = zeros(1, nev); nanid = cell(1, nev);
% for (i = 1:nev)
%      xbinnow = xbin{i}; xbinow = xbinnow(~isnan(xbinnow)); nanid{i} =isnan(xbinnow);
%      xbin{i} = xbin{i} - mean(xbinnow); mmin(i) = min(xbinnow); mmax(i) = max(xbinnow);
% end
% mm1 = max(mmin); mm2 = min(mmax);
% for (i = 1:nev)
%     ind{i} = find( (xbin{i}>=mm1) & (xbin{i}<=mm2) ); ind{i} = sort(union(ind{i}, nanid{i}));
%     nind(i) = numel(ind{i});
% end
% xrefbin = xrefbin - mean(xrefbin);
% refind = find( (xrefbin>=mm1) & (xrefbin<=mm2) ); nrefind = numel(refind);
% D1ref = zeros(mind,1);  
% if (nrefind~=mind) disp('--------------> warning: lap position dimensions do not match; firing rates stuffed with zeros'); end
% D1ref(1:nrefind) = D1refmap(refind); D1refmap = D1ref;
%%%%%align at the bin centers
midind = ones(1,nev); nind = zeros(1, nev); ind = cell(1, nev); shiftind = cell(1, nev);
for (i = 1:nev)
    %disp(['----number trajbins: ', num2str(numel(xbin{i}))]);
    %disp(['----number of traj rate points: ', num2str(numel(D1ratenow{i}))]);
     %%%%%%Warning: %%%somehow traj bin and rate map bins do not match after stopping removal
%      binnumnow = zeros(1, numel(D1ratenow{i})+1);
%      binnumnow(1) = numel(xbin{i});
%      for (j = 1:numel(evD1rate{i}))
%          binnumnow{j+1} = numel(evD1rate{i}{j}); 
%      end
%      ind{i} = 1:min(binnumnow); 
disp(evD1rate{i}')
disp(['***********'])
disp(D1ratenow{i}')
disp(['&&&&&&&'])
     ind{i} = 1:min([numel(xbin{i}) numel(evD1rate{i})]);
     nind(i) = numel(ind{i}); midind(i) = round(mean(ind{i})); 
end
[mind, ttt] = max(nind); xxbin = xbin{ttt}(ind{ttt}); fallind = 1:mind; fmidind = round(mean(fallind)); 
for (i = 1:nev)
    shiftind{i} = fallind( find( (fallind-fmidind>=1-midind(i)) & (fallind-fmidind<=nind(i)-midind(i)) ) );
end
for (k = 1:nlap)
    i = evNum(k); %%%%this lap's ev number
    if (nind(i)~=mind) disp('--------------> warning: lap position dimensions do not match; firing rates stuffed with zeros'); end
    D1rate{k} = zeros(mind,1);
    D1rate{k}(shiftind{i}) = D1ratenow{i}{k}(ind{i});
end
nrefind = numel(xrefbin); refind = 1:nrefind; midrefind = round(mean(refind)); 
refshiftind = find( (fallind-fmidind>=1-midrefind) & (fallind-fmidind<=nrefind-midrefind) );
D1ref = zeros(mind,1);  
if (nrefind~=mind) disp('--------------> warning: reference lap position dimensions do not match; firing rates stuffed with zeros'); end
D1ref(refshiftind) = D1refmap; D1refmap = D1ref;
%%%%%%The following compute the avgratemap = real summed-up map: totalspikes/totaltime
for (i = 1:nev)
    temrate = zeros(mind, 1); temtime = zeros(mind,1);
    temrate(shiftind{i}) = evD1rate{i}; evD1rate{i} = temrate;
    temtime(shiftind{i}) = evOccutime{i}; evOccutime{i} = temtime;
end
if (nev > 0)
   if (nev == 1)
      D1avgmap = evD1rate{1}; 
   else
      D1avgmap = NaN*ones(mind, 1);
      nspike = zeros(nev,mind); ntime = zeros(nev,mind); totalspike = NaN*ones(1,mind); totaltime = zeros(1,mind);
      for (k = 1:nev)
           nspike(k,:) = evD1rate{k}; ntime(k,:) = evOccutime{k};
      end
      nspike = nspike .* ntime; totaltime = sum(ntime);
      for (m = 1:mind)
          iii = find( nspike(:,m)>-10 ); 
          if (~isempty(iii)) totalspike(m) = sum(nspike(iii,m)); end
          D1avgmap(m) = totalspike(m)/totaltime(m);
      end
   end
end

function evSess = identifysession(evTime, sessionname, startT, endT)
evSess = []; minevstart = min(evTime.start); maxevend = max(evTime.ent);
iii = find( (startT <= minevstart) & (endT >= maxevend) );
if (numel(iii) == 1)
    evSess = sessionname{iii};
end

function [D2rate, nact] = get2Dratemaps(spiketime, XYbin, gridoccup, postime, XX, YY, smParm)
nx = numel(XYbin{1}); ny = numel(XYbin{2}); D2rate = NaN*ones(ny,nx); nact = NaN;
if (nx > 1) && (ny > 1)
    binsize = XYbin{1}(2) - XYbin{1}(1); sigma = round(smParm.d2sigma/binsize); %%%sigma now in number of bins
    spiketime = sort(spiketime); [postime, iii] = sort(postime); XX = XX(iii); YY = YY(iii); 
    lastpoint = 1;  %this only for saving time
    x = NaN*ones(size(spiketime)); y = NaN*ones(size(spiketime));
    for (j = 1:numel(spiketime))
        for (k = lastpoint:numel(postime)-1) %find corresponding time in position data
            if (postime(k) <= spiketime(j)) && (postime(k+1) > spiketime(j)) 
                x(j) = XX(k); y(j) = YY(k); lastpoint = k; break; 
            end
        end
    end
    %%%%%count spikes
   countx = cell(1, nx); county = cell(1, ny); count = zeros(ny, nx); %%%this asignment is for convenient plot using pcolor(x-columnx, y-rows)
   for (i = 1:nx-1) countx{i} = find((x>=XYbin{1}(i)) & (x<XYbin{1}(i+1))); end
   for (i = 1:ny-1) county{i} = find((y>=XYbin{2}(i)) & (y<XYbin{2}(i+1))); end
   countx{nx} = find(x>=XYbin{1}(nx)); county{ny} = find(y>=XYbin{2}(ny));
   %%%%% now spike count at grid (i,j) = number of common indices in countx{i} and county{j}
   for (i = 1:nx)
   for (j = 1:ny)
       common = []; count(j,i) = 0;
       if (~isempty(countx{i}) & ~isempty(county{j}))
          common = intersect(countx{i}, county{j});
          if (~isempty(common)) count(j,i)=numel(common); end  %%assign (j,i) because pcolor requires 
       end
   end
   end
   %%%%%compute firing rate
   for (j = 1:ny)
   for (i = 1:nx)
       if (gridoccup(j,i) > 0) D2rate(j,i) = count(j,i)/gridoccup(j,i); end
   end
   end
   nact = numel(find(count>0))/numel(find(gridoccup>0));
   %%now do a 2d smoothing
   if (strcmpi(smParm.d2sm, 'yes'))
      D2rate = TwoDSmooth_new(D2rate, gridoccup, sigma, sigma, smParm.d2Nsig, smParm.d2Nsig);
   end
end

function [D1rate, nact] = getlinearmap(spiketime, evTime, lappostime, lapx, xbin, occutime, smParm)
nx = numel(xbin); D1rate = NaN*ones(nx,1); nlap = numel(lappostime); count = zeros(nx,1); spikex = []; nact = NaN;
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
    D1rate = OneDSmooth_new(D1rate, occutime, sigma, smParm.d1Nsig);
    %XX = [-smParm.d1Nsig*sigma:1:smParm.d1Nsig*sigma]; 
    %wind = normpdf(XX, 0, sigma); wind = wind/sum(wind); nw = numel(wind); %will have a delay = (n-1)/2 bins
    %D1rate = conv(D1rate, wind); D1rate = D1rate((nw-1)/2+1:numel(D1rate)-(nw-1)/2);
end
end

function [spinfo, spinfotime, sparsty] = findratemapprop(rate, occu)
%%%%%%another variable: mutual information is not computed here
spinfo = NaN; spinfotime = NaN; sparsty = NaN;  
[ny, nx] = size(rate); rate = reshape(rate, nx*ny, 1); occu = reshape(occu, nx*ny, 1); mr = max(rate);
if (sum(occu) > 0)
   occu = occu / sum(occu); jjj = find(rate>-1); meanrate = sum(occu(jjj).*rate(jjj));
   spinfo = 0; spinfotime = 0; sparsty = 0;  
   if (meanrate > 0)
       rateok = rate / meanrate; iii = find(rateok > 0); 
       spinfo = sum(occu(iii).*rateok(iii).*log2(rateok(iii))); 
       spinfotime = spinfo*meanrate; sparsty = 1/sum(occu(iii).*rateok(iii).*rateok(iii));
   end
end

function baserate = computebaserate(ratenow, baserategridN)
baserate = 0;
ratenow = sort(ratenow(ratenow>-10)); nn = round(baserategridN*numel(ratenow));
if (nn > 0) baserate = mean(ratenow(1:nn)); end

function [D2refmap, XYrefbin, refsegNumout] = resolve2Drefmap(segSess, segNum, withinsegNum, D2ratenow, XYsess, XYbin, BaseSess, BaseSeg)
D2refmap = []; XYrefbin = []; refsegNumout = NaN;
ises = find( strcmpi(segSess,BaseSess) ); ibin = find( strcmpi(XYsess, BaseSess) );
if isempty(ises) || isempty(ibin)
    disp(['----------------> 2D dynamics not computed: no reference sessions found: ', BaseSess]);
else
    mapnow = D2ratenow(ises); withinsegnow = withinsegNum(ises); segNum = segNum(ises);
    refsegnum = str2num(BaseSeg);  
    if (isempty(refsegnum))
        if (strncmpi(BaseSeg, 'last', 4))
            [str, tok] = strtok(lower(BaseSeg), 'last');
            if (~isempty(str2num(str))) 
                refsegnum = (max(withinsegnow)-(str2num(str)-1)):max(withinsegnow); 
            else
                [str, tok] = strtok(str, '_');
                if (~isempty(str2num(str))) refsegnum = (max(withinsegnow)-(str2num(str)-1)):max(withinsegnow);  end
            end
        elseif (strncmpi(BaseSeg, 'first', 5))
            [str, tok] = strtok(lower(BaseSeg), 'first');
            if (~isempty(str2num(str))) 
                refsegnum = min(withinsegnow):(min(withinsegnow)+(str2num(str)-1)); 
            else
                [str, tok] = strtok(str, '_');
                if (~isempty(str2num(str))) refsegnum = min(withinsegnow):(min(withinsegnow)+(str2num(str)-1));  end
            end 
        end
    end
    if (~isempty(refsegnum)) || (strncmpi(BaseSeg, 'average', 4)) %%%if average or segments specified, do an average map
        if (~isempty(refsegnum)) 
            [~, iw, ir] = intersect( withinsegnow, refsegnum ); mapnow = mapnow(iw); segNum = segNum(iw);
        end
        refsegNumout = segNum;
        if (~isempty(mapnow))
            XYrefbin = XYbin{ibin(1)};
            D2refmap = mapnow{1}; [mm,nn] = size(D2refmap);
            for (i = 1:mm)
                for (j = 1:nn)
                    for (k = 2:numel(mapnow));
                        D2refmap(i,j) = D2refmap(i,j) + mapnow{k}(i,j);
                    end
                    D2refmap(i,j) = D2refmap(i,j)/numel(mapnow);
                end
            end
        else
           disp(['----------------> 2D dynamics not computed: no reference session segments found: ', BaseSess, '__', BaseSeg]);
        end
    elseif (strcmpi(BaseSeg, 'last')) %%if last segment of a session
        XYrefbin = XYbin{ibin(1)}; D2refmap = mapnow{numel(mapnow)}; refsegNumout = segNum(numel(mapnow));
    elseif (strcmpi(BaseSeg, 'first')) %%if first segment of a session
        XYrefbin = XYbin{ibin(1)}; D2refmap = mapnow{1}; refsegNumout = segNum(1);
    else
        disp(['----------------> 2D dynamics not computed: no reference session segments found: ', BaseSess, '__', BaseSeg]);
    end
end

function [D1refmap, xrefbin, reflapNum] = resolve1Drefmap(lapevName, lapNum, withinlapNum, D1ratenow, xev, xbin, BaseSess, BaseSeg)
D1refmap = []; xrefbin = []; reflapNum = NaN; nev = numel(xev); nlap = numel(lapNum); %%%BaseSeg = BaseLap; use the name seg for lap
evid = zeros(1, nev); lapid = zeros(1, nlap);
for (i = 1:nev)
    if (~isempty(strfind(lower(xev{i}), lower(BaseSess)))) evid(i) = 1; end
end
for (i = 1:nlap)
    if (~isempty(strfind(lower(lapevName{i}), lower(BaseSess)))) lapid(i) = 1; end
end
ises = find( lapid ); ibin = find( evid );
if (isempty(ises) || isempty(ibin))
    disp(['----------------> 1D dynamics not computed: no reference sessions found: ', BaseSess]);
else
    mapnow = D1ratenow(ises); withinlapnow = withinlapNum(ises); lapNum = lapNum(ises);
    refsegnum = str2num(BaseSeg);  
    if (isempty(refsegnum))
        if (strncmpi(BaseSeg, 'last', 4))
            [str, tok] = strtok(lower(BaseSeg), 'last');
            if (~isempty(str2num(str))) 
                refsegnum = (max(withinlapnow)-(str2num(str)-1)):max(withinlapnow); 
            else
                [str, tok] = strtok(str, '_');
                if (~isempty(str2num(str))) refsegnum = (max(withinlapnow)-(str2num(str)-1)):max(withinlapnow);  end
            end
        elseif (strncmpi(BaseSeg, 'first', 5))
            [str, tok] = strtok(lower(BaseSeg), 'first');
            if (~isempty(str2num(str))) 
                refsegnum = min(withinlapnow):(min(withinlapnow)+(str2num(str)-1)); 
            else
                [str, tok] = strtok(str, '_');
                if (~isempty(str2num(str))) refsegnum = min(withinlapnow):(min(withinlapnow)+(str2num(str)-1));  end
            end 
        end
    end
    if (~isempty(refsegnum)) || (strncmpi(BaseSeg, 'average', 4)) %%%if average or segments specified, do an average map
        if (~isempty(refsegnum)) 
           [~, iw, ir] = intersect( withinlapnow, refsegnum ); mapnow = mapnow(iw); lapNum = lapNum(iw);
        end
        reflapNum = lapNum;
        if (~isempty(mapnow))
            xrefbin = xbin{ibin(1)}; D1refmap = mapnow{1}; nmap = numel(mapnow);
            for (k = 2:nmap);
                 D1refmap = D1refmap + mapnow{k};
            end
            D1refmap = D1refmap/nmap;
        else
            disp(['----------------> 1D dynamics not computed: no reference session segments found: ', BaseSess, '__', BaseSeg]);
        end
    elseif (strcmpi(BaseSeg, 'last')) %%if last segment of a session
        xrefbin = xbin{ibin(1)}; D1refmap = mapnow{numel(mapnow)}; reflapNum = lapNum(numel(mapnow));
    elseif (strcmpi(BaseSeg, 'first')) %%if first segment of a session
        xrefbin = xbin{ibin(1)}; D1refmap = mapnow{1}; reflapNum = lapNum(1);
    else
        disp(['----------------> 1D dynamics not computed: no reference session segments found: ', BaseSess, '__', BaseSeg]);
    end
end


