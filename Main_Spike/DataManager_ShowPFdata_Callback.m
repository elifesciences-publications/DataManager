function DataManager_ShowPFdata_Callback
%%Plot the computed bhdata
hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
hfield = getappdata(hf, 'hfield'); plotparm = getappdata(hf, 'plotparm');
hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); 
spikeind = find(spikeselection == 1); %%%%%spikeind contains what will be plotted
behav = []; bhdata = [];
if (plotparm.linkbehav == 0);
    disp(['--------> no behavioral data linked']);
else
    behav = getappdata(hf, 'behav'); bhdata = getappdata(hf, 'bhdata');
end

tagmark = get(gcbo, 'Tag');
if (strcmp(tagmark, 'show1dfields'))
    plot1dfields(pinfo, data, behav, bhdata, spikeind, 'Show1Dfields');
elseif (strcmp(tagmark, 'show2dfields'))
    plot2dfields(pinfo, data, behav, bhdata, spikeind, 'Show2Dfields');
elseif (strcmp(tagmark, 'show1dspikes'))
    plot1dspikes(pinfo, data, behav, bhdata, spikeind, 'Show1DSpikes'); %%%wait until the 1D field dynam has been done
elseif (strcmp(tagmark, 'show1dphases'))
    plot1dfields(pinfo, data, behav, bhdata, spikeind, 'Show1DPhases');     
% elseif (strcmp(tagmark, 'evheaddir'))
%     plotevVelDir(behav, bhdata, dateind, 'evHeaddir');
% elseif (strcmp(tagmark, 'sesstraj'))
%     plotsesstraj(behav, bhdata, dateind, 'sessTraj');
% elseif (strcmp(tagmark, 'evtraj'))
%     plotevVelDir(behav, bhdata, dateind, 'evTraj');
% elseif (strcmp(tagmark, 'sessoccup'))
%     plotsesoccup(behav, bhdata, dateind, 'sessOccupancy');
% elseif (strcmp(tagmark, 'evoccup'))
%     plotevoccup(behav, bhdata, dateind, 'evOccupancy');
% elseif (strcmp(tagmark, 'showposdata')) %%%for all other options to plot data
%     showalldata(behav, bhdata, dateind, 'allposdata');
end
disp('************************');

function plot2dfields(pinfo, data, behav, bhdata, spikeind, tag)
if ~isfield(pinfo, 'field')
    msgbox('Fields not defined in the data base');
else
    if (numel(spikeind)>10)
        disp('--------> Too many cells selected; Only plot fields for the first 10 cells');
        spikeind = spikeind(1:10);
    end
    for (tt = 1:numel(spikeind))
        i = spikeind(tt); %%for the current cell, get event rate maps and defined fields
        clname = pinfo.general.clname{i}; sessName = pinfo.general.sessionname{i};
        D2rate = data.field.sess2DRateMaps{i}; 
        for (j = 1:numel(sessName))
             [ny,nx] = size(D2rate{j}); 
             if (nx >1) & (ny > 1)
                finaldirnow = pinfo.general.finaldir{i}; 
                xybin{1} = 5*(1:nx); xybin{2} = 5*(1:ny);  %default xybin
                if (~isempty(behav)) & (~isempty(bhdata)) %%%get xbin from behavioral data
                   posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, sessName{j}) );
                   if (numel(posid) ~= 1)
                      disp(['-------------> field not computed: no or more than 1 positon files match the session: ', finaldirnow, '___', sessnow]);
                   else
                      xybin = bhdata.sess.gridXYbin{posid}; 
                   end
                end
                %%%%%%%%%pre-treat the rate map (NaN values)
                mr = max(max(D2rate{j})); if isempty(mr)|isnan(mr)|(mr ==0 ) mr = 1; end
                rate = zeros(ny,nx);
                for (tj = 1:ny)
                for (ti = 1:nx)
                    if (isnan(D2rate{j}(tj,ti)))
                        rate(tj,ti) = -mr/63; %%%this is because the coloarmap has 64 colors, this makes sure the lowest value pixels are white
                    else
                        %if (D2rate{j}(tj,ti) == 0) % <= mr/64) 
                        %    rate(tj,ti) = 1; %mr/64;
                        %else
                            rate(tj,ti) = D2rate{j}(tj,ti);
                        %end
                    end
                end
                end
                %%%plot rate maps and xjoint
                hg = figure('Name', strcat(tag, '---', clname)); 
                hax = axes('Parent', hg, 'NextPlot', 'add', 'CLimMode', 'manual', 'CLim', [-mr/63 mr]); 
                cmap = colormap(jet); cmap(1,:) = [1 1 1]; %change background to white
                hp = pcolor(xybin{1}, xybin{2}, rate); colormap(cmap); colorbar('vert', 'peer', hax);
                xlabel ('X (pixel)'); ylabel('Y (pixel)'); 
                xlim([min(xybin{1})-10 max(xybin{1})+10]); ylim([min(xybin{2})-10 max(xybin{2})+10]);
                set(hp, 'EdgeColor', 'none', 'LineStyle', 'none', 'Parent', hax);
                %%%find and plot fields belong to this event (trajectory)
                pfsess = pinfo.field.PF2Dsess{i}; pfid = find(strcmp(pfsess, sessName{j}));
                for (k = 1:numel(pfid))
                    plocx = pinfo.field.PF2DpeakX{i}(pfid(k)); plocy = pinfo.field.PF2DpeakY{i}(pfid(k)); 
                    prate = pinfo.field.PF2DInPeakrate{i}(pfid(k));
                    strnow = ['*(',num2str(prate),'Hz)'];
                    text('Interpreter', 'none', 'Parent', hax, 'String', strnow, 'Color', [0 1 0], 'Units', 'data', 'Position', [plocx plocy]); 
                end
                str = ['Session: ', sessName{j}, '; number of fields identified: ', num2str(numel(pfid))];
                text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.97]); 
             end
        end
    end
end

function plot1dfields(pinfo, data, behav, bhdata, spikeind, tag)
if ~isfield(pinfo, 'field')
    msgbox('Fields not defined in the data base');
else
    if (numel(spikeind)>10)
        disp('----------> Too many cells selected; Only plot fields for the first 10 cells');
        spikeind = spikeind(1:10);
    end
    for (tt = 1:numel(spikeind))
        i = spikeind(tt); %%for the current cell, get event rate maps and defined fields
        clname = pinfo.general.clname{i};
        evName = pinfo.general.eventname{i}; evTime = data.events.eventtimes{i}; evType = pinfo.parm.eventtype{i}; 
        D1rate = data.field.evt1DRateMaps{i}; finaldirnow = pinfo.general.finaldir{i};
        spiketime = pinfo.parm.timeunit(i)*data.spike.spiketime{i};
        for (j = 1:numel(evTime))
             xbin = 5*(1:numel(D1rate{j})); xjoint = []; nlap = 0; %default xbin
             if (strcmp(evType{j}, 'run')) && (~isempty(D1rate{j}))
                if (~isempty(behav)) & (~isempty(bhdata)) %%%get xbin from behavioral data
                   %%%%locate event position data
                   evSess = identifysession(evTime{j}, pinfo.general.sessionname{i}, pinfo.general.sessionstartT{i}, pinfo.general.sessionendT{i});
                   posid = []; evid = [];
                   if (~isempty(evSess))
                      posid = find( strcmp(behav.general.finaldir, finaldirnow) & strcmp(behav.general.sessname, evSess) );
                   end
                   if numel(posid == 1)
                      if (isfield(behav.general, 'eventname'))
                          evid = find(strcmp(behav.general.eventname{posid}, evName{j}));
                      else
                          evid = find(strcmp(behav.behavior.eventname{posid}, evName{j}));
                      end
                   end
                   if (numel(posid)~=1)|(numel(evid)~=1)
                      disp(['-------------> field not computed: no or more than 1 positon/event files match the session: ', finaldirnow, '___', evName{j}]);
                   else
                      xbin = bhdata.event.Xbin{posid}{evid}; xjoint = bhdata.event.Xjoint{posid}{evid};
                   end
                   %%%%get lap-by-lap spike rasters
                   nlap = numel(evTime{j}.start); spikeX = cell(1, nlap); allspikeid = []; allX = [];
                   for (tj = 1:nlap)
                        lappostime = bhdata.event.LapAllPostimestamp{posid}{evid}{tj}; 
                        lapx = bhdata.event.LapAllX{posid}{evid}{tj};
                        evok.start = evTime{j}.start(tj); evok.ent = evTime{j}.ent(tj);
                        [spiketimenow, epid, spikeidnow] = SpikeEventFilter(spiketime, evok);
                        spikeX{tj} = findspikex(lappostime, lapx, spiketimenow); 
                        allspikeid = [allspikeid; spikeidnow]; allX = [allX; spikeX{tj}];
                   end
                   if strcmp(tag, 'Show1DPhases')
                        allP = data.phase.SpikePhases{i}(allspikeid);
                   end
                end
                %%%%%%%plot spike rasters
                hg = figure('Name', strcat(tag, '---', clname)); posvec = [0.08 0.05 0.90 0.90]; mr = max(D1rate{j}); if (mr ==0 ) mr = 1; end
                if strcmp(tag, 'Show1DPhases')
                   posvecnow = [posvec(1) posvec(2)+posvec(4)/3 posvec(3) posvec(4)*2/3]; 
                   haxok = axes('Parent', hg, 'Units', 'normalized', 'Position', posvecnow, 'NextPlot', 'add', 'XLim', [min(xbin)-10 max(xbin)+10], 'YLim', [-360 720]);
                   line(allX, allP, 'Parent', haxok, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
                   line(allX, allP-360, 'Parent', haxok, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
                   line(allX, allP+360, 'Parent', haxok, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
                else  
                   for (tti=1:nlap) %first tick trace is on the top.
                     mtickposvec{tti} = [posvec(1) posvec(2)+posvec(4)/3+(nlap-tti)*posvec(4)*2/3/nlap posvec(3) posvec(4)*2/3/nlap];    % multiple tick traces
                     hmtickaxes(tti) = axes('Parent', hg, 'Units', 'normalized', 'Position', mtickposvec{tti}, 'Visible', 'off', 'NextPlot', 'add');
                     xlim ([min(xbin)-10 max(xbin)+10]); %plot x-axis as linearized track
                     ylim ([0 1.2]);
                   end
                   for (tti=1:nlap) %for each episode
                    if (~isempty(spikeX{tti}))
                    for (ttj=1:numel(spikeX{tti}))
                        line([spikeX{tti}(ttj) spikeX{tti}(ttj)], [0 1], 'Parent', hmtickaxes(tti), 'LineWidth', 2, 'Color', [0 0 1]) % plot spikes as ticks
                    end
                    end
                   end
                end
                %%%plot rate maps and xjoint
                rateposvec = [posvec(1) posvec(2) posvec(3) posvec(4)/3]; 
                hax = axes('Parent', hg, 'Units', 'normalized', 'Position', rateposvec, 'NextPlot', 'add', 'XLim', [min(xbin)-10 max(xbin)+10], 'YLim', [-0.2 1.2]*mr);
                xlabel ('X (pixel)'); ylabel('Firing rate (Hz)'); 
                line(xbin, D1rate{j}, 'Parent', hax, 'LineWidth', 4, 'Color', [0 0 0]);
                for (k = 1:numel(xjoint))
                    line([xjoint(k) xjoint(k)], [0 max(D1rate{j})], 'Parent', hax, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
                end
                %%%find and plot fields belong to this event (trajectory)
                pfev = pinfo.field.PF1Devt{i}; pfid = find(strcmp(pfev, evName{j})); baserate = pinfo.field.PF1DBaseRate{i}(j);
                for (k = 1:numel(pfid))
                    sbound = pinfo.field.PF1DBoundStart{i}(pfid(k)); ebound = pinfo.field.PF1DBoundEnd{i}(pfid(k)); 
                    sloc = pinfo.field.PF1DLocStartX{i}(pfid(k)); eloc = pinfo.field.PF1DLocEndX{i}(pfid(k));
                    ploc = pinfo.field.PF1DLocPeakX{i}(pfid(k)); cloc = pinfo.field.PF1DLocComX{i}(pfid(k));
                    prate = pinfo.field.PF1DInPeakrate{i}(pfid(k));
                    line([sbound sbound], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [0 0 1]);
                    line([ebound ebound], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [0 0 1]);
                    line([sbound ebound], [prate prate], 'Parent', hax, 'LineWidth', 1, 'Color', [1 0 0]);
                    line([sloc sloc], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [1 0 0]);
                    line([eloc eloc], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [1 0 0]);
                    line([ploc ploc], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [0 1 0]);
                    line([cloc cloc], [0 prate], 'Parent', hax, 'LineWidth', 1, 'Color', [0 1 1]);
                    if strcmp(tag, 'Show1DPhases')
                       S = pinfo.phase.fieldPPSlope{i}(pfid(k)); B = pinfo.phase.fieldPPInt{i}(pfid(k));
                       xx = [sbound ebound]; yy = xx*S + B;    
                       line(xx, yy, 'Parent', haxok, 'LineWidth', 2);
                    end
                end
                line([xbin(1) xbin(numel(xbin))], [baserate baserate], 'Parent', hax, 'LineWidth', 1, 'Color', [0 0 0]);
                str = [evName{j} '; Number of fields identified: ' num2str(numel(pfid))];
                text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.94]);
                %text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.05 0.92]);   
                
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

function spikeX = findspikex(postime, posx, spiketime)
spikeX = zeros(size(spiketime));
[postime, iii] = sort(postime); posx = posx(iii);
lastpoint = 1;  %this only for saving time
for (j = 1:numel(spiketime))
    for (k = lastpoint:numel(postime)-1) %find corresponding time in position data
         if (postime(k) <= spiketime(j)) && (postime(k+1) > spiketime(j)) 
             spikeX(j) = posx(k); lastpoint = k; break; 
         end
    end
end


