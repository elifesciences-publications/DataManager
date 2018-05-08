function DataManager_ShowCrrdata_Callback
%%Plot the computed bhdata
avgplottype = 'median';   %%%%"median" for median [25% 75%] range or "mean" for mean +- se
ColorPlotType = 'centralpeakorder';
hf = gcbf; pinfo = getappdata(hf, 'pinfo'); data = getappdata(hf, 'data');
hgroup = getappdata(hf, 'hgroup'); groupselection = getappdata(hgroup, 'selection');
hfield = getappdata(hf, 'hfield'); plotparm = getappdata(hf, 'plotparm');
hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection'); 

tagmark = get(gcbo, 'Tag');
thre = 0;
if (strcmp(tagmark, 'showsessioncrr'))
    cellind = find(spikeselection == 1); %%%%%spikeind contains what will be plotted
    plotcrosscrr(pinfo, data, cellind, 'session', 'sessCrossCrr');
elseif (strcmp(tagmark, 'showevcrr'))
    cellind = find(spikeselection == 1); %%%%%spikeind contains what will be plotted
    plotcrosscrr(pinfo, data, cellind, 'ev', 'evCrossCrr');
elseif (strcmp(tagmark, 'showsessaveragecrr'))
    groupsel = find(groupselection == 1); grpind = cell(1, numel(groupsel)); grpname = cell(1, numel(groupsel));
    for (i = 1:numel(groupsel))
        grpind{i} = data.grouplist.groupindex{groupsel(i)};
        grpname{i} = data.grouplist.groupname{groupsel(i)};
    end    
    plotaveragecrosscrr(pinfo, data, grpname, grpind, 'session', 'sessAvgCrossCrr', avgplottype);
elseif (strcmp(tagmark, 'showevaveragecrr'))
    groupsel = find(groupselection == 1); grpind = cell(1, numel(groupsel)); grpname = cell(1, numel(groupsel));
    for (i = 1:numel(groupsel))
        grpind{i} = data.grouplist.groupindex{groupsel(i)};
        grpname{i} = data.grouplist.groupname{groupsel(i)};
    end  
    plotaveragecrosscrr(pinfo, data, grpname, grpind, 'ev', 'evAvgCrossCrr', avgplottype);
elseif (strcmp(tagmark, 'showevcrrcolor'))
    groupsel = find(groupselection == 1); grpind = cell(1, numel(groupsel)); grpname = cell(1, numel(groupsel));
    for (i = 1:numel(groupsel))
        grpind{i} = data.grouplist.groupindex{groupsel(i)};
        grpname{i} = data.grouplist.groupname{groupsel(i)};
    end  
    plotcrosscrrcolor(pinfo, data, grpname, grpind, 'ev', 'evCrossCrrColor', ColorPlotType);
%elseif (strcmp(tagmark, 'show1dspikes'))
%    plot1dspikes(pinfo, data, behav, bhdata, spikeind, 'Show1DSpikes'); %%%wait until the 1D field dynam has been done
% elseif (strcmp(tagmark, 'evheaddir'))
%     plotevVelDir(behav, bhdata, dateind, 'evHeaddir');

end
disp('************************');

function plotcrosscrrcolor(pinfo, data, grpname, grpind, sesstype, tag, plottype)
%%%determine normalization factor threhold if crr is in normcount mode
thre = 0; ok = 1;
if strcmp(pinfo.parm.crrmode{grpind{1}(1)}, 'normcount') || strcmp(pinfo.parm.crrmode{grpind{1}(1)}, 'basenorm')
    III=inputdlg('Threshold for normalization factor (0 for all)', 'Pair selection', 1, {'0'}, 'on');
    if (~isempty(III))
        [thre, ok] = str2num(III{1});
    else
        ok = 0;
    end
end
if ~ok
    disp('-------------> Invalid threshold; aborted.')
else

if (strcmp(sesstype, 'session'))
    snfield = 'sessionname';  typefield = 'sessType'; crrfield = 'sesscrr'; crrnormfield = 'sessCrrNorm'; 
    istr = 'Session '; 
elseif (strcmp(sesstype, 'ev'))
    snfield = 'eventname'; typefield = 'eventtype'; crrfield = 'evcrr'; crrnormfield = 'evCrrNorm';
    istr = 'Event ';
end
keyword = []; keytype = []; tit = [istr ' selection (empty for averaging all ', istr, 's']; TT = [istr ' type:'];
III=inputdlg({'Keyword'; TT}, tit, 2, {'Sleep1'; 'ripple'}, 'on');
if (~isempty(III))
   keyword = III{1}; keytype = III{2}; ook = 1;
else
   ook = 0;
end
if ook
for (tt = 1:numel(grpind))
    spikeind = grpind{tt}; nspike = numel(spikeind); clname = grpname{tt};
    %%%%%determine common sessions and common timebin
%     if (strcmp(sesstype, 'session'))
%        if isfield(pinfo.parm, 'sesskeyword') keyword = pinfo.parm.sesskeyword{spikeind(1)}; end
%        if isfield(pinfo.parm, 'sesskeytype') keytype = pinfo.parm.sesskeytype{spikeind(1)}; end
%     elseif (strcmp(sesstype, 'ev'))
%        if isfield(pinfo.parm, 'evtkeyword') keyword = pinfo.parm.evtkeyword{spikeind(1)}; end
%        if isfield(pinfo.parm, 'evtkeytype') keytype = pinfo.parm.evtkeytype{spikeind(1)}; end
%     end

    if isempty(keyword) && isempty(keytype) %%%if no key word/type specified
       sessnow = pinfo.general.(snfield){spikeind(1)}; 
       xbin = data.crr.timebin{spikeind(1)};
       for (i = 2:numel(spikeind))
           sessnow = intersect(sessnow, pinfo.general.(snfield){spikeind(i)});
           xbin = intersect(xbin, data.crr.timebin{spikeind(i)});
       end
       if (~isempty(sessnow)) && (~isempty(xbin))
        for (i = 1:numel(sessnow))
            sessName = sessnow{i}; nbin = numel(xbin); allcrr = NaN*ones(nspike, nbin); allnorm = NaN*ones(nspike,1);
            for (j = 1:nspike)
                [~,binind] = intersect(data.crr.timebin{spikeind(j)}, xbin);
                nses = find(strcmp(pinfo.general.(snfield){spikeind(j)}, sessName));
                allcrr(j,:) = data.crr.(crrfield){spikeind(j)}{nses}(binind);
                allnorm(j) = pinfo.crr.(crrnormfield){spikeind(j)}{nses};
            end
            if (thre > 0) iij = find(allnorm >= thre); allcrr = allcrr(iij,:); end 
            plotColorcrr(xbin, allcrr, nbin, tag, clname, sessName, sesstype, keyword, keytype, plottype);
        end
       else
        disp('-------------> warning: no common sessions or common time lags found');
       end
    else %%%%%if keyword or keytype are specified
       xbin = data.crr.timebin{spikeind(1)}; sessName = [];
       for (i = 2:numel(spikeind))
            xbin = intersect(xbin, data.crr.timebin{spikeind(i)});
       end
       if (~isempty(xbin))
            nbin = numel(xbin); allcrr = NaN*ones(nspike, nbin); 
            for (j = 1:nspike)
                [~,binind] = intersect(data.crr.timebin{spikeind(j)}, xbin);
                nses = numel(pinfo.general.(snfield){spikeind(j)}); sessel = ones(1, nses); normsel = ones(1, nses);
                for (ik = 1:nses)
                    if ~isempty(keyword)
                       if isempty(strfind(lower(pinfo.general.(snfield){spikeind(j)}{ik}), lower(keyword))) 
                           sessel(ik) = 0;
                       end
                    end
                    if ~isempty(keytype)
                       if ~strcmpi(pinfo.parm.(typefield){spikeind(j)}{ik}, keytype) sessel(ik) = 0; end
                    end
                    if thre > 0
                        if isempty(pinfo.crr.(crrnormfield){spikeind(j)}{ik}) || (pinfo.crr.(crrnormfield){spikeind(j)}{ik} < thre)
                            normsel(ik) = 0;
                        end
                    end
                end
                sesind = find( (sessel == 1) & (normsel == 1) );
                if numel(sesind)>= 1
                   crrnow = zeros(1, nbin); 
                   for (kk = 1:numel(sesind)) %%%%if multiple matches, take the average
                       crrnow = crrnow + data.crr.(crrfield){spikeind(j)}{sesind(kk)}(binind);
                   end
                   allcrr(j,:) = crrnow/numel(sesind);
                end
            end
            plotColorcrr(xbin, allcrr, nbin, tag, clname, sessName, sesstype, keyword, keytype, plottype);
       end 
    end
end
end
end

function plotaveragecrosscrr(pinfo, data, grpname, grpind, sesstype, tag, plottype)
%%%determine normalization factor threhold if crr is in normcount mode
thre = 0; ok = 1;
if strcmp(pinfo.parm.crrmode{grpind{1}(1)}, 'normcount') || strcmp(pinfo.parm.crrmode{grpind{1}(1)}, 'basenorm')
    III=inputdlg('Threshold for normalization factor (0 for all)', 'Pair selection', 1, {'0'}, 'on');
    if (~isempty(III))
        [thre, ok] = str2num(III{1});
    else
        ok = 0;
    end
end
if ~ok
    disp('-------------> Invalid threshold; aborted.')
else

if (strcmp(sesstype, 'session'))
    snfield = 'sessionname';  typefield = 'sessType'; crrfield = 'sesscrr'; crrnormfield = 'sessCrrNorm'; 
    istr = 'Session '; 
elseif (strcmp(sesstype, 'ev'))
    snfield = 'eventname'; typefield = 'eventtype'; crrfield = 'evcrr'; crrnormfield = 'evCrrNorm';
    istr = 'Event ';
end
keyword = []; keytype = []; tit = [istr ' selection (empty for averaging all ', istr, 's']; TT = [istr ' type:'];
III=inputdlg({'Keyword'; TT}, tit, 2, {'Sleep1'; 'ripple'}, 'on');
if (~isempty(III))
   keyword = III{1}; keytype = III{2}; ook = 1;
else
   ook = 0;
end
if ook
for (tt = 1:numel(grpind))
    spikeind = grpind{tt}; nspike = numel(spikeind); clname = grpname{tt};
    %%%%%determine common sessions and common timebin
%     if (strcmp(sesstype, 'session'))
%        if isfield(pinfo.parm, 'sesskeyword') keyword = pinfo.parm.sesskeyword{spikeind(1)}; end
%        if isfield(pinfo.parm, 'sesskeytype') keytype = pinfo.parm.sesskeytype{spikeind(1)}; end
%     elseif (strcmp(sesstype, 'ev'))
%        if isfield(pinfo.parm, 'evtkeyword') keyword = pinfo.parm.evtkeyword{spikeind(1)}; end
%        if isfield(pinfo.parm, 'evtkeytype') keytype = pinfo.parm.evtkeytype{spikeind(1)}; end
%     end

    if isempty(keyword) && isempty(keytype) %%%if no key word/type specified
       sessnow = pinfo.general.(snfield){spikeind(1)}; 
       xbin = data.crr.timebin{spikeind(1)};
       for (i = 2:numel(spikeind))
           sessnow = intersect(sessnow, pinfo.general.(snfield){spikeind(i)});
           xbin = intersect(xbin, data.crr.timebin{spikeind(i)});
       end
       if (~isempty(sessnow)) && (~isempty(xbin))
        for (i = 1:numel(sessnow))
            sessName = sessnow{i}; nbin = numel(xbin); allcrr = NaN*ones(nspike, nbin); allnorm = NaN*ones(nspike,1);
            for (j = 1:nspike)
                [~,binind] = intersect(data.crr.timebin{spikeind(j)}, xbin);
                nses = find(strcmp(pinfo.general.(snfield){spikeind(j)}, sessName));
                allcrr(j,:) = data.crr.(crrfield){spikeind(j)}{nses}(binind);
                allnorm(j) = pinfo.crr.(crrnormfield){spikeind(j)}{nses};
            end
            if (thre > 0) iij = find(allnorm >= thre); allcrr = allcrr(iij,:); end 
            plotavgcrr(xbin, allcrr, nbin, tag, clname, sessName, sesstype, keyword, keytype, plottype);
            
%             crr = NaN*ones(1, nbin); err = NaN*ones(1, nbin); nn = NaN*ones(1, nbin);
%             for (j = 1:nbin)
%                 iii = find(~isnan(allcrr(:,j))); crr(j) = mean(allcrr(iii,j)); 
%                 nn(j) = numel(iii); 
%                 if (nn(j) > 0) err(j) = std(allcrr(iii,j))/sqrt(nn(j)); end
%             end
%             mn = mean(nn); sn = std(nn)/sqrt(nbin);
%             hg = figure('Name', strcat(tag, '---', clname)); 
%             hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
%             xlabel ('time lag (s)'); ylabel('Mean correlation'); 
%             line(xbin, crr, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]);
%             Drawerrorupdown(xbin, crr, crr+err, crr-err, hax, [1 0 0]);
%             text('Interpreter', 'none', 'Parent', hax, 'String', sessName, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.94]);
%             str = ['Pair # (mean, se): ', num2str(mn), ', ', num2str(sn)];
%             text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.91]);
        end
       else
        disp('-------------> warning: no common sessions or common time lags found');
       end
    else %%%%%if keyword or keytype are specified
       xbin = data.crr.timebin{spikeind(1)}; sessName = [];
       for (i = 2:numel(spikeind))
            xbin = intersect(xbin, data.crr.timebin{spikeind(i)});
       end
       if (~isempty(xbin))
            nbin = numel(xbin); allcrr = NaN*ones(nspike, nbin); 
            for (j = 1:nspike)
                [~,binind] = intersect(data.crr.timebin{spikeind(j)}, xbin);
                nses = numel(pinfo.general.(snfield){spikeind(j)}); sessel = ones(1, nses); normsel = ones(1, nses);
                for (ik = 1:nses)
                    if ~isempty(keyword)
                       if isempty(strfind(lower(pinfo.general.(snfield){spikeind(j)}{ik}), lower(keyword))) 
                           sessel(ik) = 0;
                       end
                    end
                    if ~isempty(keytype)
                       if ~strcmpi(pinfo.parm.(typefield){spikeind(j)}{ik}, keytype) sessel(ik) = 0; end
                    end
                    if thre > 0
                        if isempty(pinfo.crr.(crrnormfield){spikeind(j)}{ik}) || (pinfo.crr.(crrnormfield){spikeind(j)}{ik} < thre)
                            normsel(ik) = 0;
                        end
                    end
                end
                sesind = find( (sessel == 1) & (normsel == 1) );
                if numel(sesind)>= 1
                   crrnow = zeros(1, nbin); 
                   for (kk = 1:numel(sesind)) %%%%if multiple matches, take the average
                       crrnow = crrnow + data.crr.(crrfield){spikeind(j)}{sesind(kk)}(binind);
                   end
                   allcrr(j,:) = crrnow/numel(sesind);
                end
            end
            plotavgcrr(xbin, allcrr, nbin, tag, clname, sessName, sesstype, keyword, keytype, plottype);
            
%             crr = NaN*ones(1, nbin); err = NaN*ones(1, nbin); nn = NaN*ones(1, nbin);
%             for (j = 1:nbin)
%                 iii = find(~isnan(allcrr(:,j))); crr(j) = mean(allcrr(iii,j)); 
%                 nn(j) = numel(iii); 
%                 if (nn(j) > 0) err(j) = std(allcrr(iii,j))/sqrt(nn(j)); end
%             end
%             mn = mean(nn); sn = std(nn)/sqrt(nbin);
%             hg = figure('Name', strcat(tag, '---', clname)); 
%             hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
%             xlabel ('time lag (s)'); ylabel('Mean correlation'); 
%             line(xbin, crr, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]);
%             Drawerrorupdown(xbin, crr, crr+err, crr-err, hax, [1 0 0]);
%             text('Interpreter', 'none', 'Parent', hax, 'String', strcat(sesstype, '_', keyword, '_', keytype), 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.94]);
%             str = ['Pair # (mean, se): ', num2str(mn), ', ', num2str(sn)];
%             text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.91]);
       end 
    end
end
end
end

function plotColorcrr(xbin, allcrr, nbin, tag, clname, sessName, sesstype, keyword, keytype, plottype)
%%%Plot individual crr using color: 
%%%Plottype: centralpeakorder - arrange all zeros at the bottom, also order according to central peak values
if strcmp(plottype, 'centralpeakorder') %%reorder allcrr
    allcrr(isnan(allcrr)) = 0; %%replacing all NaNs by 0
    [npair,~] = size(allcrr); cind = find(xbin == 0);
    abssum = ones(1, npair); peakvalue = zeros(1, npair);
    for (i = 1:npair)
        crrnow = allcrr(i,:); peakvalue(i) = crrnow(cind); abssum(i) = sum(abs(crrnow)); 
        %abssum(i) = sum(abs(crrnow(~isnan(crrnow)))); 
    end
    zeropairind = find(abssum == 0); otherpairind = setdiff([1:npair], zeropairind);
    otherpairpeakvalue = peakvalue(otherpairind); [~, iii] = sort(otherpairpeakvalue); otherpairind = otherpairind(iii);
    newpairind = [zeropairind otherpairind];
    allcrr= allcrr(newpairind,:);
end
%%%plot data
hg = figure('Name', strcat(tag, '---', clname)); 
hax = axes('Parent', hg, 'YLim', [1 npair], 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
xlabel ('time lag (s)'); ylabel('Pair number'); 
hp = pcolor(hax, xbin, 1:npair, allcrr);
set(hp, 'EdgeColor', 'none');
colormap(jet);  
hc = colorbar('vert', 'peer', hax);
%%%plot info
if isempty(keyword) || isempty(keytype) %%%if no key word/type specified
   text('Interpreter', 'none', 'Parent', hax, 'String', sessName, 'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.02 0.94]);
else
   text('Interpreter', 'none', 'Parent', hax, 'String', strcat(sesstype, '_', keyword, '_', keytype), 'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.02 0.94]);
end
text('Interpreter', 'none', 'Parent', hax, 'String', plottype, 'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.02 0.88]);

function plotavgcrr(xbin, allcrr, nbin, tag, clname, sessName, sesstype, keyword, keytype, plottype)
crr = NaN*ones(1, nbin); upcrr = NaN*ones(1, nbin); downcrr = NaN*ones(1, nbin); nn = NaN*ones(1, nbin);
for (j = 1:nbin)
    iii = find(~isnan(allcrr(:,j))); nn(j) = numel(iii);
    if nn(j) > 0
       if strcmp(plottype, 'mean')
          crr(j) = mean(allcrr(iii,j)); err = std(allcrr(iii,j))/sqrt(nn(j)); upcrr(j) = crr(j)+err; downcrr(j) = crr(j)-err; 
       elseif strcmp(plottype, 'median')
          crr(j) = median(allcrr(iii,j)); upcrr(j) = prctile(allcrr(iii,j), 75); downcrr(j) = prctile(allcrr(iii,j), 25); 
       end
    end
end
%%%plot data
hg = figure('Name', strcat(tag, '---', clname)); 
hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
xlabel ('time lag (s)'); ylabel('Average correlation'); 
if strcmp(plottype, 'mean')
   Drawerrorupdown(xbin, crr, upcrr, downcrr, hax, [1 0 0]);
elseif strcmp(plottype, 'median')
   hcurve = area(hax, xbin, [downcrr' (upcrr-downcrr)']); set(hcurve(1), 'FaceColor', [1 1 1]); set(hcurve(2), 'FaceColor', [0.8 0.8 0.8]);  
end
line(xbin, crr, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]);
%%%plot info
if isempty(keyword) || isempty(keytype) %%%if no key word/type specified
   text('Interpreter', 'none', 'Parent', hax, 'String', sessName, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.94]);
else
   text('Interpreter', 'none', 'Parent', hax, 'String', strcat(sesstype, '_', keyword, '_', keytype), 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.94]);
end
mn = mean(nn); sn = std(nn)/sqrt(nbin);
str = ['Pair # (mean, se): ', num2str(mn), ', ', num2str(sn)];
text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.91]);
if strcmp(plottype, 'mean')
    str = ['Plot type: mean +- se'];
elseif strcmp(plottype, 'median')
    str = ['Plot type: median [25% 75%] range'];
end
text('Interpreter', 'none', 'Parent', hax, 'String', str, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.88]);

function plotcrosscrr(pinfo, data, spikeind, sesstype, tag)
if (numel(spikeind)>10)
    disp('--------> Too many cells selected; Only plot fields for the first 10 cells');
    spikeind = spikeind(1:10);
end
for (tt = 1:numel(spikeind))
     i = spikeind(tt); %%for the current cell pair, get crr and timebin
     clname = pinfo.general.clname{i}; xbin = data.crr.timebin{i};
     if strcmp(sesstype, 'session')
         sessName = pinfo.general.sessionname{i}; crr = data.crr.sesscrr{i}; 
         P1time = pinfo.crr.sess1stPtime{i}; P1crr = pinfo.crr.sess1stPcrr{i};
         P2time = pinfo.crr.sess2ndPtime{i}; P2crr = pinfo.crr.sess2ndPcrr{i};
     elseif strcmp(sesstype, 'ev')
         sessName = pinfo.general.eventname{i}; crr = data.crr.evcrr{i};
         P1time = pinfo.crr.ev1stPtime{i}; P1crr = pinfo.crr.ev1stPcrr{i};
         P2time = pinfo.crr.ev2ndPtime{i}; P2crr = pinfo.crr.ev2ndPcrr{i};
     end
     for (j = 1:numel(sessName))
         %%%plot crrs and timebin
         hg = figure('Name', strcat(tag, '---', clname)); 
         hax = axes('Parent', hg, 'NextPlot', 'add'); %, 'XLim', [min(xbin) max(xbin)], 'YLim', [-0.5 1]);
         xlabel ('time lag (s)'); ylabel('Correlation'); 
         if ~isempty(crr{j})
            line(xbin, crr{j}, 'Parent', hax, 'LineWidth', 2, 'Color', [0 0 0]);
            %%%find and plot peaks identified
            line(P1time{j}, P1crr{j}, 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30,...
                  'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
            line(P2time{j}, P2crr{j}, 'Parent', hax, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30,...
                  'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]);
         end
         text('Interpreter', 'none', 'Parent', hax, 'String', sessName{j}, 'Color', [0 0 0], 'Units', 'normalized', 'Position', [0.02 0.94]);

    end
end



