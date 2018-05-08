function [pinfo, data] = DataManager_FindCrrAllSpikePairs(pinfo, spikepinfo, spikedata, cellind, selectiontype)

%Find all spike pairs: only cells with its cellind in spikepinfo assigned 1 are selected 
%Fields assigned here:
pinfo.general.recarea = []; pinfo.general.finaldir = []; pinfo.general.datedir = []; pinfo.general.animalname = [];  %{str}
%pinfo.general.genotype = []; %{str} %pinfo.general.age = []; %{str}
pinfo.general.clname = []; pinfo.general.TTname = []; pinfo.general.parmfile = []; pinfo.general.wavefile = []; %{str};
       data.spike.spiketime = []; data.crr.pairind = []; data.events.eventtimes = []; %{{ev.start, ev.ent, ev.ref, ev.marker}}
pinfo.general.sessionname = []; % {{session1; session2;}};
pinfo.general.sessionstartT = []; pinfo.general.sessionendT = []; pinfo.general.sessionlength = []; % {[l1 l2]};
pinfo.general.eventname = []; %{{str1; str2; ...}}
pinfo.general.crrtype = [];

data.spike.spiketime = spikedata.spike.spiketime; %direct transfer of spike times
data.events.eventtimes = spikedata.events.eventtimes;
subf = fieldnames(spikepinfo.general);
npair = 0; %start count of spike pairs found -- count day by day
datedir = unique(spikepinfo.general.datedir); % for unique datedir
ndatedir = numel(datedir); allN = 0; dsrate = 1;
disp(['---------> number of dates ', num2str(ndatedir)]);
for (i = 1:ndatedir)
    iii = find( strcmp(spikepinfo.general.datedir, datedir{i}) );
    for (m = 1:numel(iii))
        for (n = m:numel(iii))
            s1ind = iii(m); s2ind = iii(n); 
            if isingroups(cellind, s1ind, s2ind, selectiontype)
            allN = allN + 1;
            ses1 = spikepinfo.general.sessionname{s1ind}; ses2 = spikepinfo.general.sessionname{s2ind};
            [sesnow, sesin1, ~] = intersect(ses1, ses2);
            ev1 = spikepinfo.general.eventname{s1ind}; ev2 = spikepinfo.general.eventname{s2ind};
            [evnow, evin1, ~] = intersect(ev1, ev2);
            if (~isempty(sesnow)) && (mod(allN, dsrate)==0)
                npair = npair + 1; data.crr.cellind{npair} = [s1ind s2ind];
                if (m==n) 
                    pinfo.general.crrtype{npair} = 'auto';
                else
                    pinfo.general.crrtype{npair} = 'cross';
                end
                for (j = 1:numel(subf))
                    if (strcmp(subf{j}, 'datedir'))
                        pinfo.general.(subf{j}){npair} = datedir{i};
                    elseif (strcmp(subf{j}, 'wavefile')) | (strcmp(subf{j}, 'gain'))
                    elseif (strcmp(subf{j}, 'sessionname'))
                        pinfo.general.(subf{j}){npair} = ses1; %sesnow;
                    elseif (strcmp(subf{j}, 'sessionstartT'))
                        %disp(spikepinfo.general.sessionstartT{s1ind})
                        %sesin1
                        pinfo.general.(subf{j}){npair} = spikepinfo.general.sessionstartT{s1ind}; %(sesin1);
                    elseif (strcmp(subf{j}, 'sessionendT'))
                        pinfo.general.(subf{j}){npair} = spikepinfo.general.sessionendT{s1ind}; %(sesin1);
                    elseif (strcmp(subf{j}, 'sessionlength'))
                        pinfo.general.(subf{j}){npair} = spikepinfo.general.sessionlength{s1ind}; %(sesin1);
                    elseif (strcmp(subf{j}, 'eventname'))
                        pinfo.general.(subf{j}){npair} = ev1; %evnow;
                    elseif (~strncmpi(subf{j}, 'clust', 4))
       pinfo.general.(subf{j}){npair} = strcat(spikepinfo.general.(subf{j}){s1ind}, '__', spikepinfo.general.(subf{j}){s2ind});
                    elseif (strcmp(subf{j}, 'clustMaxFitErr'))
                        pinfo.general.(subf{j}){npair} = [spikepinfo.general.(subf{j}){s1ind} spikepinfo.general.(subf{j}){s2ind}];
                    elseif (strcmp(subf{j}, 'clustMaxCutoffI'))
                        pinfo.general.(subf{j}){npair} = [spikepinfo.general.(subf{j}){s1ind} spikepinfo.general.(subf{j}){s2ind}];
                    end
                end
            end
            %disp(pinfo.general.clname{npair});
            end
        end
    end
end
disp(['---------> number of cell pairs ', num2str(npair)]);

function ok = isingroups(cellind, s1ind, s2ind, selectiontype) 
ok = 0;
if strcmp(selectiontype, 'Allgroups')
    cellok = cell2mat(cellind);
    if (~isempty(find(cellok == s1ind))) && (~isempty(find(cellok == s2ind)))
        ok = 1;
    end
elseif strcmp(selectiontype, 'Withingroups')
    for (i = 1:numel(cellind))
        cellok = cellind{i};
        if (~isempty(find(cellok == s1ind))) && (~isempty(find(cellok == s2ind)))
           ok = 1; break
        end
    end
elseif strcmp(selectiontype, 'Crossgroups')
    for (i = 1:numel(cellind))
        cellok1 = cellind{i};
        for j = i+1:numel(cellind)
            cellok2 = cellind{j};
            if ((~isempty(find(cellok1 == s1ind))) && (~isempty(find(cellok2 == s2ind)))) || ((~isempty(find(cellok1 == s2ind))) && (~isempty(find(cellok2 == s1ind))))
                ok = 1; break
            end
        end
    end
end


