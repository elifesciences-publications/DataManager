function [pinfo, data] = DataManager_FindCrrAllSpikePairs(pinfo, spikepinfo, spikedata, spikeselection)

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

% %%%%get EEG event data and append to the spike general data structure
% subf = fieldnames(spikepinfo.general);
% if (~isempty(EEGfile))
% [pp, nn, ee] = fileparts(EEGfile);
% if (strcmp(ee, '.eegdb'))
%     disp('---------> add EEG events');
%     [eegev, nev] = geteegevent(EEGfile);
% if (nev > 0)
%    for (i = 1:numel(subf))
%             [A,B] = size(spikepinfo.general.(subf{i}));
%             if (A ==1)
%                 if (isfield(eegev, subf{i}))
%                     spikepinfo.general.(subf{i}) = [spikepinfo.general.(subf{i}) eegev.(subf{i})];
%                 else
%                      ttt = cell(1, nev); spikepinfo.general.(subf{i}) = [spikepinfo.general.(subf{i}) ttt];
%                 end
%             elseif (B ==1)
%                 if (isfield(eegev, subf{i}))
%                    spikepinfo.general.(subf{i}) = [spikepinfo.general.(subf{i}); eegev.(subf{i})];
%                 else
%                      ttt = cell(nev, 1); spikepinfo.general.(subf{i}) = [spikepinfo.general.(subf{i}); ttt];
%                 end
%             end
%    end
%    spikedata.spike.spiketime = [spikedata.spike.spiketime eegev.data];
%    spikedata.events.eventtimes = [spikedata.events.eventtimes eegev.eventtimes];
% end
% end
% end
% % for (i = 1:numel(spikepinfo.general.clname))
% %      disp(spikepinfo.general.clname{i});
% %      disp(spikepinfo.general.genotype{i});
% %      disp(spikepinfo.general.datedir{i});
% % end

data.spike.spiketime = spikedata.spike.spiketime; %direct transfer of spike times
data.events.eventtimes = spikedata.events.eventtimes;
subf = fieldnames(spikepinfo.general);
npair = 0; %start count of spike pairs found -- count day by day
datedir = unique(spikepinfo.general.datedir); % for unique datedir
ndatedir = numel(datedir);
disp(['---------> number of dates ', num2str(ndatedir)]);
for (i = 1:ndatedir)
    iii = find( strcmp(spikepinfo.general.datedir, datedir{i}) );
    for (m = 1:numel(iii))
        for (n = m:numel(iii))
            s1ind = iii(m); s2ind = iii(n); 
            if (spikeselection(s1ind) == 1) && (spikeselection(s2ind) == 1)
            ses1 = spikepinfo.general.sessionname{s1ind}; ses2 = spikepinfo.general.sessionname{s2ind};
            [sesnow, sesin1, ~] = intersect(ses1, ses2);
            ev1 = spikepinfo.general.eventname{s1ind}; ev2 = spikepinfo.general.eventname{s2ind};
            [evnow, evin1, ~] = intersect(ev1, ev2);
            if (~isempty(sesnow))
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

% function [eegev, nev] = geteegevent(EEGfile)
% nev = 0; eegev = [];
%     S = load(EEGfile, '-mat'); %load the file
%     eeg = S.eeg; eegdata = S.eegdata; S = [];
%     if (isfield(eeg, 'ripple'))
%         [eegevnow, nevnow] = getrippleev(eeg, eegdata, 'ripple');
%         nev = nev + nevnow; eegev = combineeegev(eegev, eegevnow);
%     end
%     if (isfield(eeg, 'theta'))
%         [eegevnow, nevnow] = getrippleev(eeg, eegdata, 'theta');
%         nev = nev + nevnow; eegev = combineeegev(eegev, eegevnow);
%     end
%     if (isfield(eeg, 'spindle'))
%         [eegevnow, nevnow] = getrippleev(eeg, eegdata, 'spindle');
%         nev = nev + nevnow; eegev = combineeegev(eegev, eegevnow);
%     end
% 
% function eegev = combineeegev(eegev, eegevnow)
% if (~isempty(eegevnow))
%     if (isempty(eegev))
%         eegev = eegevnow;
%     else
%         nev = numel(eegevnow.clname);
%         if (nev >0)
%             subf = fieldnames(eegev);
%             for (i = 1:numel(subf))
%                 [A,B] = size(eegev.(subf{i}));
%                 if (A ==1)
%                     if (isfield(eegevnow, subf{i}))
%                         eegev.(subf{i}) = [eegev.(subf{i}) eegevnow.(subf{i})];
%                     else
%                         ttt = cell(1, nev); eegev.(subf{i}) = [eegev.(subf{i}) ttt];
%                     end
%                 elseif (B ==1)
%                     if (isfield(eegevnow, subf{i}))
%                         eegev.(subf{i}) = [eegev.(subf{i}); eegevnow.(subf{i})];
%                     else
%                         ttt = cell(nev,1); eegev.(subf{i}) = [eegev.(subf{i}); ttt];
%                     end
%                 end
%             end
%         end
%    end
% end
% 
% function [ev, nev] = getrippleev(eeg, eegdata, mark)
% nev = 0; ev = []; neeg = numel(eeg.general.recarea); checkstr = cell(1, neeg);
% for (i = 1:neeg)
%     checkstr{i} = strcat(eeg.general.finaldir{i}, filesep, eeg.general.eegTT{i});
% end
% iii = find(strcmp(eeg.parm.band, mark));
% if (~isempty(iii))
%     subf = fieldnames(eeg.general);
%     for (i = 1:numel(subf))
%         eeg.general.(subf{i}) = eeg.general.(subf{i})(iii);
%     end
%     subf = fieldnames(eeg.(mark));
%     for (i = 1:numel(subf))
%         eeg.(mark).(subf{i}) = eeg.(mark).(subf{i})(iii);
%     end
%     eegdata.event.eventtimes = eegdata.event.eventtimes(iii);
%     checkstr = checkstr(iii);
% end
% timestr = findtimestr(mark); %%%these are the variables to export
% unistr = unique(checkstr);
% for (i = 1:numel(unistr)) %%%this is for combining different EEG event files (during different sessions) together to make a spike file for the whole day
%     iii = find(strcmp(checkstr, unistr{i})); 
%     iscomputed = zeros(1, numel(iii)); 
%     for (j = 1:numel(iii))
%         if (isfield(eeg.(mark), 'sessNum'))
%            if ~isempty(eeg.(mark).sessNum{iii(j)}) iscomputed(j) = 1; end
%         end
%     end
%     iii = iii(iscomputed==1); %%%now this contains all sessions within the same day in which the band-specific events are defined 
%     if (~isempty(iii))
%     for (k = 1:numel(timestr))
%         nev = nev + 1; ev.data{nev} = []; ev.eventtimes{nev} = cell(0,1); 
%         ev.sessionname{nev} = cell(0,1); ev.eventname{nev} = cell(0,1); 
%         ev.sessionstartT{nev} = []; ev.sessionendT{nev} = []; ev.sessionlength{nev} = [];
%         for (j = 1:numel(iii))
%             [A,B] = size(eeg.(mark).(timestr{k}){iii(j)});
%             if (A ==1)
%                ev.data{nev} = [ev.data{nev} eeg.(mark).(timestr{k}){iii(j)}]; 
%             elseif (B == 1)
%                ev.data{nev} = [ev.data{nev}; eeg.(mark).(timestr{k}){iii(j)}]; 
%             end
%             ev.sessionname{nev} = [ev.sessionname{nev}; eeg.general.sessname{iii(j)}];
%             ev.sessionstartT{nev} = [ev.sessionstartT{nev}; eeg.general.sessstartT{iii(j)}];
%             ev.sessionendT{nev} = [ev.sessionendT{nev}; eeg.general.sessendT{iii(j)}];
%             ev.sessionlength{nev} = [ev.sessionlength{nev}; eeg.general.sesslength{iii(j)}];
%             ev.eventname{nev} = [ev.eventname{nev} eeg.general.eventname{iii(j)}];
%             ev.eventtimes{nev} = [ev.eventtimes{nev} eegdata.event.eventtimes{iii(j)}];
%         end
%         [datename, tok] = strtok(eeg.general.datedir{iii(1)}, '_'); 
%         ev.clname{nev} = strcat(datename, '_', eeg.general.eegTT{iii(1)}, '_', mark, '_', timestr{k}); %{str};
%         ev.TTname{nev} = eeg.general.eegTT{iii(1)};
%         ev.parmfile{nev} = unistr{i}; ev.wavefile{nev} = [];
%         ev.gain{nev}=[]; 
%         subf = fieldnames(eeg.general);
%         for (j = 1:numel(subf))
%             if (~(isfield(ev, subf{j}))) || (numel(ev.(subf{j})) < nev) %%%if field not defined yet as above
%                ev.(subf{j}){nev} = eeg.general.(subf{j}){iii(1)};
%             end
%         end
%     end
%     end
% end
% 
% function timestr = findtimestr(mark)
% timestr = [];
% if (strcmp(mark, 'ripple'))
%     timestr = {'sessStartT'; 'sessEndT'; 'sessPeakT'};
% elseif (strcmp(mark, 'theta'))
%     timestr = {'maxTime'; 'minTime'};
% elseif (strcmp(mark, 'spindle'))
%     timestr = {'maxTime'; 'minTime'; 'sessStartT'; 'sessEndT'; 'sessPeakT'};
% end
