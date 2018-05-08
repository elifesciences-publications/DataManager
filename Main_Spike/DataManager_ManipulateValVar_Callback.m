function DataManager_ManipulateValVar_Callback
%%Manipulate variables and values of a database

hf = gcbf; dbtype = getappdata(hf, 'dbtype'); celllistvar = getappdata(hf, 'celllistvar'); celllisttitle = getappdata(hf, 'celllisttitle');
%%%%%%%added options for all types of databases
pinfostr = 'pinfo'; datastr = 'data'; catok = 'general'; varok = 'parmfile';
switch dbtype
    case '.eegdb'
        pinfostr = 'eeg'; datastr = 'eegdata'; catok = 'general'; varok = 'eegfile'; %%also determine a variable (one in pinfo.general) to identify cells
    case '.behavdb'
        pinfostr = 'behav'; datastr = 'bhdata'; catok = 'general'; varok = 'sessID'; 
    case '.seqdb'
        catok = 'general'; varok = 'tmpID';     
end
pinfo = getappdata(hf, pinfostr); data = getappdata(hf, datastr); 
tagmark = get(gcbo, 'Tag');
hgroup = getappdata(hf, 'hgroup');hfield = getappdata(hf, 'hfield');
ok = 1; plotparm = getappdata(hf, 'plotparm'); vv = plotparm.showdetail; %if vv=1, shoe details of computation
cc = plotparm.compute; %if cc=0, assign parameters for computation only; if 1, do the real computation
ow = plotparm.overwrite; %if ow=1, plot into the current database
[writefilename, okk] = getoutputfile(hf, ow, dbtype);

%%%group selected
groupselection = getappdata(hgroup, 'selection'); 
%%%field selected
oldfield = getappdata(hf, 'hfield'); for (tth = 1:numel(oldfield)) oldFselection{tth} = getappdata(oldfield(tth), 'selection'); end

%get selected cellind
cellind = []; grpind = find(groupselection == 1); 
%also for changeparm, need to use current selected spikeind for versatility
hspike = getappdata(hf, 'hspike'); spikeselection = getappdata(hspike, 'selection');
currentind = find(spikeselection == 1);
for (kk = 1:numel(grpind)) cellind = union(cellind, data.grouplist.groupindex{grpind(kk)}); end
%%%now do the manipulations
if ((~isempty(cellind))||(~isempty(currentind))) && okk
    if (strcmp(tagmark, 'removecells'))
        [pinfo,data, ok] = removecells(pinfo,data, grpind, catok, varok); %%%remove cells from a database
    elseif (strcmp(tagmark, 'removevars'))
        [pinfo,data, ok] = removevars(pinfo,data, hfield); %%%remove variables for all cells in a database
    elseif (strcmp(tagmark, 'renamevars'))
        [pinfo,data, ok] = renamevars(pinfo,data, hfield); %%%remove variables for all cells in a database
    elseif (strcmp(tagmark, 'copyvars'))
        [pinfo,data, ok] = copyvars(pinfo,data, hfield); %%%copy variables for all cells in a database
    elseif (strcmp(tagmark, 'squeezevals'))
        [pinfo,data, ok] = squeezevals(pinfo,data, cellind, hfield); %%%change (not re-compute) values of selected variables for selected cell groups
    elseif (strcmp(tagmark, 'crosssqueezevals'))
        [pinfo,data, ok] = crosssqueezevals(pinfo,data, cellind, hfield);
    elseif (strcmp(tagmark, 'changeparm'))
        [pinfo,data, ok] = changeparm(pinfo,data, currentind, hfield);
    elseif (strcmp(tagmark, 'createNewVar'))
        [pinfo,data, ok] = createnewvar(pinfo,data, cellind, hfield, catok, varok);
    elseif (strcmp(tagmark, 'reassignvals'))
        [pinfo,data, ok] = reassignvals(pinfo,data, cellind, hfield);
    end
    %%%%%%%%save and plot the new database
    if (~isempty(pinfo.(catok).(varok)))
        if (ok)
            if (strcmp(dbtype, '.eegdb'))
                eeg = pinfo; eegdata = data;  if (ow ==0) save(writefilename, 'eeg', 'eegdata'); end
                plotparm.linkspike = 0; plotparm.linkeeg = 1; plotparm.linkbehav = 0;
            elseif (strcmp(dbtype, '.behavdb'))
                behav = pinfo; bhdata = data; if (ow==0) save(writefilename, 'behav', 'bhdata'); end
                plotparm.linkspike = 0; plotparm.linkeeg = 0; plotparm.linkbehav = 1;
            elseif (strcmp(dbtype, '.spikedb'))
                if (ow ==0) save(writefilename, 'pinfo', 'data'); end
                plotparm.linkspike = 1; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
            else
                if (ow ==0) save(writefilename, 'pinfo', 'data'); end
                plotparm.linkspike = 0; plotparm.linkeeg = 0; plotparm.linkbehav = 0;
            end
            if (ow)
               iii = get(hf, 'Children'); delete(iii(~strcmp(get(iii, 'Type'), 'uimenu'))); 
               fname = get(hf, 'Name'); ftt = strfind(fname, '__'); 
               newname = fname(1:ftt-1); set(hf, 'Name', newname); 
               hmain = hf; 
            else
               hmain = DataManager_DataManager_Callback; 
            end
            setappdata(hmain,'plotparm', plotparm);
            DataManager_PlotSpikeDatabase(hmain, pinfo, data, celllisttitle, celllistvar, dbtype);
            set(hmain, 'Name', strcat(get(hmain, 'Name'), '__', writefilename));
            setappdata(hmain, pinfostr, pinfo); setappdata(hmain, datastr, data); pinfo = []; data = [];
            %%%%%%%%%%%%re-set all the selections and update display -- can only do this if overwrite
            %%%%%group selections
            if ow
            hgroup = getappdata(hmain, 'hgroup'); 
            hgselection = getappdata(hgroup, 'selection'); oldselection = groupselection; htext = getappdata(hgroup, 'htext');
            if (numel(hgselection) == numel(oldselection))
                iii = find(hgselection ~= oldselection);
                setappdata(hgroup, 'selection', oldselection);
                for (jjj = 1:numel(iii))
                    if (oldselection(iii(jjj)) == 0) & (iii(jjj)<=numel(htext))
                        set(htext(iii(jjj)), 'Color', [0 0 0]); 
                    elseif (iii(jjj)<=numel(htext))
                        set(htext(iii(jjj)), 'Color', [1 0 0]); 
                    end
                end
            end
            %%%%%%field selections
            hfield = getappdata(hmain, 'hfield'); 
            if (numel(hfield) == numel(oldfield))
               for (tth = 1:numel(hfield))
                    hgselection = getappdata(hfield(tth), 'selection'); htext = getappdata(hfield(tth), 'htext');
                    if (numel(hgselection) == numel(oldFselection{tth}))
                        iii = find(hgselection ~= oldFselection{tth});
                        setappdata(hfield(tth), 'selection', oldFselection{tth});
                        for (jjj = 1:numel(iii))
                             if (oldFselection{tth}(iii(jjj)) == 0)&(iii(jjj)<=numel(htext)) 
                                 set(htext(iii(jjj)), 'Color', [0 0 0]); 
                             elseif (iii(jjj)<=numel(htext))
                                 set(htext(iii(jjj)), 'Color', [1 0 0]); 
                             end
                        end
                    end
               end  
            end
            setappdata(hmain, 'hgroup', hgroup); setappdata(hmain, 'hfield', hfield);
            DataManager_UpdateDisplay(hmain);
            end
        end
    else
        disp('-------------> no cells in the database!');
    end
else
    disp('--------------> no groups/cells selected, groups do not contain any cells, or action cancelled');
end
%disp('**********************');

function [pinfo,data, ok] = changeparm(pinfo,data, cellind, hfield)
ok = 1;
input = questdlg('Are you sure to change the selected parameter for the selected group of cells (some variable values will be invalid anymore)?', 'Change parameters', 'No');
if (strcmp(input, 'Yes'))
    fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle);
    for (i = 1:nfield) 
        if (strcmp(fieldtitle{i}, 'parm'))
            subfieldnow = fieldnames(pinfo.(fieldtitle{i}));
            fieldselection = getappdata(hfield(i), 'selection'); fieldind = find(fieldselection == 1);
            seleparmname = subfieldnow(fieldind); nsel = numel(fieldind); 
            if (nsel > 0)
            for (tt = 1:nsel)
                def{tt} = []; typ{tt} = [];  pp{tt} = ['Enter ', seleparmname{tt}, ':']; 
                if iscell(pinfo.parm.(seleparmname{tt}))
                    if iscell(pinfo.parm.(seleparmname{tt}){cellind(1)})
                       typ{tt} = 'cell'; def{tt} = rearangecell(pinfo.parm.(seleparmname{tt}){cellind(1)});
                    elseif ischar(pinfo.parm.(seleparmname{tt}){cellind(1)})
                       typ{tt} = 'str'; def{tt} = pinfo.parm.(seleparmname{tt}){cellind(1)};
                    else
                       typ{tt} = 'numM'; def{tt} = rearrangeMnum(pinfo.parm.(seleparmname{tt}){cellind(1)});
                    end
                else
                    typ{tt} = 'numS'; def{tt} = num2str(pinfo.parm.(seleparmname{tt})(cellind(1))); 
                end
            end
            III=inputdlg(pp, 'Input for selected parmameters', nsel, def, 'on');
            if (~isempty(III))
                for (tt = 1:nsel)
                if (~isempty(typ{tt}) && (~isempty(III{tt}))) 
                    for (tt = 1:nsel)
                        if (strcmp(typ{tt}, 'numS'))
                            pinfo.parm.(seleparmname{tt})(cellind) = ones(numel(cellind),1)*str2num(III{tt});
                        elseif (strcmp(typ{tt}, 'cell'))
                            for (kk = 1:numel(cellind))
                                pinfo.parm.(seleparmname{tt}){cellind(kk)} = arrangebackcell(III{tt});
                            end
                        elseif (strcmp(typ{tt}, 'numM'))
                            for (kk = 1:numel(cellind))
                                pinfo.parm.(seleparmname{tt}){cellind(kk)} = arrangebackMnum(III{tt})';
                            end 
                        elseif (strcmp(typ{tt}, 'str'))
                            for (kk = 1:numel(cellind))
                                pinfo.parm.(seleparmname{tt}){cellind(kk)} = III{tt};
                            end 
                        end
                    end
                    %disp('*****entered the inner loop');
                elseif isempty(III{tt})
                    for (tt = 1:nsel)
                        if (strcmp(typ{tt}, 'numS')) %%%number single
                            pinfo.parm.(seleparmname{tt})(cellind) = NaN*ones(numel(cellind),1);
                        else
                            for (kk = 1:numel(cellind))
                                pinfo.parm.(seleparmname{tt}){cellind(kk)} = [];
                            end 
                        end
                    end
                else
                    disp('---------> unknown parameter type or parameter not assigned'); ok = 0;
                end
                end
            else
                ok = 0; disp('---------> no input or action cancelled'); 
            end
            else
                ok = 0; disp('---------> no parameters selected'); 
            end
        end
    end 
else
    ok = 0;
end

function [pinfo,data, ok] = crosssqueezevals(pinfo,data, cellind, hfield)
ok = 1; worksubfield = fieldnames(pinfo.work);
[sel,ok] = listdlg('ListString', worksubfield, 'SelectionMode', 'single', 'PromptString', 'Select a lead variable to squeeze');
if (ok)
    squfield = worksubfield{sel}; coind = setdiff( 1:numel(worksubfield), sel); cofield = worksubfield(coind);
    [sel,ok] = listdlg('ListString', cofield, 'SelectionMode', 'multiple', 'PromptString', 'Select one or more variables to co-squeeze');
    if (ok)
        %cosqufield = worksubfield(sel);
        cosqufield = cofield(sel);
        vartype = resolvetypeformultiple(pinfo.work.(squfield)(cellind));
        if (strcmp(vartype, 'cellcell')) || (strncmpi(vartype, 'vectorcell',10))
            [flag, sel] = askvectorkeyword(strcat('work.',squfield), vartype, 'crossvalue');
            for (kkk = 1:numel(cellind))
                 if iscell(pinfo.work.(squfield))
                    valthen = pinfo.work.(squfield){cellind(kkk)}; nval = numel(valthen);
                    [valnow,ind] = squeezenow(valthen, flag, sel, vartype);
                    pinfo.work.(squfield){cellind(kkk)} = valnow;
                 else
                    valthen = pinfo.work.(squfield)(cellind(kkk)); nval = numel(valthen);
                    [valnow,ind] = squeezenow(valthen, flag, sel, vartype);
                    pinfo.work.(squfield)(cellind(kkk)) = valnow;
                 end
                 %if (~isempty(ind))
                     for (i = 1:numel(cosqufield))
                         if (iscell(pinfo.work.(cosqufield{i})))
                             valthen = pinfo.work.(cosqufield{i}){cellind(kkk)};
                             %if (numel(valthen) == nval)
                             if isempty(ind) | (max(ind) <= numel(valthen))
                                 pinfo.work.(cosqufield{i}){cellind(kkk)} = valthen(ind);
                             else
                                 ind = ind(ind<=numel(valthen));
                                 pinfo.work.(cosqufield{i}){cellind(kkk)} = valthen(ind);
                                 disp('----------> warning: co-squeeze indices out of range');
                             end
                         else
                             valthen = pinfo.work.(cosqufield{i})(cellind(kkk));
                             %if (numel(valthen) == nval)
                             if isempty(ind) | (max(ind)<=numel(valthen))    
                                 pinfo.work.(cosqufield{i})(cellind(kkk)) = valthen(ind);
                             else
                                 ind = ind(ind<=numel(valthen));
                                 pinfo.work.(cosqufield{i})(cellind(kkk)) = valthen(ind);
                                 disp('----------> warning: co-squeeze indices out of range');
                             end
                         end
                     end
                 %else
                 %    disp('----------> can not co-squeeze: squeeze mode does not generate an index');
                 %end
            end
        end
    end
end

function [pinfo,data, ok] = createnewvar(pinfo,data, cellind, hfield, catok, varok)
choice = {'Construct from selected work variables'; 'New manual assignment'};
[sel,ok] = listdlg('ListString', choice, 'SelectionMode', 'single', 'PromptString', 'Create a new variable');
if ok
if (sel == 2)
        pp = {'Enter variable name'; 'Enter variable value'}; def = {'var1'; '1'};
        III=inputdlg(pp, 'Create a new variable', 2, def, 'on');
        if (~isempty(III))
            varname = III{1}; 
            varval = III{2}; if (~isempty(str2num(varval))) varval = str2num(varval); end 
            if (~isfield(pinfo.work, varname)) pinfo.work.(varname) = cell(1, numel(pinfo.(catok).(varok))); end 
            for (i = 1:numel(cellind)) pinfo.work.(varname){cellind(i)} = varval; end
        else
            disp('----------> action canceled');
        end
elseif (sel == 1)
   fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle); 
   for (i = 1:nfield) 
   if (strcmp(fieldtitle{i}, 'work'))
         subfieldnow = fieldnames(pinfo.(fieldtitle{i}));
         fieldselection = getappdata(hfield(i), 'selection'); fieldind = find(fieldselection == 1);
         selvar = subfieldnow(fieldind); nsel = numel(fieldind); str = [];
         if (nsel > 0)
            vartypenow = cell(1, nsel);
            for (tt = 1:nsel) 
                str = [str, selvar{tt}, ';'];
                vartypenow{tt} = resolvetypeformultiple(pinfo.work.(selvar{tt})(cellind));
            end
            vartype = vartypenow{1}; ook = 1;
            for (tt = 2:nsel)
                if (~strcmp(vartypenow{tt}, vartype)) ook = 0; end
            end
            if (ook)
            pp = {'Enter the new variable name'}; 
            prompstr = ['Define from these work variables: ', str]; 
            def = {'var1'};
            III = inputdlg(pp, prompstr, 1, def, 'on');
            if (~isempty(III))
                varname = III{1}; 
                promptstr = 'Definition';
                [flag, sel] = askvectorkeyword(promptstr, vartype, 'crossvar'); 
                if (~isfield(pinfo.work, varname)) pinfo.work.(varname) = cell(1, numel(pinfo.(catok).(varok))); end 
                for (j = 1:numel(cellind))
                    if (nsel == 2)
                        nn1 = numel(pinfo.work.(selvar{1}){cellind(j)});  nn2 = numel(pinfo.work.(selvar{2}){cellind(j)});
                        if (nn1 == 1) | (nn2 == 1)
                            if (nn1 == 1) t1now = 1; t2now = 2; end
                            if (nn2 == 1) t1now = 2; t2now = 1; end
                            if strncmpi(flag, 'divide',3) & strncmpi(vartype, 'vectorcell',10)
                                pinfo.work.(varname){cellind(j)} = pinfo.work.(selvar{t2now}){cellind(j)} ./ pinfo.work.(selvar{t1now}){cellind(j)};
                            elseif strncmpi(flag, 'minus',3) & strncmpi(vartype, 'vectorcell',10)
                                pinfo.work.(varname){cellind(j)} = pinfo.work.(selvar{t2now}){cellind(j)} - pinfo.work.(selvar{t1now}){cellind(j)};
                            elseif strncmpi(flag, 'plus',3) & strncmpi(vartype, 'vectorcell',10)
                                pinfo.work.(varname){cellind(j)} = pinfo.work.(selvar{t2now}){cellind(j)} + pinfo.work.(selvar{t1now}){cellind(j)};
                            elseif strncmpi(flag, 'multiply',3) & strncmpi(vartype, 'vectorcell',10)
                                pinfo.work.(varname){cellind(j)} = pinfo.work.(selvar{t2now}){cellind(j)} .* pinfo.work.(selvar{t1now}){cellind(j)};
                            end
                        else
                            nn = numel(pinfo.work.(selvar{1}){cellind(j)});
                            for (tt = 2:nsel) nn = min([nn numel(pinfo.work.(selvar{tt}){cellind(j)})]); end
                            if (nsel ==2) & strncmpi(flag, 'divide',3) & strncmpi(vartype, 'vectorcell',10)
                                pinfo.work.(varname){cellind(j)} = pinfo.work.(selvar{1}){cellind(j)}(1:nn) ./ pinfo.work.(selvar{2}){cellind(j)}(1:nn);
                            elseif (nsel ==2) & strncmpi(flag, 'minus',3) & strncmpi(vartype, 'vectorcell',10)
                                pinfo.work.(varname){cellind(j)} = pinfo.work.(selvar{1}){cellind(j)}(1:nn) - pinfo.work.(selvar{2}){cellind(j)}(1:nn);
                            elseif (nsel ==2) & strncmpi(flag, 'plus',3) & strncmpi(vartype, 'vectorcell',10)
                                pinfo.work.(varname){cellind(j)} = pinfo.work.(selvar{1}){cellind(j)}(1:nn) + pinfo.work.(selvar{2}){cellind(j)}(1:nn);
                            elseif (nsel ==2) & strncmpi(flag, 'multiply',3) & strncmpi(vartype, 'vectorcell',10)
                                pinfo.work.(varname){cellind(j)} = pinfo.work.(selvar{1}){cellind(j)}(1:nn) .* pinfo.work.(selvar{2}){cellind(j)}(1:nn);
                            end
                        end
                    else %%%if more than 2 variables selected
                        nn = numel(pinfo.work.(selvar{1}){cellind(j)});
                        for (k = 1:nn)
                        if (strncmpi(vartype, 'vectorcell',10))
                           valnow = NaN*ones(nsel,1);
                           for (tt = 1:nsel) valnow(tt) = pinfo.work.(selvar{tt}){cellind(j)}(k); end
                           pinfo.work.(varname){cellind(j)}(k) = squeezenow(valnow, flag, sel, vartype);
                        else
                           valnow = cell(nsel,1);
                           for (tt = 1:nsel) valnow{tt} = pinfo.work.(selvar{tt}){cellind(j)}(k); end 
                           pinfo.work.(varname){cellind(j)}{k} = squeezenow(valnow, flag, sel, vartype);
                        end
                        end
                    end
                end
            else
                disp('----------> action canceled');
            end
            else
                disp('----------> variable types do not match');
            end
         else
            disp('----------> no selected work variables');
         end
   end
   end
end
end

function [pinfo,data, ok] = squeezevals(pinfo,data, cellind, hfield)
ok = 1; 
input = questdlg('Are you sure to squeeze the values of the selected variables?', 'Squeeze values', 'No');
if (strcmp(input, 'Yes'))
    fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle);
    for (i = 1:nfield) %search selected sub-fields in each big field
        if (strcmp(fieldtitle{i}, 'work'))
         subfield{i} = fieldnames(pinfo.(fieldtitle{i}));
         fieldselection = getappdata(hfield(i), 'selection'); fieldselectindex = find( fieldselection == 1);
         for (j = 1:numel(fieldselectindex) )
             vartype = resolvetypeformultiple(pinfo.(fieldtitle{i}).(subfield{i}{fieldselectindex(j)})(cellind));
             if (strcmp(vartype, 'cellcell')) | (strncmpi(vartype, 'vectorcell',10))
                [flag, sel] = askvectorkeyword(strcat(fieldtitle{i},'.',subfield{i}{fieldselectindex(j)}), vartype, 'crossvalue');
                 for (kkk = 1:numel(cellind))
                     if iscell(pinfo.(fieldtitle{i}).(subfield{i}{fieldselectindex(j)}))
                         val = pinfo.(fieldtitle{i}).(subfield{i}{fieldselectindex(j)}){cellind(kkk)};
                         pinfo.(fieldtitle{i}).(subfield{i}{fieldselectindex(j)}){cellind(kkk)} = squeezenow(val, flag, sel, vartype);
                     else
                         val = pinfo.(fieldtitle{i}).(subfield{i}{fieldselectindex(j)})(cellind(kkk));
                         pinfo.(fieldtitle{i}).(subfield{i}{fieldselectindex(j)})(cellind(kkk)) = squeezenow(val, flag, sel, vartype);
                     end
                 end
             end
         end
        end
    end  
else
    ok = 0;
end

function [out, ind] = squeezenow(val, flag, sel, vartype, varargin)
out = []; ind = [];
if (strncmpi(flag, 'all', 3)) | (strncmpi(flag, 'default', 3))
    out = val; ind = 1:numel(val); %%%if flag = default = all: already assigned
elseif (strncmpi(flag, 'manual',3))
    iii = intersect( 1:numel(val), sel );
    out = val(iii); ind = iii;
elseif (strncmpi(flag, 'equal', 3))
    ind = [];
    if isnumeric(sel)
       for (i = 1:numel(sel))
           iii = find(val == sel(i)); ind = union(ind, iii);
       end
    else
       for (ijk = 1:numel(val))
           DT = strfind(lower(val{ijk}), lower(sel));
           if ~iscell(DT)
               if ~isempty(DT) ind = union(ind, ijk); end
           else
               for (ttj = 1:numel(DT))
                   if ~isempty(DT{ttj}) ind = union(ind, ijk); break; end
               end
           end
           %if ~isempty(strfind(lower(val{ijk}), lower(sel))) ind = union(ind, ijk); end
       end
       %iii = find(strncmpi(val, sel, numel(sel))); ind = union(ind, iii);
    end
    out = val(ind);
elseif (strncmpi(flag, 'notequal', 5))
    ind = 1:numel(val);
    if isnumeric(sel)
       for (i = 1:numel(sel))
           iii = find(val == sel(i)); ind = setdiff(ind, iii);
       end
    else
       for (ijk = 1:numel(val))
           DT = strfind(lower(val{ijk}), lower(sel));
           if ~iscell(DT)
               if ~isempty(DT) ind = setdiff(ind, ijk); end
           else
               for (ttj = 1:numel(DT))
                   if ~isempty(DT{ttj}) ind = setdiff(ind, ijk); break; end
               end
           end
           %if ~isempty(strfind(lower(val{ijk}), lower(sel))) ind = setdiff(ind, ijk); end
       end
       %iii = find(strncmpi(val, sel, numel(sel))); ind = setdiff(ind, iii);
    end
    out = val(ind);
elseif (strncmpi(flag, 'anyone',3))
    aa = rand; 
    if (aa==0) 
        out = val(1); ind = 1;
    else
        ind = ceil(rand*numel(val));
        out = val(ind);
    end
elseif (strcmp(vartype, 'vectorcell1'))
    if (nargin < 5)  %%%if mm, ss not provided, compute locally
        mm = mean(val(~isnan(val))); ss = std(val(~isnan(val)));
    else
        mm = varargin{1}; ss = varargin{2};
    end
    if (strncmpi(flag, 'mean',3))
        out = mean(val(~isnan(val)));
    elseif (strncmpi(flag, 'min',3))
        [out,ind] = min(val);
    elseif (strncmpi(flag, 'max',3))
        [out,ind] = max(val);    
    elseif (strncmpi(flag, 'median',3))
        out = median(val(~isnan(val)));
    elseif (strncmpi(flag, 'larger',3))
        ind = find(val>=sel); out = val(ind);
    elseif (strncmpi(flag, 'smaller',3))
        ind = find(val<=sel); out = val(ind);
    elseif (strncmpi(flag, 'std',3))
        out = std(val(~isnan(val)));
    elseif (strcmp(flag, 'norm-mean'))
        out = val - mm;
    elseif (strcmp(flag, 'norm/mean'))
        if (mm ~=0) out = val/mm; end
    elseif (strcmp(flag, 'norm/std'))
        if (ss ~=0) out = val/ss; end
    elseif (strcmp(flag, 'norm-mean/std'))
        if (ss ~=0) out = (val-mm)/ss; end
    elseif (strcmp(flag, '+'))
        out = val + sel(1);
    elseif (strcmp(flag, '-'))
        out = val - sel(1);
    elseif (strcmp(flag, '*'))
        out = val * sel(1);
    elseif (strcmp(flag, '/'))
        out = val / sel(1);
    elseif (strcmp(flag, '-/'))
        out = (val - sel(1))/sel(2);
    end
elseif (strcmp(vartype, 'vectorcell2'))
    tt = NaN*ones(size(val));
    for (i = 1:numel(val))
        if (~isempty(val{i}))
            tt(i) = val{i};
        end
    end
    if (nargin < 5)  %%%if mm, ss not provided, compute locally
        mm = mean(tt(~isnan(tt))); ss = std(tt(~isnan(tt)));
    else
        mm = varargin{1}; ss = varargin{2};
    end
    if (strncmpi(flag, 'mean',3))
        out = mean(tt(~isnan(tt)));
    elseif (strncmpi(flag, 'min',3))
        [out,ind] = min(tt);
    elseif (strncmpi(flag, 'max',3))
        [out,ind] = max(tt);    
    elseif (strncmpi(flag, 'median',3))
        out = median(tt(~isnan(tt)));
    elseif (strncmpi(flag, 'larger',3))
        ind = find(tt>=sel); out = tt(ind);
    elseif (strncmpi(flag, 'smaller',3))
        ind = find(tt<=sel); out = tt(ind);
    elseif (strncmpi(flag, 'std',3))
        out = std(tt(~isnan(tt)));
    elseif (strcmp(flag, 'norm-mean'))
        out = val - mm;
    elseif (strcmp(flag, 'norm/mean'))
        if (mm ~=0) out = val/mm; end
    elseif (strcmp(flag, 'norm/std'))
        if (ss ~=0) out = val/ss; end
    elseif (strcmp(flag, 'norm-mean/std'))
        if (ss ~=0) out = (val-mm)/ss; end
    elseif (strcmp(flag, '+'))
        out = tt + sel(1);
    elseif (strcmp(flag, '-'))
        out = tt - sel(1);
    elseif (strcmp(flag, '*'))
        out = tt * sel(1);
    elseif (strcmp(flag, '/'))
        out = tt / sel(1);
    elseif (strcmp(flag, '-/'))
        out = (tt - sel(1))/sel(2);
    end
end

function [flag, sel] = askvectorkeyword(prompstr, type, crosstype)
flag = []; sel = 0;  %%%crosstype = 'crossvar', 'crossvalue', 'crosscell'
switch crosstype
    case 'crossvalue'
         if (strncmpi(type, 'vectorcell', 10))
 pp = {'Enter squeezing method (min; max; mean; std; median; +-*/; norm-mean; norm/mean; norm/std; norm-mean/std; all; any; manual; larger; smaller; equal; notequal'; 'If .+-*/manual/larger/smaller/equal/notequal. is selected, enter one or more value indices (1,2...):'};
         elseif (strcmp(type, 'cellcell'))
 pp = {'Enter squeezing method (all; any; manual; equal; notequal)'; 'If .manual/equal/notequal. selected, Enter a value:'};
         end
    case 'crossvar'
         if (strncmpi(type, 'vectorcell', 10))
 pp = {'Enter defining method (plus; minus; multiply; divide; min; max; mean; median; any'; 'Do not enter below (default 0):'};
         elseif (strcmp(type, 'cellcell'))
 pp = {'any; equal'; 'If .equal. selected, Enter a value:'};
         end
    case 'normalize'
         if (strncmpi(type, 'vectorcell', 10))
 pp = {'Enter normalization method (+-*/; norm-mean; norm/mean; norm/std; norm-mean/std'; 'If .+-*/. selected, enter a constant (1,2...):'};
         end
end
def = {'any'; '0'};
III=inputdlg(pp, prompstr, 2, def, 'on');
if (~isempty(III))
     flag = III{1};  
     if (strncmpi(flag, 'equal', 3)) && strcmp(type, 'cellcell') 
         sel = III{2}; 
     elseif (strncmpi(flag, 'notequal', 5)) && strcmp(type, 'cellcell') 
         sel = III{2};
     elseif (~isempty(str2num(III{2})))
         sel = str2num(III{2}); 
     end    
end

function ttype = resolvetypeformultiple(feature)
%disp(feature{1});
ttype = 'none';
if (isnumeric(feature))
    ttype = 'vector';
elseif (iscellstr(feature))
    ttype = 'stringcell';
elseif (iscell(feature)) %probe the real value in each individual element and assign array type
    for (i = 1:numel(feature))
        if (~isempty(feature{i}))
            if (iscellstr(feature{i}))
                ttype = 'cellcell';
            elseif (isnumeric(feature{i}))
                ttype = 'vectorcell1';
            elseif (ischar(feature{i}))
                if (isempty(str2num(feature{i}))) %if not a numeric string
                   ttype = 'stringcell';
                else
                   ttype = 'numericstrcell';
                end
            elseif iscell(feature{i})
                for (j = 1:numel(feature{i}))
                    if (~isempty(feature{i}{j}))
                        if (isnumeric(feature{i}{j}))
                            ttype = 'vectorcell2';
                        else
                            ttype = 'cellcell';
                        end
                        break
                    end
                end
            end
            break
        end
    end
end

function [pinfo,data, ok] = reassignvals(pinfo,data, cellind, hfield)
fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle); ok = 1;
for (i = 1:nfield) 
     if (strcmp(fieldtitle{i}, 'work'))
         subfieldnow = fieldnames(pinfo.(fieldtitle{i}));
         fieldselection = getappdata(hfield(i), 'selection'); fieldind = find(fieldselection == 1);
         seleparmname = subfieldnow(fieldind); nsel = numel(fieldind); 
         if (nsel > 0)    
             choice = {'Normalize from existing values'; 'Text to numbers'; 'New manual assignment'};
             [sel,ook] = listdlg('ListString', choice, 'SelectionMode', 'single', 'PromptString', 'Reassigning method');
             if ook & (sel == 1) %%%if (global) normalize across cells
                for (tt = 1:nsel)
                    vartype = resolvetypeformultiple(pinfo.work.(seleparmname{tt})(cellind));
                    if (strncmpi(vartype, 'vectorcell',10))
                        [flag, sel] = askvectorkeyword(strcat('work.',seleparmname{tt}), vartype, 'normalize');
                        if ~iscell(pinfo.work.(seleparmname{tt}))
                           val = pinfo.work.(seleparmname{tt})(cellind); mm = mean(val(~isnan(val))); ss = std(val(~isnan(val)));
                           pinfo.work.(seleparmname{tt})(cellind) = squeezenow(val, flag, sel, vartype, mm, ss);
                        else
                           val = [];
                           for (kkk = 1:numel(cellind))
                               val = [val; pinfo.work.(seleparmname{tt}){cellind(kkk)}];
                           end
                           mm = mean(val(~isnan(val))); ss = std(val(~isnan(val)));
                           val = [];
                           for (kkk = 1:numel(cellind))
                               val = pinfo.work.(seleparmname{tt}){cellind(kkk)};
                               pinfo.work.(seleparmname{tt}){cellind(kkk)} = squeezenow(val, flag, sel, vartype, mm, ss);
                           end
                        end
                    end
                end
             elseif ook & (sel == 2)
                for (tt = 1:nsel)
                    for (kkk = 1:numel(cellind))
                         val = pinfo.work.(seleparmname{tt}){cellind(kkk)}; 
                         if ~isnumeric(val)
                            valout = [];
                            if (iscell(val))
                                for (ppp = 1:numel(val))
                                     valnow = findnumberoutoftext(val{ppp});
                                     valout = [valout valnow];
                                end
%                                for (ppp = 1:numel(val))
%                                    ccc = regexp(val{ppp}, '\d', 'match'); %%%this get all the numbers out
%                                    strnow = [];
%                                    for (qqq = 1:numel(ccc))
%                                        strnow = strcat(strnow, ccc{qqq});
%                                    end
%                                    valout(ppp) = str2num(strnow);
%                               end
                            else
                               valout = findnumberoutoftext(val);
%                                strnow = [];
%                                    for (qqq = 1:numel(ccc))
%                                        strnow = strcat(strnow, ccc{qqq});
%                                    end
%                                    valout = str2num(strnow);
                            end
                            pinfo.work.(seleparmname{tt}){cellind(kkk)} = valout;
                         end
                    end
                end
             elseif ook & (sel == 3) %%%if manual assignment
             for (tt = 1:nsel)
                def{tt} = []; typ{tt} = [];  pp{tt} = ['Enter value for ', seleparmname{tt}, ':']; 
                if iscell(pinfo.work.(seleparmname{tt}))
                    if iscell(pinfo.work.(seleparmname{tt}){cellind(1)})
                       typ{tt} = 'cell'; def{tt} = rearangecell(pinfo.work.(seleparmname{tt}){cellind(1)});
                    elseif ischar(pinfo.work.(seleparmname{tt}){cellind(1)})
                       typ{tt} = 'str'; def{tt} = pinfo.work.(seleparmname{tt}){cellind(1)};
                    else
                       typ{tt} = 'numM'; def{tt} = rearrangeMnum(pinfo.work.(seleparmname{tt}){cellind(1)});
                    end
                else
                    typ{tt} = 'numS'; def{tt} = num2str(pinfo.work.(seleparmname{tt})(cellind(1))); 
                end
            end
            III=inputdlg(pp, 'Values for selected variables', nsel, def, 'on');
            if (~isempty(III))
                if (~isempty(typ{tt}) & (~isempty(III{tt}))) 
                    for (tt = 1:nsel)
                        if (strcmp(typ{tt}, 'numS'))
                            pinfo.work.(seleparmname{tt})(cellind) = ones(numel(cellind),1)*str2num(III{tt});
                        elseif (strcmp(typ{tt}, 'cell'))
                            for (kk = 1:numel(cellind))
                                pinfo.work.(seleparmname{tt}){cellind(kk)} = arrangebackcell(III{tt});
                            end
                        elseif (strcmp(typ{tt}, 'numM'))
                            for (kk = 1:numel(cellind))
                                pinfo.work.(seleparmname{tt}){cellind(kk)} = arrangebackMnum(III{tt})';
                            end 
                        elseif (strcmp(typ{tt}, 'str'))
                            for (kk = 1:numel(cellind))
                                pinfo.work.(seleparmname{tt}){cellind(kk)} = III{tt};
                            end 
                        end
                    end
                    %disp('*****entered the inner loop');
                else
                    disp('---------> unknown variable type or no input'); ok = 0;
                end
            else
                ok = 0; disp('---------> no input or action cancelled'); 
            end
            end
        else
                ok = 0; disp('---------> no work varaibles selected'); 
         end
    end
end 

function [pinfo,data, ok] = copyvars(pinfo,data, hfield)
ok = 1; 
fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle);
for (i = 1:nfield) 
        subfield = fieldnames(pinfo.(fieldtitle{i})); 
        fieldselection = getappdata(hfield(i), 'selection'); 
        for (j = 1:numel(subfield))
            if (j<=numel(fieldselection))&&(fieldselection(j) == 1) %%if a subfield is selected, ask for where to copy
                str = [fieldtitle{i} '.' subfield{j}];
                prompstr = ['Copy variable: ']; pp{1} = ['Enter a category to copy variable ', str, ' to'];
                def = {'work'};
                III=inputdlg(pp, prompstr, 1, def, 'on');
                if (~isempty(III))
                   catnow = III{1}; 
                   if (isfield(pinfo, catnow))
                       if (isempty(pinfo.(catnow)))
                       %pinfo.(catnow) = setfield(pinfo.(catnow), subfield{j}, pinfo.(fieldtitle{i}).(subfield{j}));
                           pinfo.(catnow)(1).(subfield{j}) = pinfo.(fieldtitle{i}).(subfield{j});
                       else
                           nstr = numel(subfield{j}); %%%%this is to rename the copied field if already exist
                           nfnames = fieldnames(pinfo.(catnow));
                           nfound = numel(find(strncmpi(nfnames, subfield{j}, nstr)));
                           if (nfound > 0) %%%if already a fields
                               newname = strcat(subfield{j},num2str(nfound+1));
                               pinfo.(catnow).(newname) = pinfo.(fieldtitle{i}).(subfield{j});
                           else
                              pinfo.(catnow).(subfield{j}) = pinfo.(fieldtitle{i}).(subfield{j});
                           end
                       end
                       
                   else
                       disp(['----------> warning: entered category ', catnow, ' not found in the database']);
                   end
                else
                   ok = 0; break
                end
            end
        end
end 

function [pinfo,data, ok] = renamevars(pinfo,data, hfield)
ok = 1; 
fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle);
for (i = 1:nfield) 
    if (strcmp(fieldtitle{i}, 'work'))
        subfield = fieldnames(pinfo.(fieldtitle{i}));
        fieldselection = getappdata(hfield(i), 'selection');
        for (j = 1:numel(subfield))
            if (fieldselection(j) == 1) %%if a subfield is selected, ask for where to copy
                str = [fieldtitle{i} '.' subfield{j}];
                prompstr = ['Rename variable: ']; pp{1} = ['Enter a new name for varaible: ', str];
                def{1} = subfield{j};
                III=inputdlg(pp, prompstr, 1, def, 'on');
                if (~isempty(III))
                   newnam = III{1}; 
                   if (~isfield(pinfo.(fieldtitle{i}), newnam)) | (~strcmp(newnam, subfield{j}))
                       pinfo.(fieldtitle{i}).(newnam) = pinfo.(fieldtitle{i}).(subfield{j});
                       pinfo.(fieldtitle{i}) = rmfield(pinfo.(fieldtitle{i}), subfield{j});
                   else
                       disp(['----------> warning: field already exist ', newnam]);
                   end
                else
                   ok = 0; break
                end
            end
        end
    end
end 

function [pinfo,data, ok] = removevars(pinfo,data, hfield)
ok = 1; 
input = questdlg('Are you sure to remove the selected variables?', 'Remove variables', 'No');
if (strcmp(input, 'Yes'))
    fieldtitle = fieldnames(pinfo); nfield = numel(fieldtitle);
    for (i = 1:nfield) 
        subfield = fieldnames(pinfo.(fieldtitle{i}));
        fieldselection = getappdata(hfield(i), 'selection');
        for (j = 1:numel(subfield))
            if (fieldselection(j) == 1) %%if a subfield is selected, delete it
                pinfo.(fieldtitle{i}) = rmfield(pinfo.(fieldtitle{i}), subfield{j});
            end
        end
    end 
else
    ok = 0;
end

function [pinfo,data, ok] = removecells(pinfo,data, groupind, catok, varok)
ok = 1;
if (~isempty(groupind))
           %%%%first remove group lists: the rest in restind
           allind = 1:numel(data.grouplist.groupname);
           restind = setdiff(allind, groupind);
           %%%%now keep the remaining cells
           rmcellind = []; oldclname = [];
           for (i = 1:numel(groupind)) %%%get the removed cell indices
               oldindexnow = data.grouplist.groupindex{groupind(i)};
               rmcellind = union(rmcellind, oldindexnow);
           end
           restcellind = setdiff(1:numel(pinfo.(catok).(varok)), rmcellind);
           for (i = 1:numel(restind))
               oldcellindnow = data.grouplist.groupindex{restind(i)};
               oldparmfile{i} = pinfo.(catok).(varok)(oldcellindnow);
           end
           excep = []; pinfo = keepspikepinfo(pinfo, restcellind, excep);
           excep = 'grouplist'; data = keepspikepinfo(data, restcellind, excep);
           %%%%re-assign the groupindex
           ngrpp = 0; groupnow = [];
           for (i = 1:numel(restind))
               inddd = []; ncell = 0;
               for (j = 1:numel(oldparmfile{i}))
                   kkk = find(strcmp(pinfo.(catok).(varok), oldparmfile{i}{j}));
                   if (numel(kkk)>1)
                       disp('-------------> group indices may be out of order');
                   elseif (numel(kkk) == 1)
                       ncell = ncell + 1; inddd(ncell) = kkk;
                   end
               end
               if (~isempty(inddd))
                   ngrpp = ngrpp + 1;
                   groupnow.groupindex{ngrpp} = inddd;
                   groupnow.groupname{ngrpp} = data.grouplist.groupname{restind(i)}; 
                   groupnow.grouptype{ngrpp} = data.grouplist.grouptype{restind(i)}; 
                   groupnow.groupcrit{ngrpp} = data.grouplist.groupcrit{restind(i)};
               end
           end
           if (isempty(groupnow))
               data.grouplist = [];data.grouplist.groupname{1} = 'List0'; data.grouplist.groupindex{1} = 1:numel(pinfo.general.parmfile);
               data.grouplist.grouptype{1} = 'Manual'; data.grouplist.groupcrit{1} = [];
           else
               data.grouplist = groupnow;
           end
else
           disp('-----------> no cell groups selected'); ok = 0;

end
   
function pinfo = keepspikepinfo(pinfo, spikeind, excep)
%keep spikeind in the data structure pinfo
%%%%%get all the fields in pinfo
fieldlist = fieldnames(pinfo);
nfield = numel(fieldlist);
for (iii = 1:nfield)
    if (~strcmp(fieldlist{iii}, excep)) && isstruct(pinfo.(fieldlist{iii}))
           subfieldlist = fieldnames(pinfo.(fieldlist{iii})); %%all the subfields
           for (jjj = 1:numel(subfieldlist))
               pinfo.(fieldlist{iii}).(subfieldlist{jjj}) = pinfo.(fieldlist{iii}).(subfieldlist{jjj})(spikeind);
           end
    end
end

function [F, B, P] = findcolors(str)
F = []; B = []; P = [];
[F, tok] = strtok(str, ';');
if (numel(tok)>1)
    tok = tok(2:numel(tok)); [B, tok] = strtok(tok, ';');
    if (numel(tok) > 1)
        tok = tok(2:numel(tok)); [P, tok] = strtok(tok, ';');
    end
end

function [writefilename, okk] = getoutputfile(hf, ow, dbtype)
okk = 1; writefilename = [];
if (ow == 0)
   [fname, pname] = uiputfile(fullfile(cd, strcat('*', dbtype)), 'Write the new spike database to:');
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

function def = rearangecell(val) % val could be a single or multiple strings
if iscell(val) %multiple strings
    def = [];
    for (i = 1:numel(val))
        def = [def ';' val{i}];
    end
else %single string
    def = val;
end

function val = arrangebackcell(def)
val = []; 
if (~isempty(strfind(def, ';')))
    nv = 1; [str, tok] = strtok(def, ';'); val{nv} = str;
    while (~isempty(tok))
        [str,tok] = strtok(tok, ';'); nv = nv + 1; val{nv} = str;
    end
else
    val = def;
end

function def = rearrangeMnum(val) %re-arrange a multiple number array
nv = numel(val); def = [];
for (i = 1:nv)
    def = [def ';' num2str(val(i))];
end

function val = arrangebackMnum(def)
val = []; 
if (~isempty(strfind(def, ';')))
    nv = 1; [str, tok] = strtok(def, ';'); val(nv) = str2num(str);
    while (~isempty(tok))
        [str,tok] = strtok(tok, ';'); nv = nv + 1; val(nv) = str2num(str);
    end
else
    val = str2num(def);
end

function valout = findnumberoutoftext(val)
valout = [];
if ~isempty(val)
numind = regexp(val, '\d'); %%%this gives all the digit character indices
iii = find(diff(numind)>1); %%%find the number breaks
if isempty(iii)
   valout = str2num(val(numind));
else
   valout = NaN*ones(1, numel(iii)+1);
   valout(1) = str2num(val(1:iii(1)));
   for (ijk = 2:numel(iii))
       valout(ijk) = str2num(val(numind(iii(ijk-1)+1):numind(iii(ijk))));
   end
   valout(numel(iii)+1) = str2num(val((numind(iii(numel(iii))+1)):numind(numel(numind))));
end
end