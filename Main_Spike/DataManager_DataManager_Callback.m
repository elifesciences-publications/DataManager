function hmain = DataManager_DataManager_Callback
%%This is a big standalone function for analyzing spike and EEG data

[MCroot, MCname, DAname, DEname, ANname, MAname] = CurrentVersion;

hf = figure('Name', MAname, 'NumberTitle', 'off', 'NextPlot', 'add',...
    'MenuBar', 'figure', 'Units', 'normalized', 'Position', [0.05 0.2 0.9 0.7]);
uimenu('Parent', hf, 'Label', '          ');

hfile = uimenu('Parent', hf, 'Label', '  Session  ');
    uimenu(hfile, 'Label', 'Generate New SpikeDatabase', 'Callback', 'DataManager_GenerateSpikeDatabase_Callback');
    uimenu(hfile, 'Label', 'Generate New BehavDatabase', 'Callback', 'DataManager_GenerateBehavDatabase_Callback');

    uimenu(hfile, 'Label', 'Combine Database', 'Callback', 'DataManager_CombineDatabase_Callback', 'Separator', 'on');

    uimenu(hfile, 'Label', 'Open Database', 'Callback', 'DataManager_OpenDatabase_Callback', 'Separator', 'on');
    uimenu(hfile, 'Label', 'SaveAs', 'Callback', 'DataManager_SaveAs_Callback', 'Separator', 'on');
    uimenu(hfile, 'Label', 'SaveSelect', 'Callback', 'DataManager_SaveSelect_Callback');
    uimenu(hfile, 'Label', 'SaveStripped', 'Callback', 'DataManager_SaveStripped_Callback');
    uimenu(hfile, 'Label', 'Save', 'Callback', 'DataManager_Save_Callback');
    uimenu(hfile, 'Label', 'Close', 'Callback', 'DataManager_Close_Callback');
    
hrecompute = uimenu('Parent', hf, 'Label', '  Recompute  ');    
    uimenu(hrecompute, 'Label', 'ReComputeDataBase_GoodLapsOnly', 'Callback', 'DataManager_ReComputeDatabase_GoodLaps');
    uimenu(hrecompute, 'Label', 'ReComputeDataBase_NoStops', 'Callback', 'DataManager_ReComputeDatabase_NoStops');
    
hcross = uimenu('Parent', hf, 'Label', '  Cross Analysis  ');  
    hcrr = uimenu(hcross, 'Label', 'Cross/auto-correlation', 'Callback', 'DataManager_GenerateCrrDatabase_Callback');

    

hppro = uimenu('Parent', hf, 'Label', '  Projects  ');
      
    hRTT = uimenu('Parent', hppro, 'Label', 'Rett Mice');
       uimenu(hRTT, 'Label', 'Identify Pairs with Overlapping Place Fields', 'Callback', 'DataManager_RettMice_OverlappingPlaceCellPairs');
       uimenu(hRTT, 'Label', 'Quantify Event Participation', 'Callback', 'DataManager_RettMice_EventParticipation');
       uimenu(hRTT, 'Label', 'Revise Place Field Properties', 'Callback', 'DataManager_RettMice_RevisePlaceFieldProperties');   
if nargout == 1
   hmain = hf;
end