function [pinfo, data] = DataManager_ComputeCrrInit_Callback(pinfo, data, cellind, vv)
%compute initial variables of a database
if (~isempty(cellind))
   
   disp('----------> compute cross-/auto-correlations');
   [pinfo, data] = DataManager_FindCrr(pinfo, data, cellind, vv);
  
%get behavioral/fields information -- not done yest!
%%%%%%%%%disp('-------> get behavioral properties for all spikes');
%pinfo = DataManager_FindBehaviorInfo(pinfo);
        
%%get spatial information -- not done yet!
%pinfo = Manage_FindSpatialInfo(pinfo);
  
%get interneuron/pyramidal cell classification -- NOT DONE yet!, don't do it at the beginning: it requires user interference 
% disp('-------> classify neuron type');
% pinfo = Manage_FindNeuronType(pinfo);

%get phase information - NOT done yet!!! 
%disp('-------> get phase properties for all spikes');
%pinfo = Manage_FindPhaseInfo(pinfo); There Maybe a new program doing this - chech TrajSwitch - phase programs
end 


