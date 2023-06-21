function mtabData_nM = convertMoles(negTransitions, posTransitions, mtabNames, mtabData, volume_mL)

%% Convert from ng added to nM concentrations
% This process involves reloading the transition list files and matching
% names to molecular weights (MW). For DT5, I had to manually correct 
% many, many names. That's why the data file loaded at the beginning 
% contains two variables for mtabNames.
% 
% Edited to include inital sample volume - this enables to code to convert
% from ng added (the standard unit used in the lab for standard curves) to
% nM by first converting to ng/mL (based on your input volume) and then
% using the MWs from the Transition Lists to convert each compound to nM. 
% 20230620 BMG 

posInfo = readtable(posTransitions);
posInfo(posInfo.isParent == 0,:) = [];
MWp = [posInfo(:,1), posInfo(:,14)]; %please confirm column 14 is the MW or the code will not execute properly
negInfo = readtable(negTransitions);
negInfo(negInfo.isParent == 0,:) = [];
MWn = [negInfo(:,1), negInfo(:,14)]; %%please confirm column 14 is the MW or the code will not execute properly
MW = [MWp;MWn];
clear MWp MWn
MW = unique(MW, 'rows');
clear posInfo negInfo

% Making both compound name columns into strings and removing the neg/pos
% identifier.
MW.CompoundName = string(MW.CompoundName);
mtabNamesAgnostic = strrep(strrep(mtabNames, ' pos', ''),' neg','');


% Time to index where each unique molecule is found in mtabNames.
[~, iNames] = ismember(mtabNamesAgnostic, MW.CompoundName);
iNames(iNames==0) = [];
% Use those indices to sort MW values.
if sum(mtabNamesAgnostic == MW.CompoundName(iNames)) ~= length(mtabNames)
    disp("name mismatch")
    return
end
MWtoConvert = MW.StdMW(iNames);
% convert from ng added to ng/L to nM
mtabData_nM = ((mtabData./volume_mL).*1000)./MWtoConvert;

end