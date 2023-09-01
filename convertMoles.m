function mtabData_conc = convertMoles(Transitions, mtabNames, mtabData, units, volume_mL)

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

tInfo = readtable(Transitions);

posInfo = tInfo;
posInfo(posInfo.isParent == 0,:) = [];
pLog = strcmp(posInfo.ionMode,'positive');
posInfo = posInfo(pLog,:);
MWp = table([string(posInfo.("MoleculeListName")),posInfo.StdMW]);
MWp = splitvars(MWp,'Var1','NewVariableNames',{'CompoundName','StdMW'});
MWp.StdMW = str2double(MWp.StdMW);

clear pLog

negInfo = tInfo;
negInfo(negInfo.isParent == 0,:) = [];
nLog = strcmp(negInfo.ionMode,'negative');
negInfo = negInfo(nLog,:);
MWn = table([string(negInfo.("MoleculeListName")),negInfo.StdMW]);
MWn = splitvars(MWn,'Var1','NewVariableNames',{'CompoundName','StdMW'});
MWn.StdMW = str2double(MWn.StdMW);

clear nLog

MW = [MWp;MWn];

clear MWp MWn 

MW = unique(MW, 'rows');
clear posInfo negInfo

% Making both compound name columns into strings and removing the neg/pos
% identifier.
MW.CompoundName = string(MW.CompoundName);
MW.CompoundName = stripName(MW.CompoundName);
mtabNamesAgnostic = strrep(strrep(mtabNames, ' pos', ''),' neg','');
mtabNamesAgnostic = stripName(mtabNamesAgnostic);

% Time to index where each unique molecule is found in mtabNames.
[test, iNames] = ismember(mtabNamesAgnostic, MW.CompoundName);
iNames(iNames==0) = [];
% Use those indices to sort MW values.
if sum(strcmp(mtabNamesAgnostic,MW.CompoundName(iNames))) ~= length(mtabNames)
    disp("name mismatch")
    return
end

MWtoConvert = MW.StdMW(iNames);

% convert from mass to concentration (e.g., ng to nM)

mtabData_conc = ((mtabData./volume_mL).*1000)./MWtoConvert;

end