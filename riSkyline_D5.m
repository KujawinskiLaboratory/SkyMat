% Noah Germolus 06 May 2021
% This script is designed to be a wrapper for considerSkyline.m, and both
% are based on the considerMAVEN/riMAVEN code by Krista Longnecker. 
% The objective of these combined files is to take output from Skyline
% (peak areas from UPLC-Orbitrap data) and convert it to concentrations by
% using a standard curve as a ratio (light/heavy).

%%
clear

%%
codeDir = 'C:\Users\brianna.garcia\Documents\GitHub\SkyMat'; % User must set code base directory

cd(codeDir) %Change directory to work from base code directory
addpath(genpath(codeDir)) %Add all folders and subfolders within the code diretory to your path

%% Set filenames
fileBase = 'SkyMat_testing_3isotopes'; % Set this, don't mess with the automatic date system.
today = datestr(datetime('now'),'.yyyy.mm.dd');
NameOfFile = string([fileBase,today,'_D5.mat']);

%% Set the sequence file here.
wDir = '\Example_Dataset\Example_Input'; %User must set directory where the  sequence file is location
fName = 'SkyMat_3isotopes_test_pos_and_neg.xlsx'; %User must set the sequence file name here
sampleInfoFile = string([wDir filesep fName]);

clear wDir
%% Set the location and names of the quantification tables exported from Skyline
sDir = 'Example_Dataset\Example_Input'; %User must set directory where the quantification table file(s) are location
dfile_pos = string([sDir filesep 'SkyMat_3isotopes_test_pos.csv']); %User must set the file name for the positive mode quantification table
dfile_neg = string([sDir filesep 'SkyMat_3isotopes_test_neg.csv']); %User must set the file name for the negative mode quantification table
clear sDir

%% Set directory for where SkyMat codes are - this will create an output folder for your results
oDir = string([fileBase,today,'_Output']); %Creates an output folder based on your previously defined fileBase and todays date
oFolder = strcat(oDir,filesep,'D5'); %Create a D5 output folder within the base output directory created above
mkdir(oFolder); 
addpath(genpath(oDir)) %Add this new output folder and all subfolders to the path

cd(oFolder) %Change directory to your output folder so all downstream processing outputs collect in this folder

clear oDir 
%% ConsiderSkyline processing for positive mode.

units = 'ng'; %set unit for standard curve 
% acceptable units are ng, pg, ng/mL, and pg/mL, note that the units are case sensitive

[pos_D5.sNames, pos_D5.kgd] = considerSkyline(dfile_pos, sampleInfoFile,...
    'pos','heavyD5',2, units, oFolder);

%% ConsiderSkyline processing for negative mode.

[neg_D5.sNames, neg_D5.kgd] = considerSkyline(dfile_neg, sampleInfoFile,...
 'neg','heavyD5',2, units, oFolder);

%% Save temporary file before merging data 

save('temp_D5');

%% MERGING DATA FROM TWO MODES
clear fName today

% The traditional approach here is to take both metabolite lists, positive
% and negative mode, and keep both sets of data and append the ion mode to
% the metabolite name. In the future, there may be a different routine that
% automatically calibrates each metabolite at each sample using both the
% mode and isotope that give the tightest error bounds, eliminating this
% step. 
mtabNames_D5 = sort(cat(1,[neg_D5.kgd.names + " neg"],[pos_D5.kgd.names + " pos"]));
if length(unique(mtabNames_D5)) ~= length(mtabNames_D5)
    error('Something is wrong - duplicate names in the list of metabolites')
end

% For the pooled samples (and perhapds others), I will have duplicate sets 
% of names with either _pos or _neg appended; 
tInfo_D5 = readtable(sampleInfoFile);
clear sampleInfoFile

% Before I dive into the unknowns, remove anything that has goodData = 0
% This step does take place within considerSkyline, but we're re-reading
% the sample info file here to create permanent variables for the
% workspace, which need to be re-pruned.
k = find(tInfo_D5.goodData==0);
tInfo_D5(k,:) = [];
clear k

% Parse out the names. Use this to figure out the unique samples and setup
% a new matrix that I can propagate with the metabolites from both positive
% and negative ion mode. Bit of a hack, and growing worse.
% NPG 20 Sept 2023: I think this whole section might need to be removed.
% We're adding extra columns for parsing out sample metadata, which is
% something I do in downstream processing or have straight-up in the sample
% info table. It doesn't really do much good to have this "hack" present in
% what's supposed to be the basic processing script. 
nrow = size(tInfo_D5,1);
tInfo_D5.cName = repmat({''},nrow,1);

% First, go through and iterate through the pooled samples
% to provide numbers for these (otherwise will have duplicate
% names). Need to do separately for both modes.
s = contains(tInfo_D5.SampleName,'pool') & contains(tInfo_D5.SampleName,'pos');
ks = find(s==1);
for a = 1:length(ks)
    t = tInfo_D5.SampleName(ks(a));
    tInfo_D5.SampleName(ks(a)) = strcat('pool',num2str(a,'%02.f'),'_',t); %YZ 03.31.2023 added '%02.f'
    tInfo_D5.cName(ks(a)) = {strcat('pool',num2str(a,'%02.f'))};
    clear t
end
clear a ks 

s = contains(tInfo_D5.SampleName,'pool') & contains(tInfo_D5.SampleName,'neg');
ks = find(s==1);
for a = 1:length(ks)
    t = tInfo_D5.SampleName(ks(a));
    tInfo_D5.SampleName(ks(a)) = strcat('pool',num2str(a,'%02.f'),'_',t); %YZ 03.31.2023 added '%02.f'
    tInfo_D5.cName(ks(a)) = {strcat('pool',num2str(a,'%02.f'))};
    clear t
end
clear a ks 

% Now find the Unknown...should have the same number for positive and
% negative ion mode.
s = strcmp(tInfo_D5.SampleType,'Unknown');
sp = strcmp(tInfo_D5.ionMode,'pos');
ksp = (find(s==1 & sp==1));
sn = strcmp(tInfo_D5.ionMode,'neg');
ksn = (find(s==1 & sn==1));

if ~isequal(length(ksp),length(ksn))
    error('Something wrong, these should be the same length')
end
clear s sp sn ksp ksn


% examples of additional columns used in the BIOS-SCOPE project
% tInfo_D5.cruise = repmat({''},nrow,1);
% tInfo_D5.cast = zeros(nrow,1);
% tInfo_D5.niskin = zeros(nrow,1);
% tInfo_D5.depth = zeros(nrow,1);
% tInfo_D5.addedInfo = repmat({'none'},nrow,1);

for a = 1:nrow
    if strcmp(tInfo_D5.SampleType{a},'Unknown') %only do unknowns      
        one = tInfo_D5.SampleName{a};
        r_pooled = regexp(one,'pool');
            if r_pooled
                %put the type of this pooled sample into 'addedInfo'
                tInfo_D5.addedInfo(a) = {'pooled'};
            else
                %actual sample
                tInfo_D5.addedInfo(a) = {'sample'}; %redundant...'
                if contains(one, " pos")
                tInfo_D5.cName(a) = {erase(one, " pos")};
                elseif contains(one," neg")
                tInfo_D5.cName(a) = {erase(one, " neg")};
                %fprintf('here')
                end 
            end
        clear one r_* under
    end
end
clear a nrow

% NPG 20 Sept 2023: This used to take five lines. Not sure why. But, this
% makes a table with the sample names as the first column. 
sInfo_D5 = table(unique(tInfo_D5.cName), 'VariableNames',{'cName'});

if isequal(sInfo_D5.cName(1),{''})
    sInfo_D5(1,:) = [];
end

% Preallocate double-type matrix for metabolite data.
mtabData_D5 = zeros(size(mtabNames_D5,1),size(sInfo_D5,1));
mtabData_D5_filtered = zeros(size(mtabNames_D5,1),size(sInfo_D5,1));

% Need to track some additional details; namely which file came from which
% ion mode.
mtabDetails_D5 = table();

% Get the index for rows for positive AND negative mtabs and reorder. 
kgdNames = [pos_D5.kgd.names + " pos";neg_D5.kgd.names + " neg"]; 
[c idx_New idx_Old] = intersect(mtabNames_D5,kgdNames);
all_LOD = [pos_D5.kgd.LOD;neg_D5.kgd.LOD]; 
LOD_D5 = all_LOD(idx_Old);
all_LOQ = [pos_D5.kgd.LOQ;neg_D5.kgd.LOQ]; 
LOQ_D5 = all_LOQ(idx_Old);
all_r2 = [pos_D5.kgd.r2_line;neg_D5.kgd.r2_line];
r2_line_D5 = all_r2(idx_Old);


clear c idx_New idx_Old all_LOD kgdNames all_LOQ

[c idx_posNew idx_posOld] = intersect(mtabNames_D5,pos_D5.kgd.names + " pos");
[c idx_negNew idx_negOld] = intersect(mtabNames_D5,neg_D5.kgd.names + " neg");


mtabDetails_D5.mode(idx_posNew,1) = {'pos'};
mtabDetails_D5.mode(idx_negNew,1) = {'neg'};

sInfo_D5.runOrder_pos(:,1) = 0;
sInfo_D5.runOrder_neg(:,1) = 0;

sInfo_D5.FileName_pos(:,1) = {''};
sInfo_D5.FileName_neg(:,1) = {''};

% This section takes the ordered sample names and metabolite names and
% reshuffles the mode-specific calibrated measurements into a single
% matrix. It also contains some of the metadata hack from earlier that
% should be removed (commented lines).
for a = 1:size(sInfo_D5,1)
    s = strcmp(sInfo_D5.cName(a),tInfo_D5.cName);
    ks = find(s==1);
    % The section starts by searching a sample name, anticipating both a
    % pos and neg mode for each sample. 
    if length(ks) ~= 2
        error('Something is wrong, should be two of each')
        % If you get this error, check to see if your goodData column is
        % properly pruning your data so that there's the same number of
        % files for each mode, AND that all your cNames are actually the
        % same across modes--typos happen. 
    end
    
    for aa = 1:2
        %propagate sInfo_D5 with the cast/depth/etc. information, only do once
%         if aa == 1
%             sInfo_D5.type(a) = tInfo_D5.type(ks(aa));
%             sInfo_D5.cName(a) = tInfo_D5.cName(ks(aa));
%             sInfo_D5.cruise(a) = tInfo_D5.cruise(ks(aa));
%             sInfo_D5.cast(a) = tInfo_D5.cast(ks(aa));
%             sInfo_D5.niskin(a) = tInfo_D5.niskin(ks(aa));
%             sInfo_D5.depth(a) = tInfo_D5.depth(ks(aa));
%             sInfo_D5.addedInfo(a) = tInfo_D5.addedInfo(ks(aa));
%         end
        % Two cases, because depending on the ionMode, we're shifting data
        % from a different struct into the data matrices.
        im = tInfo_D5.ionMode{ks(aa)};
        if isequal(im,'pos')
            tName = tInfo_D5.FileName(ks(aa));
            RunOrder = tInfo_D5.runOrder(ks(aa));
            sInfo_D5.FileName_pos(a,1) = tName;
            sInfo_D5.runOrder_pos(a,1) = str2num(string(RunOrder));

            [c ia tIdx] =intersect(tName,pos_D5.sNames);
            mtabData_D5(idx_posNew,a) = pos_D5.kgd.goodData(idx_posOld,tIdx);
            mtabData_D5_filtered(idx_posNew,a) = pos_D5.kgd.goodData_filtered(idx_posOld,tIdx);
            clear c ia tIdx tName
            
        elseif isequal(im,'neg')
            tName = tInfo_D5.FileName(ks(aa));
            sInfo_D5.FileName_neg(a,1) = tName;
            RunOrder = tInfo_D5.runOrder(ks(aa));
            sInfo_D5.runOrder_neg(a,1) = str2num(string(RunOrder));

            [c ia tIdx] =intersect(tName,neg_D5.sNames);
            mtabData_D5(idx_negNew,a) = neg_D5.kgd.goodData(idx_negOld,tIdx);
            mtabData_D5_filtered(idx_negNew,a) = neg_D5.kgd.goodData_filtered(idx_negOld,tIdx);
            clear c ia tIdx tName
        else 
            error('Something wrong')
        end
        clear im RunOrder
    end
    clear aa s ks        
end
clear a

clear idx_*

clear r s
 
clear a dfile_neg dfile_pos neg_info pos_info sampleInfoFile_neg ...
    sampleInfoFile_pos
 
save(NameOfFile)

%% Use the convertMoles.m function to convert from mass to concentration 
%(e.g., pg to pM)
% input variables for function include:
% tDir - directory where your transition list is found that includes
% columns for isParent and StdMW
%tFile - the name of the Transition list file in .csv format.
%mtabNames - this can be either _D5, _D5, or the _filtered version of
%those
%units - acceptable units are ng, pg, ng/mL, and pg/mL, note that the units are case sensitive
%volume in mL - for example here '25' as a numeric input

tDir = 'Example_Input'; %User must set in the directory location where the transition list is located
tFile = string([tDir filesep 'TransitionList_SkyMat_Example.xlsx']);

mtabData_D5_conc = convertMoles(tFile, mtabNames_D5, mtabData_D5, units, 25);
mtabData_D5_conc_filtered = convertMoles(tFile, mtabNames_D5, mtabData_D5_filtered, units, 25);

LOD_D5_conc = convertMoles(tFile, mtabNames_D5, LOD_D5, units, 25);
LOQ_D5_conc = convertMoles(tFile, mtabNames_D5, LOQ_D5, units, 25);

save(NameOfFile)

%% Create wide-format compiled table with metabolite names,LOD, and LOQ

%update units variable to make compatible for saving to csvfile
units = strrep(units,'/','_per_');

%determine concentration units based on unit input
if strcmp(units,"ng") || strcmp(units,"ng_per_mL") 
     conc_units = "nM"; 
elseif strcmp(units,"pg") || strcmp(units,"pg_per_mL") 
        conc_units = "pM";
end 

% Save the unfiltered dataset converted to concentration
conc_table_D5 = splitvars(table(mtabNames_D5,r2_line_D5,LOD_D5_conc,LOQ_D5_conc,mtabData_D5_conc));
conc_table_D5.Properties.VariableNames = [{'Metabolite','r2_line','LOD','LOQ'},sInfo_D5.cName'] ;

writetable(conc_table_D5, strcat(fileBase,"_",conc_units,'_concTable_D5.csv'));

%Save the dataset with values < LOD filtered out and converted to concentration
conc_table_D5_filtered = splitvars(table(mtabNames_D5,r2_line_D5,LOD_D5_conc,LOQ_D5_conc,mtabData_D5_conc_filtered));
conc_table_D5_filtered.Properties.VariableNames = [{'Metabolite','r2_line','LOD','LOQ'},sInfo_D5.cName'] ;

writetable(conc_table_D5_filtered, strcat(fileBase,"_",conc_units,'_concTable_D5_filtered.csv'));

%Save updated MATLAB file
save(NameOfFile)