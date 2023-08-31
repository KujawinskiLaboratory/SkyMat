% Noah Germolus 06 May 2021
% This script is designed to be a wrapper for considerSkyline.m, and both
% are based on the considerMAVEN/riMAVEN code by Krista Longnecker. 
% The objective of these combined files is to take output from Skyline
% (peak areas from UPLC-Orbitrap data) and convert it to concentrations by
% using a standard curve as a ratio (light/heavy).

%
% Yuting Zhu start editing the script 3/23/23

clear

%% Set filenames
fileBase = 'testing_Compilation_3isotopes'; % Set this, don't mess with the automatic date system.
today = datestr(datetime('now'),'.yyyy.mm.dd');
NameOfFile = string([fileBase,today,'_D5.mat']);

%% Set the sequence file here.
wDir = '/Volumes/KujLab/Brianna/Projects/Skyline_MATLAB_codeCompilation/testing/3_isotopes_Skyline/raw_file_example';
fName = 'CMP_Yuting_Exomtab_pos_redo_022323.xlsx';
sampleInfoFile = string([wDir filesep fName]);

clear wDir

%% Set the location and names of the quantification tables exported from Skyline
sDir = '/Volumes/KujLab/Brianna/Projects/Skyline_MATLAB_codeCompilation/testing/3_isotopes_Skyline/Skyline';
dfile_pos = string([sDir filesep 'Skyline_pos_test_3isotopes_QuantTable_manipulated.csv']);
dfile_neg = string([sDir filesep 'CINAR_quant_neg.csv']); 
clear sDir

%% Set directory for where SkyMat codes are - this will create an output folder for your results
oDir = '/Volumes/whoi/dept/mcg/KujLab/Brianna/Projects/Skyline_MATLAB_codeCompilation/testing/3_isotopes_Skyline/SkyMat';
oFolder = string([oDir filesep 'Output']);
mkdir(oFolder);

cd(oFolder);

clear oDir 

%% ConsiderSkyline processing for positive mode.

units = 'ng'; %set unit for standard curve (e.g., ng or pg)

[pos_D5.sNames, pos_D5.kgd] = considerSkyline(dfile_pos, sampleInfoFile,...
    'pos','heavyD5',2, units, oFolder);

%% ConsiderSkyline processing for negative mode.

[neg_D5.sNames, neg_D5.kgd] = considerSkyline(dfile_neg, sampleInfoFile,...
 'neg','heavyD5',2, units, oFolder);

%% MERGING DATA FROM TWO MODES
clear fName fileBase today

mtabNames_D5 = sort(cat(1,[neg_D5.kgd.names + " neg"],[pos_D5.kgd.names + " pos"]));
if length(unique(mtabNames_D5)) ~= length(mtabNames_D5)
    error('Something is wrong - duplicate names in the list of metabolites')
end

%for the pooled samples (and perhapds others), I will have duplicate sets 
%of names with either _pos or _neg appended; 
tInfo_D5 = readtable(sampleInfoFile);
clear sampleInfoFile

%before I dive into the unknowns, remove anything that has goodData = 0
k = find(tInfo_D5.goodData==0);
tInfo_D5(k,:) = [];
clear k

%%first, go through and iterate through the pooled samples
%%to provide numbers for these (otherwise will have duplicate
%%names)
%%NOTE: update names of pooled samples here
%s = strcmp(tInfo.Sample_Name,'pool');

%YZ 03.31.2023 use contains instead of strcmp so can find the pooled in my dataset
s = contains(tInfo_D5.Sample_Name,'pool') & contains(tInfo_D5.Sample_Name,'pos');
ks = find(s==1);
for a = 1:length(ks)
    t = tInfo_D5.Sample_Name(ks(a));
    tInfo_D5.Sample_Name(ks(a)) = strcat('pool',num2str(a,'%02.f'),t); %YZ 03.31.2023 added '%02.f'
    clear t
end
clear a ks a

s = contains(tInfo_D5.Sample_Name,'pool') & contains(tInfo_D5.Sample_Name,'neg');
ks = find(s==1);
for a = 1:length(ks)
    t = tInfo_D5.Sample_Name(ks(a));
    tInfo_D5.Sample_Name(ks(a)) = strcat('pool',num2str(a,'%02.f'),t); %YZ 03.31.2023 added '%02.f'
    clear t
end
clear a ks a

%now find the Unknown...should have the same number for positive and
%negative ion mode bc have pruned out the different QC samples already
s = strcmp(tInfo_D5.Sample_Type,'Unknown');
sp = strcmp(tInfo_D5.ionMode,'pos');
ksp = (find(s==1 & sp==1));
sn = strcmp(tInfo_D5.ionMode,'neg');
ksn = (find(s==1 & sn==1));

if ~isequal(length(ksp),length(ksn))
    error('Something wrong, these should be the same length')
end
clear s sp sn ksp ksn

%%parse out the names. Use this to figure out the unique samples and setup
%%a new matrix that I can propagate with the metabolites from both positive
%%and negative ion mode. Bit of a hack, and growing worse.
nrow = size(tInfo_D5,1);
tInfo_D5.type = repmat({''},nrow,1);
tInfo_D5.cName = repmat({''},nrow,1);
%examples of additional columns used in the BIOS-SCOPE project
% tInfo.cruise = repmat({''},nrow,1);
% tInfo.cast = zeros(nrow,1);
% tInfo.niskin = zeros(nrow,1);
% tInfo.depth = zeros(nrow,1);
% tInfo.addedInfo = repmat({'none'},nrow,1);

for a = 1:nrow
    if strcmp(tInfo_D5.Sample_Type{a},'Unknown') %only do unknowns      
        one = tInfo_D5.Sample_Name{a};
        r_pooled = regexp(one,'pool');
        r_spiked = regexp(one,'NOSPIKE');

            if r_spiked
                %pooled with 500 ng/ml spike
                if 0
                    %keep it
                    tInfo_D5.type(a) = {'spiked'};
                    tInfo_D5.cName(a) = {'spiked'};
                else
                    %skip
                end
            elseif r_pooled
                %pooled sample
                tInfo_D5.type(a) = {'pooled'};
                %put the number of this pooled sample into 'addedInfo'
                r_nL = regexp(one,'p'); %lower case
                r_nU = regexp(one,'C'); %upper case
                %tInfo.addedInfo(a) = {one(r_nL+1 : r_nU-1)};
                tInfo_D5.addedInfo(a) = {'pooled'};
                tInfo_D5.cName(a) = {one(1:r_nU-1)};
            else
                %actual sample
                tInfo_D5.addedInfo(a) = {'sample'}; %redundant...
                tInfo_D5.cName(a) = {one(1:end-4)};
                %fprintf('here')
            end
        clear one r_* under
    end
end
clear a nrow

sInfo_D5 = table;
sInfo_D5.cName = unique(tInfo_D5.cName);
%the first row of this will be empty, delete that
if isequal(sInfo_D5.cName(1),{''})
    sInfo_D5(1,:) = [];
end


%now make an empty matrix for the data...will be all numbers so no need for
%special format
mtabData_D5 = zeros(size(mtabNames_D5,1),size(sInfo_D5,1));
%need to track some additional details:
mtabDetails_D5 = table();

%get the index for rows for positive AND negative mtabs:

kgdNames = [pos_D5.kgd.names + " pos";neg_D5.kgd.names + " neg"]; 
[c idx_New idx_Old] = intersect(mtabNames_D5,kgdNames);
all_LOD = [pos_D5.kgd.LOD;neg_D5.kgd.LOD]; 
LOD_D5 = all_LOD(idx_Old);
all_LOQ = [pos_D5.kgd.LOQ;neg_D5.kgd.LOQ]; 
LOQ_D5 = all_LOQ(idx_Old);


clear c idx_New idx_Old all_LOD kgdNames all_LOQ

[c idx_posNew idx_posOld] = intersect(mtabNames_D5,pos_D5.kgd.names + " pos");
[c idx_negNew idx_negOld] = intersect(mtabNames_D5,neg_D5.kgd.names + " neg");


mtabDetails_D5.mode(idx_posNew,1) = {'pos'};
mtabDetails_D5.mode(idx_negNew,1) = {'neg'};

sInfo_D5.runOrder_pos(:,1) = 0;
sInfo_D5.runOrder_neg(:,1) = 0;

sInfo_D5.FileName_pos(:,1) = {''};
sInfo_D5.FileName_neg(:,1) = {''};

for a = 1:size(sInfo_D5,1)
    s = strcmp(sInfo_D5.cName(a),tInfo_D5.cName);
    ks = find(s==1);
    if length(ks) ~= 2
        error('Something is wrong, should be two of each')
    end
    
    %some variant of this:
    for aa = 1:2
        %propagate sInfo with the cast/depth/etc. information, only do once
%         if aa == 1
%             sInfo.type(a) = tInfo.type(ks(aa));
%             sInfo.cName(a) = tInfo.cName(ks(aa));
%             sInfo.cruise(a) = tInfo.cruise(ks(aa));
%             sInfo.cast(a) = tInfo.cast(ks(aa));
%             sInfo.niskin(a) = tInfo.niskin(ks(aa));
%             sInfo.depth(a) = tInfo.depth(ks(aa));
%             sInfo.addedInfo(a) = tInfo.addedInfo(ks(aa));
%         end
        
        im = tInfo_D5.ionMode{ks(aa)};
        if isequal(im,'pos')
            tName = tInfo_D5.File_Name(ks(aa));
            sInfo_D5.FileName_pos(a,1) = tName;

            [c ia tIdx] =intersect(tName,pos_D5.sNames);
            mtabData_D5(idx_posNew,a) = pos_D5.kgd.goodData(idx_posOld,tIdx);
            clear c ia tIdx tName
            
        elseif isequal(im,'neg')
            tName = tInfo_D5.File_Name(ks(aa));
            sInfo_D5.FileName_neg(a,1) = tName;

            [c ia tIdx] =intersect(tName,neg_D5.sNames);
            mtabData_D5(idx_negNew,a) = neg_D5.kgd.goodData(idx_negOld,tIdx);
            clear c ia tIdx tName
        else 
            error('Something wrong')
        end
        clear im
    end
    clear aa s ks        
end
clear a

clear idx_*

clear r s

for a = 1: size(sInfo_D5,1)
    %do positive ion mode first
    gc = sInfo_D5{a,'FileName_pos'}{:}; %added {:} to deal with table output
    t = regexp(gc,'_');
    if ~isempty(t)
        sInfo_D5.runOrder_pos(a,1) = str2num(gc(t(end)+1:end));
    else
        sInfo_D5.runOrder_pos(a,1) = NaN;
    end
    clear gc t
    
    %then negative ion mode
    gc = sInfo_D5{a,'FileName_neg'}{:}; %added {:} to deal with table output
    t = regexp(gc,'_');
    if ~isempty(t)
        sInfo_D5.runOrder_neg(a,1) = str2num(gc(t(end)+1:end));
    else
        sInfo_D5.runOrder_neg(a,1) = NaN;
    end
    clear gc t
end
clear a
 
clear a dfile_neg dfile_pos neg_info pos_info sampleInfoFile_neg ...
    sampleInfoFile_pos
 
save(NameOfFile)


