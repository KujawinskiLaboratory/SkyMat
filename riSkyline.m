%% Noah Germolus 06 May 2021
% This script is designed to be a wrapper for considerSkyline.m, and both
% are based on the considerMAVEN/riMAVEN code by Krista Longnecker. 
% The objective of these combined files is to take output from Skyline
% (peak areas from UPLC-Orbitrap data) and convert it to concentrations by
% using a standard curve as a ratio (light/heavy).

clear

% Set filenames
fileBase = 'DT2_NPG'; % Set this, don't mess with the automatic date system.
today = datestr(datetime('now'),'.yyyy.mm.dd');
NameOfFile = string(['DT2_NPG',today,'.mat']);

% There is not a peak quality parameter I've found yet in Skyline.
% setQuality = 

% Set the sequence file here.
wDir = 'F:/Metabolomics/Derivatization2/2021_0405_FB2/sequence_fromMethods';
fName = 'mtab_Noah_FB2_040521_asv.xlsx';
sampleInfoFile = string([wDir filesep fName]);
clear wDir

% Where are the lists of exported exported data from Skyline?
sDir = 'F:/Metabolomics/Derivatization2/SkylineExports';
dfile_pos = string([sDir filesep 'Quant_pos_17May2021.csv']);
dfile_neg = string([sDir filesep 'Quant_neg_17May2021.csv']); 
clear sDir

% Unclear whether we need an SRM list here but that's where it is in
% riMAVEN

% This is also a place for divergence: Skyline doeesn't export the same
% janky csvs that MAVEN does, so we will use high-level import functions.

% Move onto the processing for positive mode.
[pos.sNames, pos.kgd] = considerSkyline(dfile_pos, sampleInfoFile, 'pos');
[neg.sNames, neg.kgd] = considerSkyline(dfile_neg, sampleInfoFile, 'neg');

clear fName fileBase today

% MERGING DATA FROM TWO MODES

mtabNames = sort(cat(1,[neg.kgd.names + " neg"],[pos.kgd.names + " pos"]));
if length(unique(mtabNames)) ~= length(mtabNames)
    error('Something is wrong - duplicate names in the list of metabolites')
end


%for the pooled samples (and perhapds others), I will have duplicate sets 
%of names with either _pos or _neg appended; 
tInfo = readtable(sampleInfoFile);
clear sampleInfoFile

%%first, go through and iterate through the pooled samples
%%to provide numbers for these (otherwise will have duplicate
%%names)
%%NOTE: update names of pooled samples here
s = strcmp(tInfo.Sample_Name,'NOPOOLS');
ks = find(s==1);
for a = 1:length(ks)
    t = tInfo.Sample_Name(ks(a));
    tInfo.Sample_Name(ks(a)) = strcat('p',num2str(a),t);
    clear t
end
clear a ks a

%before I dive into the unknowns, remove anything that has goodData = 0
k = find(tInfo.goodData==0);
tInfo(k,:) = [];
clear k

%now find the Unknown...should have the same number for positive and
%negative ion mode bc have pruned out the different QC samples already
s = strcmp(tInfo.sType,'sample');
sp = strcmp(tInfo.ionMode,'pos');
ksp = (find(s==1 & sp==1));
sn = strcmp(tInfo.ionMode,'neg');
ksn = (find(s==1 & sn==1));

if ~isequal(length(ksp),length(ksn))
    error('Something wrong, these should be the same length')
end
clear s sp sn ksp ksn

%%parse out the names. Use this to figure out the unique samples and setup
%%a new matrix that I can propagate with the metabolites from both positive
%%and negative ion mode. Bit of a hack, and growing worse.
nrow = size(tInfo,1);
tInfo.type = repmat({''},nrow,1);
tInfo.cName = repmat({''},nrow,1);
%examples of additional columns used in the BIOS-SCOPE project
% tInfo.cruise = repmat({''},nrow,1);
% tInfo.cast = zeros(nrow,1);
% tInfo.niskin = zeros(nrow,1);
% tInfo.depth = zeros(nrow,1);
% tInfo.addedInfo = repmat({'none'},nrow,1);

for a = 1:nrow
    if strcmp(tInfo.sType{a},'sample') %only do unknowns      
        one = tInfo.Sample_Name{a};
        r_pooled = regexp(one,'NOPOOL');
        r_spiked = regexp(one,'NOSPIKE');

            if r_spiked
                %pooled with 500 ng/ml spike
                if 0
                    %keep it
                    tInfo.type(a) = {'spiked'};
                    tInfo.cName(a) = {'spiked'};
                else
                    %skip
                end
            elseif r_pooled
                %pooled sample
                tInfo.type(a) = {'pooled'};
                %put the number of this pooled sample into 'addedInfo'
                r_nL = regexp(one,'p'); %lower case
                r_nU = regexp(one,'P'); %upper case
                %tInfo.addedInfo(a) = {one(r_nL+1 : r_nU-1)};
                tInfo.addedInfo(a) = {'pooled'};
                tInfo.cName(a) = {one(1:r_nU-1)};
            else
                %actual sample
                tInfo.addedInfo(a) = {'sample'}; %redundant...
                tInfo.cName(a) = {one(1:end-4)};
                %fprintf('here')
            end
        clear one r_* under
    end
end
clear a nrow

sInfo = table;
sInfo.cName = unique(tInfo.cName);
%the first row of this will be empty, delete that
if isequal(sInfo.cName(1),{''});
    sInfo(1,:) = [];
end


%now make an empty matrix for the data...will be all numbers so no need for
%special format
mtabData = zeros(size(mtabNames,1),size(sInfo,1));
%need to track some additional details:
mtabDetails = table();

%get the index for rows for positive AND negative mtabs:
[c idx_posNew idx_posOld] = intersect(mtabNames,pos.kgd.names + " pos");
[c idx_negNew idx_negOld] = intersect(mtabNames,neg.kgd.names + " neg");

mtabDetails.mode(idx_posNew,1) = {'pos'};
mtabDetails.mode(idx_negNew,1) = {'neg'};

sInfo.runOrder_pos(:,1) = 0;
sInfo.runOrder_neg(:,1) = 0;

sInfo.FileName_pos(:,1) = {''};
sInfo.FileName_neg(:,1) = {''};

for a = 1:size(sInfo,1)
    s = strcmp(sInfo.cName(a),tInfo.cName);
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
        
        im = tInfo.ionMode{ks(aa)};
        if isequal(im,'pos')
            tName = tInfo.File_Name(ks(aa));
            sInfo.FileName_pos(a,1) = tName;

            [c ia tIdx] =intersect(tName,pos.sNames);
            mtabData(idx_posNew,a) = pos.kgd.goodData(idx_posOld,tIdx);
            clear c ia tIdx tName
            
        elseif isequal(im,'neg')
            tName = tInfo.File_Name(ks(aa));
            sInfo.FileName_neg(a,1) = tName;

            [c ia tIdx] =intersect(tName,neg.sNames);
            mtabData(idx_negNew,a) = neg.kgd.goodData(idx_negOld,tIdx);
            clear c ia tIdx tName
        else 
            error('Something wrong')
        end
        clear im
    end
    clear aa s ks        
end
clear a

clear idx_* tInfo

for a = 1: size(sInfo,1)
    %do positive ion mode first
    gc = sInfo{a,'FileName_pos'}{:}; %added {:} to deal with table output
    t = regexp(gc,'_');
    if ~isempty(t)
        sInfo.runOrder_pos(a,1) = str2num(gc(t(end)+1:end));
    else
        sInfo.runOrder_pos(a,1) = NaN;
    end
    clear gc t
    
    %then negativeion mode first
    gc = sInfo{a,'FileName_neg'}{:}; %added {:} to deal with table output
    t = regexp(gc,'_');
    if ~isempty(t)
        sInfo.runOrder_neg(a,1) = str2num(gc(t(end)+1:end));
    else
        sInfo.runOrder_neg(a,1) = NaN;
    end
    clear gc t
end
clear a
 
%making decisions about pos/neg mtabs; this is done either based on past
%experience with compounds, decisions made within a project, by plotting
%the positive and negative ion mode data together, or looking at the plots
%in MAVEN. The first time this code is run, I allow everything to go
%through and then I later rerun the code to make decisions to leave one
%version for each compound
if 0  
    toDelete = {'2deoxyinosine neg';'2deoxyguanosine neg';...
        '3deoxyguanosine neg';'NAD neg';'Sucrose_341 neg';...
        'Trehalose_341 neg';'adenine neg';...
        'adenosine 5''-monophosphate neg';'biotin neg';...
        'desthiobiotin neg';'folic acid neg';'folinic acid neg';...
        'guanosine pos';'n-acetyl muramic acid neg';...
        'pantothenic acid neg';'s-(5''-adenosyl)-L-homocysteine neg';...
        'xanthine neg';'xanthosine neg';'5deoxyadenosine neg';...
        'cytidine neg';'hemin II';'2deoxycytidine neg';...
        'Trehalose pos';'aspartic acid neg';...
        'D-fructose 1_6-bisphosphate'};

    [c ia ib] = intersect(toDelete,mtabNames);

    mtabNames(ib,:)=[];
    mtabDetails(ib,:)=[];
    mtabData(ib,:)=[];
    clear c ia ib toDelete
end


save(NameOfFile)


