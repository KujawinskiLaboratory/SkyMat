%% Noah Germolus 17 May 2021
% Now that I have a pipeline to get here from Skyline, it's time to take a
% look at individual metabolites.

clear
clc
%setDefaultFigs

%% Part 0: Combining files
% I want to have both sets of quantitifed metabolites (filter and
% dissolved) in one dataset.
load('DT2_NPG_filters.2021.05.18.mat')
sInfoFilters = sInfo;
clear NameOfFile
fulldata = struct;
fulldata.filters.neg = neg; fulldata.filters.pos = pos;
clear pos neg
mtabDataFilters = mtabData; 
mtabNamesFilters = mtabNames;
clear mtabData sInfo

load('DT2_NPG_dissolved.2021.05.19.mat')
fulldata.dissolved.pos = pos;
fulldata.dissolved.neg = neg;
clear pos neg

[~, ia, ib] = intersect(sInfo, sInfoFilters);
[~, inota] = setdiff(sInfo, sInfoFilters);
[~, inotb] = setdiff(sInfoFilters, sInfo);
sInfoBig = [sInfo(inota,:); sInfoFilters(inotb,:); sInfoFilters(ib,:)];
sInfoBig.matrix(:) = "filter";
sInfoBig.matrix([1:6,19:27]) = "BATS"; %Matrix with which the injection was quantified
mtabDataBig = [mtabDataFilters(:,inotb), mtabData(:,inota), mtabDataFilters(:,2:7), mtabData(:,[8:10,12:17])];

% All the cal curve stuff is NaN by my own design soooo
sInfoBig(1:12,:)=[]; mtabDataBig(:,1:12)=[];
clear sInfo ia ib mtabData sInfoFilters mtabDataFilters inota inotb 
NameOfFile = "DT2_Working_NPG_2021.05.20.mat";
save(NameOfFile)

%% Part 1: visual
% Here I just want to load up the file and show things for the filters;
% expected vs. measured concentrations. 

filtersamples = ismember(sInfoBig.matrix,'filter');
sInfoSmall = sInfoBig(filtersamples,:);
sInfoSmall.group = [1;1;1;2;2;2];
sInfoSmall.known = 72.5*ones(6,1);
G = findgroups(sInfoSmall.group);
for ii=1:length(mtabNames)
    data = mtabDataBig(ii,filtersamples)';
    avgs = splitapply(@mean, data, G);
    stds = splitapply(@std, data, G);
    exp = splitapply(@mean, sInfoSmall.known, G);
    figure
    errorbar(1:length(avgs),avgs,stds, 'LineStyle', 'none',...
        'Marker', 'x', 'MarkerSize', 3)
    hold on
    scatter(1:length(avgs), exp)
    title(mtabNames(ii))
    ylabel('Concentration, ng/mL')
    xlim([0.5,2.5])
    ylim([0,max([avgs + stds + 5;exp + 5])])
    xticklabels({'','buffered','','unbuffered'})
    legend({'Rep Mean +/- \sigma', 'Target'}, 'Location', 'southeast')
end

dissolvedSamples = ismember(sInfoBig.matrix, "BATS");
sInfoSmall = sInfoBig(dissolvedSamples, :);
sInfoSmall.group = [1;1;1;2;2;2;3;3;3];
sInfoSmall.known = [0;0;0;10;10;10;1;1;1];
G = findgroups(sInfoSmall.group);
for ii=1:length(mtabNames)
    data = mtabDataBig(ii,dissolvedSamples)';
    avgs = splitapply(@mean, data, G);
    stds = splitapply(@std, data, G);
    exp = splitapply(@mean, sInfoSmall.known, G);
    figure
    errorbar(1:length(avgs),avgs,stds, 'LineStyle', 'none',...
        'Marker', 'x', 'MarkerSize', 3)
    hold on
    scatter(1:length(avgs), exp)
    title(mtabNames(ii))
    ylabel('Concentration, ng/mL')
    xlim([0.5,3.5])
    ylim([0,max([avgs + stds + 1;exp])])
    xticklabels({'','VSW 0ng/ml','','VSW 10ng/ml','','BATS 1ng/ml'})
    legend({'Rep Mean +/- \sigma', 'Target'}, 'Location', 'best')
end