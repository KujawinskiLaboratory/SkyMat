%% Noah Germolus 17 May 2021
% Now that I have a pipeline to get here from Skyline, it's time to take a
% look at individual metabolites.

clear
clc
setDefaultFigs

%% Part 0: Combining files
% I want to have both sets of quantitifed metabolites (filter and
% dissolved) in one dataset.
load('DT2_NPG_filters.2021.05.18.mat')
sInfoBig = sInfo;
clear NameOfFile
fulldata = struct;
fulldata.filters.neg = neg; fulldata.filters.pos = pos;
clear pos neg
mtabDataBig = mtabData; 
clear mtabData sInfo

load('DT2_NPG_dissolved.2021.05.19.mat')

[~, ia, ib] = intersect(sInfo, sInfoBig);
[~, inota] = setdiff(sInfo, sInfoBig);
[~, inotb] = setdiff(sInfoBig, sInfo);
sInfoBig = [sInfo(inota,:); sInfoBig(inotb,:); sInfoBig(ib,:)];
sInfoBig.matrix(:) = "filter";
sInfoBig.matrix([7:12,19:27]) = "BATS"; 
mtabData = [mtabData(:,inota), mtabData(:,ia)];
mtabDataBig = [mtabDataBig(:,inotb), mtabDataBig(:,ib)]


%% Part 1: visual
% Here I just want to load up the file and show things for the filters;
% expected vs. measured concentrations. 

filtersamples = ismember(sInfo.quantMatrix,'filter');
sInfoSmall = sInfo(filtersamples,:);
G = findgroups(sInfoSmall.group);
for ii=1:length(mtabNames)
    data = mtabData(ii,filtersamples)';
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
    ylim([0,max([avgs + stds + 1;exp])])
    xticklabels({'','buffered','','unbuffered'})
    legend({'Rep Mean +/- \sigma', 'Target'}, 'Location', 'southeast')
end