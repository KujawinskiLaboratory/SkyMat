%% Noah Germolus 17 May 2021
% Now that I have a pipeline to get here from Skyline, it's time to take a
% look at individual metabolites.

clear
clc
load('DT2_NPG.2021.05.17.mat')
setDefaultFigs

%% Part 1: visual
% Here I just want to load up the file and show things for the filters;
% expected vs. measured concentrations. 

filtersamples = ismember(sInfo.matrix,'filter');
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