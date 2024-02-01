%% Noah Germolus 16 Aug 2023
% This script is what to run after calibrating Zoop2 metabolites with
% SkyMat. Its intentions are threefold:
% 1. Combine two sets of calibrated concentrations by determining which
% curve (13C or D5 label normalized) to use on a per-sample basis. 
% 2. Calculate LOD/LOQ to use based on the confluence of the two standard
% curves. 
% 3. Convert numbers to molar basis. 

%% Part 1: Loading files. 
% Two files are available and their variables are the same except for an
% appended set of characters, "_D5" or "_C13" depending on which isotope
% was used in the calibration. 

% Added the line below because the riSkyline scripts pull you into whatever
% your working data directory is. We want to actually remain there, but
% may need to remind MATLAB where this script is without resetting your
% working directory. 
addpath('.')

clear
clc
load('Example_Dataset/Example_Output\C13\SkyMat_testing_3isotopes.2024.02.01_C13.mat')
load('Example_Dataset/Example_Output\D5\SkyMat_testing_3isotopes.2024.02.01_D5.mat')

outdir = "Example_Dataset/Example_Output\Combined";
if ~exist("outdir", "dir")
    mkdir(outdir)
    disp("Output directory for combined files created:")
    disp(outdir)
end
filename = "SkyMat_testing_3isotopes_OneMode.mat";

%% Part 2: Curve Metrics.
% So, I was going to arbitrarily select a range where D5 calibrates high
% concentrations and 13C calibrates low concentrations (so, for example, we
% just use the 13C ratio curve data below 700 pg/mL) but this isn't good
% enough for me. Plus, what if the values you get are different in a way
% the code can't handle--say, a sample has 500 pg/mL tryptophan using the
% D5 curve and 800 using the 13C? So, let's calculate the 95% prediction
% interval on each curve and find where one overtakes the other. 

% Generalized formula for 95% PI on a regression yhat +/-
% t(0.025,n-2)*SE
% Where yhat is the predicted value, t(...) is the 2.5%ile of
% the t-distribution with n-2 degrees of freedom, and SE is the standard
% error of the estimate.
% We won't use yhat, because the interval, not the prediction +/- the
% interval, is what we want. 
% SE is tricky, so I've wired considerSkyline to export variables that
% will help us parametrize the whole thing as a continuous function.

% The prediction interval for a certain level of confidence is necessarily
% a smooth function, and we will need two functions, one above and one 
% below the line of calibration.
% One could argue that you might just take [highest slope and highest
% intercept], for example, as a condition for the upper bound of the
% prediction interval, but this doesn't work--the regression is constrained
% tightly at its midpoint and looser at the ends. 

% This got complicated when I realized the way we use our calibration
% curves doesn't lend itself to this approach: concentration is the x-axis,
% and so we actually have to invert the interval through the use of some
% fancy triangle math to get the interval we're after. Furthermore, because
% the normal prediction interval is about the y-direction, it doesn't exist
% for the x-direction depending on which y-based PI you're looking at. Both
% upper and lower need to be calculated, then del-x calculated based on
% congruency and symmetry. (It's symmetric about x and y)

% PI function for y-direction.
predFunc = @(x,A,B,C,xhat,slope,intercept) [B.*sqrt(A).*sqrt(C + (x-xhat).^2),...
     B.*sqrt(A).*sqrt(C + (x-xhat).^2)];
% Inverts the above function to find the delx interval.
XPf = @(x,A,B,C,xhat,slope,intercept) predFunc(x,A,B,C,xhat,slope,intercept)./slope;
% Traces the function above to the x-value where the estimate is valid
xfind = @(x,xpf_out) [x+xpf_out(:,1), x-xpf_out(:,2)];


% I can only perform this analysis for metabolites with a valid curve in
% both isotopes.
mtabNames_both = intersect(mtabNames_D5, mtabNames_C13);
mtabNames_all = union(mtabNames_D5,mtabNames_C13);
mtabData_all = zeros(size(mtabNames_all,1),size(mtabData_C13_conc,2));

LOQ = zeros(size(mtabNames_all,1));
LOD = zeros(size(mtabNames_all,1));

PI95 = nan.*zeros(length(mtabNames_all), size(mtabData_all,2));

w = waitbar(0,'','Name','Checking metabolites for prediction interval choice...');
for ii = 1:length(mtabNames_all)
    waitbar(ii/length(mtabNames_all),w, mtabNames_all(ii));
    % Two initial cases, where the data is only good for one isotope. Just
    % place the full dataset into the final matrix.
    mi13C = find(mtabNames_C13 == mtabNames_all(ii));
    miD5 = find(mtabNames_D5 == mtabNames_all(ii));
    if ~ismember(mtabNames_all(ii),mtabNames_both) && ismember(mtabNames_all(ii),mtabNames_C13)
        mtabData_all(ii,:) = mtabData_C13_conc(mi13C,:);
        LOD(ii) = LOD_C13_conc(mi13C);
        LOQ(ii) = LOQ_C13_conc(mi13C);
    elseif ~ismember(mtabNames_all(ii),mtabNames_both) && ismember(mtabNames_all(ii),mtabNames_D5)
        mtabData_all(ii,:) = mtabData_D5_conc(miD5,:);
        LOD(ii) = LOD_D5_conc(miD5);
        LOQ(ii) = LOQ_D5_conc(miD5);
    else
        % The case where the two must be compared

        xmax = max([mtabData_C13_conc(mi13C,:),mtabData_D5_conc(miD5,:)]);
        if xmax>1e5
            xmax=1e5;
        end
        xt = [0:xmax/5000:max(1,round(2*xmax))]';
        if contains(mtabNames_all(ii), " pos")
            mii13C = find(string(pos_C13.kgd.names) == strrep(mtabNames_all(ii)," pos",""));
            miiD5 = find(string(pos_D5.kgd.names) == strrep(mtabNames_all(ii)," pos",""));
            PI_13C = XPf(xt,pos_C13.kgd.A(mii13C),...
                pos_C13.kgd.B(mii13C),pos_C13.kgd.C(mii13C),...
                pos_C13.kgd.xM(mii13C),pos_C13.kgd.slope(mii13C),...
                pos_C13.kgd.intercept(mii13C));
            x_13C = xfind(xt,PI_13C);
            PI_D5 = XPf(xt,pos_D5.kgd.A(miiD5),...
                pos_D5.kgd.B(miiD5),pos_D5.kgd.C(miiD5),...
                pos_D5.kgd.xM(miiD5),pos_D5.kgd.slope(miiD5),...
                pos_D5.kgd.intercept(miiD5));
            x_D5 = xfind(xt,PI_D5);
        elseif contains(mtabNames_all(ii), " neg")
            mii13C = find(string(neg_C13.kgd.names) == strrep(mtabNames_all(ii)," neg",""));
            miiD5 = find(string(neg_D5.kgd.names) == strrep(mtabNames_all(ii)," neg",""));
            PI_13C = XPf(xt,neg_C13.kgd.A(mii13C),...
                neg_C13.kgd.B(mii13C),neg_C13.kgd.C(mii13C),...
                neg_C13.kgd.xM(mii13C),neg_C13.kgd.slope(mii13C),...
                neg_C13.kgd.intercept(mii13C));
            x_13C = xfind(xt,PI_13C);
            PI_D5 = XPf(xt,neg_D5.kgd.A(miiD5),...
                neg_D5.kgd.B(miiD5),neg_D5.kgd.C(miiD5),...
                neg_D5.kgd.xM(miiD5),neg_D5.kgd.slope(miiD5),...
                neg_D5.kgd.intercept(miiD5));
            x_D5 = xfind(xt,PI_D5);
        end
        % I will use the lower function--that is, the interval below the
        % calibration curve, to calculate PIs for the concentration values,
        % as this covers the lower part of the curve and I will default to
        % whatever lowers the error earliest for low concentrations.
        [x_D, iD] = sort([x_D5(:,1);x_D5(:,2)]);
        [x_C, iC] = sort([x_13C(:,1);x_13C(:,2)]);
        [x_D,iDu,~] = unique(x_D);
        [x_C,iCu,~] = unique(x_C);
        PI_D = [PI_D5(:,1);PI_D5(:,2)];
        PI_D = PI_D(iD);
        PI_D = PI_D(iDu);
        PI_C = [PI_13C(:,1);PI_13C(:,2)];
        PI_C = PI_C(iC);
        PI_C = PI_C(iCu);

        % PI_D = PI_D5(:,2); PI_C = PI_13C(:,2);
        % x_D = x_D5(:,2); x_C = x_13C(:,2);

        % Primacy of 13C: I am going to use comparisons to 13C here, where
        % positive conditions refer to 13C being "better" than D5. I write
        % this mostly for me as I code.

        % Because the back-casted x values don't necessarily line up to
        % each other, this creates a problem. I have three options.
        % 1. Find a way to make the backcasting functions obsolete. I
        % tried, but the simple equivalency of functions method doesn't
        % work so well when I'm doing this interval inversion.
        % 2. Create another linear grid using xt from earlier and use
        % spline interpolation to make the points equivalent. This is sort
        % of cheating, but since these are smooth, sort-of-quadratic
        % functions, I'm not worried about introducing error through this
        % method.
        % 3. Fit a function to the each dataset and find their
        % intersection. This is easier down the line (to find zeros in
        % MATLAB); however, I tried fitting a couple quadratics and it just
        % doesn't give me the accuracy I need.

        % table of interpolated values.
        error1 = (["Not enough valid points to interpolate all conf. curves for "+mtabNames_all(ii)]);
        try vqC = interp1(x_C(x_C>0),PI_C(x_C>0),xt,"nearest");
        catch 
            disp(error1)
        end
        vqD = interp1(x_D(x_D>0),PI_D(x_D>0),xt,"nearest");

        CTI = table(xt,vqC,vqD);

        % First, check if one curve is always better than the other one.
        check1 = (vqC<vqD); % "is the 13C curve better? for which points?"
        if sum(check1)==size(CTI,1) % If all points are better with 13C
            mtabData_all(ii,:) = mtabData_C13_conc(mi13C,:);
            LOD(ii) = LOD_C13_conc(mi13C);
            LOQ(ii) = LOQ_C13_conc(mi13C);
        elseif sum(check1)==0 % If all points are better with D5
            mtabData_all(ii,:) = mtabData_D5_conc(miD5,:);
            LOD(ii) = LOD_D5_conc(miD5);
            LOQ(ii) = LOQ_D5_conc(miD5);
        elseif check1(1) == 1 % if the above are false but 13C starts out better
            kb = find(check1==0);
            xcrit = xt(kb(1));
            kc = mtabData_C13_conc(mi13C,:)<xcrit;
            kd = mtabData_D5_conc(miD5,:)<xcrit;
            if kc~=kd
                %not sure what to do here yet. There's a possibility that
                %a sample, calibrated with both curves, could read
                %different concentrations on different sides of the
                %critical point between the two curves. For now, I will
                %accept the 13C-based critical point as law.
                disp("calibrated values straddle a critical point")
            end
            mtabData_all(ii,kc) = mtabData_C13_conc(mi13C,kc);
            mtabData_all(ii,~kc) = mtabData_D5_conc(miD5,~kc);
            LOD(ii) = LOD_C13_conc(mi13C);
            LOQ(ii) = LOQ_C13_conc(mi13C);

        else % If D5 starts out better and gets worse
            kb = find(check1==1);
            xcrit = xt(kb(1));
            kc = mtabData_C13_conc(mi13C,:)<xcrit;
            kd = mtabData_D5_conc(miD5,:)<xcrit;
            if kc~=kd
                %not sure what to do here yet. There's a possibility that
                %a sample, calibrated with both curves, could read
                %different concentrations on different sides of the
                %critical point between the two curves. For now, I will
                %accept the 13C-based critical point as law.
                disp("calibrated values straddle a critical point")
            end
            mtabData_all(ii,kd) = mtabData_C13_conc(mi13C,kd);
            mtabData_all(ii,~kd) = mtabData_D5_conc(miD5,~kd);
            LOD(ii) = LOD_D5_conc(miD5);
            LOQ(ii) = LOQ_D5_conc(miD5);

        end

        minInt = min([vqC,vqD],[],2);
        mRound = round(mtabData_all(ii,:),2)';
        [Lia, Locb] = ismember(mRound, round(xt,2));
        PI95(ii,Lia) = minInt(Locb(Locb>0))';

        if 1
            % Nobody is obligated to use this, but if you want to view the
            % results of which intervals are better for what
            % concentrations, this will plot all of them. 
            waitbar((ii+0.5)/length(mtabNames_all),w, "plotting");
            set(groot,'defaultFigureVisible','off')
            figure
            plot(xt, real(vqC), "LineWidth",2,"Color","r");
            hold on
            plot(xt, real(vqD), "LineWidth",2,"Color","b")
            title(mtabNames_all(ii))
            xlabel(["Metabolite Concentration "+units])
            ylabel(["Calibration Prediction Interval "+units])
            xline(LOD(ii),'--k','LOD','DisplayName','LOD',"HandleVisibility","off")
            xline(LOQ(ii),':k','LOQ','DisplayName','LOQ',"HandleVisibility","off")
            xline(mtabData_all(ii,:), "-g")
            legend({"^{13}C Curve", "D_5 Curve", "samples"})
            exportgraphics(gca, [outdir+ filesep+ "PredictionIntervals.pdf"], 'Append',  true)
            hold off
            close(gcf)
            set(groot,'defaultFigureVisible','on')
        end
    end

end
close all
delete(w)
clear w

%% Sorting out the ion modes.

uniqueNames = unique(strrep(strrep(mtabNames_all, " neg", ""), " pos", ""));
mtabNames_OneMode = uniqueNames;
mtabData_OneMode = zeros(length(uniqueNames),size(PI95,2));
LOQ_OneMode = zeros(length(uniqueNames),1);
LOD_OneMode = zeros(length(uniqueNames),1);
var_OneMode = zeros(length(uniqueNames),size(PI95,2));
modeUsed = string(zeros(length(uniqueNames),size(PI95,2)));

for ii= 1:length(uniqueNames)
    PosName = uniqueNames(ii) + " pos";
    NegName = uniqueNames(ii) + " neg";
    iallp = find(mtabNames_all == PosName);
    ialln = find(mtabNames_all == NegName);
    if isempty(iallp)
        mtabData_OneMode(ii,:) = mtabData_all(ialln,:);
        LOQ_OneMode(ii,1) = LOQ(ialln);
        LOD_OneMode(ii,1) = LOD(ialln);
        var_OneMode(ii,:) = PI95(ialln,:);
        modeUsed(ii,:) = "-";
        continue
    elseif isempty(ialln)
        mtabData_OneMode(ii,:) = mtabData_all(iallp,:);
        LOQ_OneMode(ii,1) = LOQ(iallp);
        LOD_OneMode(ii,1) = LOD(iallp);
        var_OneMode(ii,:) = PI95(iallp,:);
        modeUsed(ii,:) = "+";
        continue
    end
    varp = PI95(iallp,:); varn = PI95(ialln, :);
    PosBetter = (varp<varn);
    mtabData_OneMode(ii,PosBetter) = mtabData_all(iallp,PosBetter);
    mtabData_OneMode(ii,~PosBetter) = mtabData_all(ialln, ~PosBetter);
    LOQ_OneMode(ii,1) = min([LOQ(iallp),LOQ(ialln)]);
    LOD_OneMode(ii,1) = min([LOD(iallp),LOD(ialln)]);
    var_OneMode(ii,PosBetter) = varp(PosBetter);
    var_OneMode(ii, ~PosBetter) = varn(~PosBetter);
    modeUsed(ii,PosBetter) = "+";
    modeUsed(ii, ~PosBetter) = "-";
end




%% Part 2.1: Cleanup Step 1
% Here I want to get rid of all the stuff I generated in the previous
% section, and create a streamlined file that only contains my LOD, LOQ,
% full mtab data, names, and one set of sample info. 

mtabData = mtabData_OneMode;
mtabNames = mtabNames_OneMode;
mtabNames = strrep(mtabNames, "′","'");
sInfo = sInfo_C13; % These two are the same between isotopes.
tInfo = tInfo_C13;
LOD = LOD_OneMode;
LOQ = LOQ_OneMode;
var = var_OneMode;



save([outdir + filesep + filename],"var", "mtabData", "tInfo", "sInfo","mtabNames", "LOQ", "LOD", "conc_units" )

clear

