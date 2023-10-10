% considerSkyline V0 : This function is called by riSkyline.m and is
% responsible for quantification of metabolites via a five-point standard
% curve. It calculates LOD and LOQ as well. Quantification is based on a
% ratio to heavy spike (13C6) with the exception of any metabolite with a
% '0' in its name. These were not derivatized by benzoyl chloride and
% therefore would not have a heavy label. 

function [sampleNames, keepGoodData] = considerSkyline(...
    exportedSkyline, sampleInfoFile, ionMode, SILISType,nSILIS, units, oFolder)

% INPUT %
% 1. exportedSkyline: name of the *.csv exported from Skyline. Must contain
%   at least columns for peak area (both light and heavy), sample type, sample
%   filename, molecule name, and "analyte concentration" (this is for known
%   concentrations like standards and spikes).
% 2. sampleInfoFile: the modified sequence file containing file names and
%   condition information (*.xslx)
% 3. ionMode: can be "neg" or "pos"
% 4. SILISType: can be 'heavyD5' or 'heavyC13'
% 5. nSILIS: can be 1 or 2. This function gets run twice per isotope (once
%   per ionMode), so even though you specify a single isotope type in the 
%   function call, specify '2' for this parameter if you've got both 
%   labels in your dataset.
% 6. units: unit of standard curve (e.g., 'ng' or 'pg') 
% 7. oFolder: Output folder location. This variable is defined in the code and will
%   create the folder and change directory to this folder. 

% OUTPUT %
% 1. sampleNames: A string vector containing the sample names for each
%   calibrated measurement.
% 2. keepGoodData: A struct containing information such as calibrated
%   concentrations, calibration metrics, and flags.

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
% Skyline has a special way of exporting empty cells, which we change to
% NaN here.
data = readtable(exportedSkyline, 'TreatAsEmpty', "#N/A");
info = readtable(sampleInfoFile); 
%reduce info to only samples in the current analysis.
info = info(strcmp(ionMode, info.ionMode),:); 

% Now we have read the two files necessary for quantification. Since we're
% quantifying based on precursor ions, this is a way to pick that out.
data.precursor = (data.FragmentIonType == "precursor" &...
    data.IsotopeLabelType == "light");

compoundList = table(unique(data.MoleculeName, 'stable'),'VariableNames', {'names'});

% CHECKING CONFIRM IONS
% Skyline does not output data rows where it doesn't find the ion, unlike
% our prior approach to MAVEN, where ions that aren't found simply don't
% appear in the export. The transition won't have any exported data if it
% is entirely absent, but will export rows for samples in which the
% transition is found. 
% *This means that for a given ion, #samples does not necessarily equal
% #records. 
% Here I will check for the confirm ion by evaluating if the confirm had a 
% nonzero peak area in *all* non-blank samples. I have not made it a
% requirement for calibration that a confirm is present, as Skyline isn't
% as good at picking up the MS2 fragments at low concentrations. Best to
% check in MZMine and not rely on this step. 

% This section, rather than printing to the cmd window, will for each
% compound write to the diary textfile if a compound is missing MS2 peaks. 
diaryFilename = append(ionMode,"_",SILISType,"_considerSkyline_flags.txt");
if exist(append(oFolder, filesep, diaryFilename), 'file') == 2
    fprintf('The file "%s" already exists in the working directory.\n', diaryFilename);
    fprintf('Please delete the file before proceeding.\n');
    error('Diary file already exists. Code execution halted.');
else
    % Create the diary file or perform other operations
    diary(diaryFilename);
    fprintf('Diary file "%s" created.\n', diaryFilename)

    for a = 1:length(compoundList.names)
    quant = (strcmp(compoundList.names(a), data.MoleculeName) &...
        ~strcmp("blank", data.SampleType) &...
        data.precursor==1);
    confirm = (strcmp(compoundList.names(a), data.MoleculeName) &...
        ~strcmp("blank", data.SampleType) &...
        data.precursor==0);
        if sum(confirm)<sum(quant)
            missing = 100*(sum(quant)-sum(confirm))/sum(quant);
            disp(string([compoundList.names(a)+" is missing confirm ions in "+...
                string(missing)+" percent of injections."]))
        end
    end

diary off

end 

% Get rid of any files you don't want analyzed. We add a column in the
% sequence file called "goodData" (boolean) and would typically only import
% the "good" files into Skyline, but this is a place to make sure we only
% have what we need. 
pruneData = 1; % Generally should always be 1.
if pruneData
    badNames = info.FileName(info.goodData ==0) + ".raw";
    info(info.goodData==0, :) = [];
    for k = 1:length(badNames)
        data(strcmp(badNames(k),data.FileName), :) = [];
    end
    clear k
end
clear badNames pruneData

% How many possible standards are there, and where will they be?
% KL 10/20/2016 now need to look for positive or negative set bc UPLC data
% has two standard curves
switch ionMode
    case 'neg'
        kStandard = (strcmp('std', info.sType) & strcmp("neg", info.ionMode));           
    case 'pos'
        kStandard = (strcmp('std', info.sType) & strcmp("pos", info.ionMode));           
end

% Set up indices and names to point to the standards and samples.
standardNames = info.FileName(kStandard);
nStandards = sum(kStandard);
kSample = strcmp('Unknown',info.SampleType);
sampleNames = info.FileName(kSample);
nSamples =sum(kSample);
clear kStandard kSample

% Set up goodData, which is (numCompounds x numSamples)
goodData(1:length(compoundList.names),nSamples) = NaN; 
goodDataError = goodData;

% Add variables to the compound list concerning standard curve and goodness
% of fit.
warning('off', 'stats:dataset:subsasgn:DefaultValuesAddedVariable');
compoundList.r2_line = zeros(length(compoundList.names),1);
compoundList.slope = zeros(length(compoundList.names),1);
compoundList.intercept = zeros(length(compoundList.names),1);
compoundList.SDslope = zeros(length(compoundList.names),1);
compoundList.SDintercept = zeros(length(compoundList.names),1);
compoundList.nPoints = zeros(length(compoundList.names),1);
compoundList.LOD = zeros(length(compoundList.names),1);
compoundList.LOQ = zeros(length(compoundList.names),1);


% Require five points for calibration (set lower for now because of small
% curves).
nRequired = 5;

% This PDF will contain graphs of the standard curves as fitted, as well as
% prediction intervals. 
% NPG 19 Sept 2023: Really should redo these sections so that it doesn't
% stop the code completely if the file exists. 
curveFilename = append(ionMode,"_",SILISType, "_mtabs_stdCurves.pdf");
if exist(append(oFolder, filesep, curveFilename), 'file') == 2
    fprintf('The file "%s" already exists in the working directory.\n', curveFilename);
    fprintf('Please delete the file before proceeding.\n');
    error('PDF file already exists. Code execution halted.');
else
    % Create the diary file or perform other operations
    fprintf('Calibration curve PDF file "%s" created.\n', curveFilename)

% Go through one compound at a time, (1) make the standard curve, (2) use that to
% calculate the concentrations for each sample.
for a = 1:length(compoundList.names)
    
    clear xdata ydata
    
    % This index pulls all observations of a molecule or its SIL-IS as
    % specified in the inputs. 
    k = (strcmp(compoundList.names(a),data.MoleculeName) &...
        (data.precursor==1 | string(data.IsotopeLabelType) == SILISType));
    
    % Call the rest only if the molecule is actually in the exported data.
    if sum(k)>0
        % smallDS is the small dataset extracted for just this molecule.
        smallDS = data(k,:);
        clear k
        
        % Match record numbers and unite the sample type column without
        % merging the two datasets.
        [~, ia, ib] =intersect(smallDS.FileName,[info.FileName+".raw"]);
        smallDS.sType(ia,1) = info.sType(ib,1);
        clear c ia ib  
        
        % Calculating light-heavy ratios
        smallDS.LHR = zeros(height(smallDS),1);
        
        % Hey, wait, do we need this if statement?
        if nSILIS == 2
            heavyArea = smallDS.Area(smallDS.IsotopeLabelType==...
            convertCharsToStrings(SILISType));
        elseif nSILIS == 1
            heavyArea = smallDS.Area(smallDS.IsotopeLabelType==convertCharsToStrings(SILISType));
        end
        
        % Find light peak areas and take the ratios. You can get Skyline to
        % do this for you, but we do it in the code here.
        lightArea = smallDS.Area(smallDS.IsotopeLabelType=="light");
        if length(heavyArea) == length(lightArea)
            smallDS.LHR(smallDS.IsotopeLabelType=="light") = ...
                lightArea./heavyArea;
            smallDS.LHR(isinf(smallDS.LHR))=NaN;
        else 
            disp(['There was a mismatch in the number of'...
                ' light and heavy ion measurements for '...
                compoundList.names{a}])
            continue
        end
        
        % Now that we've calculated those ratios, I'm going to delete the
        % heavy measurements. 
        smallDS(string(smallDS.IsotopeLabelType)==SILISType,:)=[];
        
        [~, idxDS, idxStandards] = intersect(smallDS.FileName,...
            [standardNames + ".raw"]);
        clear c

        % Set standard curve analyte concentration values based on what was
        % entered in Skyline
        setStandardConcentrations = smallDS.AnalyteConcentration;
        setStandardConcentrations(isnan(setStandardConcentrations))=[];      
        xdata = setStandardConcentrations;
        ydata(1:length(xdata),1) = NaN;
        
        % Get all possible values from the standard curve
        ydata(idxStandards) = smallDS.LHR(idxDS);
        xdata = cat(1,xdata); 
        ydata = cat(1,ydata); 
        clear idxDS idxStandards 
        
        % Remember, will also have cases where no data were found for 
        % select samples so need to setup the spacers in there.
        % "tData" is basically the data to be fed into the calibrated
        % equation. 
        [~, ia, ib] = intersect(smallDS.FileName,[sampleNames+".raw"]);
        tData(1:length(sampleNames),1) = NaN;
        tData(ib) = smallDS.LHR(ia);
        clear c ia ib
        
        su = strcmp(info.sType,'rep');
        ksu = find(su==1);
        [c, ~, ib] = intersect([info.FileName(ksu)+".raw"],smallDS.FileName);
        % This is probably redundant but checks to see if the export was 
        % blank. "unknownsOnly" omits pools.
        if ~isempty(c)
            tData_unknownsOnly = smallDS.LHR(ib);
        else
            tData_unknownsOnly = NaN;
        end
        clear su ksu c ib
        
        % We need to make the calibration curve as close to the range of
        % the data, and the first step is finding what the LHR of our max
        % sample is. 
        m = nanmax(tData_unknownsOnly);
        % If all the unknowns fail the quality check, this next step
        % will fail. Haven't seen this until now (6/26/2018)
        if isnan(m)
            %easiest to make kMax empty
            kMax = [];
        else
            kMax = find(ydata <= m);
        end
        clear tData_unknownsOnly

diary on    % Record what happens here, for example if something is above 
            % the standard curve in one or more samples (happens more than
            % you'd think).

        if ~isempty(kMax)
            % Have at least one point on the curve.
            if isequal(kMax(end),length(ydata))
                % Already at the end of the standard curve...so use all the points
                % do nothing...but send up a flag since the data are above
                % the standard curve.
                disp([compoundList.names{a} ' is above the standard curve'])
            elseif isequal(kMax(end)+1,length(ydata))
                % Only one more above the points in the standard curve, use all the
                % points.
            elseif isequal(kMax,1)
                % Data are at the low end of the standard curve, but let's require
                % more points above my data to get a reasonable curve...
                xdata = xdata(1:nRequired);
                ydata = ydata(1:nRequired);
            elseif length(kMax)+2  < nRequired
                % Use the number of points sent in nRequired.
                ydata = ydata(1:nRequired);
                xdata = xdata(1:nRequired);
            elseif length(kMax) + 1 < nRequired
                % Use the standard curve to one point beyond the range of my
                % samples
                ydata = ydata(1:kMax(end)+1);
                xdata = xdata(1:kMax(end)+1);
            else
                % Use the standard curve to one point beyond the range of my
                % samples
                ydata = ydata(1:kMax(end)+1);
                xdata = xdata(1:kMax(end)+1);
                
            end
        elseif isempty(kMax)
            %all of the points in the standard curve are higher than what was
            %measured in the samples
            ydata = ydata(1:nRequired);
            xdata = xdata(1:nRequired);
        end
        clear kMax m
       
diary off 

clear diaryFilename

        % Need at least three points to make a curve AND get the error 
        % estimates. These next couple lines actually cause a lot of
        % errors. Usually, this a result of the sequence file being
        % formatted wrong such that xdata and ydata don't end up the same
        % length.
        try
            show = [xdata ydata];
        catch
            fprintf('here')
        end
        % For whatever reason we screen NaNs out again. There probably are
        % none by this point. 
        i = isnan(show);
        sfmi = sum(i,2);
        k = find(sfmi==0);
        xdata = xdata(k);
        ydata = ydata(k);
        clear show i sfmi k

        % This will be helpful bc will show where I had <2 points.
        % Remember that this also takes into account the rules I set above about
        % how wide to make the standard curve.
        compoundList.nPoints(a) = length(ydata);
        
        if length(xdata)>2 % Don't make a line with fewer than 3 points.
            
            % See getErrors subfunction at end. This and useErrors together
            % perform a regression, get goodness-of-fit, and then use it to
            % calibrate unknown concentrations.
            dataOut = getErrors(xdata,ydata);
            [calcError, calcConc] = useErrors(dataOut,tData); %then calculate the concentrations


            % If the slope is negative, this is garbage. If it's good,
            % store information about both calculated concentrations and
            % the cal curve in memory. 
            if dataOut.slope > 0
                goodData(a,:) = calcConc;
                goodDataError(a,:) = calcError; %can get percent by calcError./calcConc

                compoundList.slope(a) = dataOut.slope;
                compoundList.intercept(a) = dataOut.intercept;
                compoundList.SDslope(a) = dataOut.SDslope;
                compoundList.SDintercept(a) = dataOut.SDintercept;
                compoundList.r2_line(a) = dataOut.r2;
                compoundList.LOD(a) = 3.3*(dataOut.SDintercept./dataOut.slope);
                compoundList.LOQ(a) = 10*(dataOut.SDintercept./dataOut.slope);


                if 1 % You can turn off the plotting if you like.
                    % This plots the standard curve. Ironically, it refits it
                    % completely using MATLAB's built in functions.
                    set(groot,'defaultFigureVisible','off')
                    figure
                    plot(fitlm(xdata,ydata));
                    hold on
                    plot(calcConc,tData,'ok','DisplayName','Samples')
                    title(string(compoundList.names{a}) + " " + string(ionMode) + " " + string(SILISType))
                    xlabel(append('Standard Concentration Added (', units, ")"))
                    ylabel('Peak Ratio (light/heavy)')
                    text(.95,.97, "R^2 =" + string(dataOut.r2),'Units','normalized')
                    xline(compoundList.LOD(a),'--g','LOD','DisplayName','LOD')
                    xline(compoundList.LOQ(a),'--b','LOQ','DisplayName','LOQ')
                    exportgraphics(gca, curveFilename, 'Append',  true)
                    hold off
                    close(gcf)
                    set(groot,'defaultFigureVisible','on')

                end

            else
                % slope is less than zero
                goodData(a,:) = NaN;
                goodDataError(a,:) = NaN;
                compoundList.slope(a) = NaN;
                compoundList.intercept(a) = NaN;
                compoundList.SDslope(a) = NaN;
                compoundList.SDintercept(a) = NaN;
                compoundList.r2_line(a) = NaN;
                compoundList.LOD(a) = NaN;
                compoundList.LOQ(a) = NaN;

            end
            
        else
            % not enough points to make a standard curve
            goodData(a,:) = NaN;
            goodDataError(a,:) = NaN;
            compoundList.slope(a) = NaN;
            compoundList.intercept(a) = NaN;
            compoundList.SDslope(a) = NaN;
            compoundList.SDintercept(a) = NaN;
            compoundList.r2_line(a) = NaN;
            compoundList.LOD(a) = NaN;
            compoundList.LOQ(a) = NaN;


        end
        clear dataOut calcError calcConc tData
        
    end 
    
end

end

clear a compound xdata ydata smallDS

% Replace values less than the calculated LOD with NaN         
goodData_filtered = goodData;

for a = 1:length(compoundList.names)
    k = find(goodData_filtered(a,:)<= compoundList.LOD(a));
    goodData_filtered(a,k) = NaN;
    
    clear k
end 

ds1 = table(goodData);
ds2 = table(goodData_filtered);
ds3 = table(goodDataError);
keepingAll = cat(2,compoundList,ds1,ds2,ds3);
clear ds1 ds2 ds3 compoundList goodData goodDataError

% %remove some unneeded variables:
% keepingAll.indexMain = [];
% keepingAll.indexConfirm = [];

%perhaps do a little pruning to provide a dataset array with only the data
%that are good from the criteria above and are not overlapping with
%zero...
i = isnan(keepingAll.SDintercept);
k = find(i~=1);

keepGoodData = keepingAll(k,:);
%here, we need to consider on a sample by sample basis and not make
%decisions based on the entire set for each compound

for a = 1:size(keepGoodData,1)
    for aa = 1:size(keepGoodData.goodData,2)
        tD = keepGoodData.goodData(a,aa);
        tE = keepGoodData.goodData_filtered(a,aa);
        tF = keepGoodData.goodDataError(a,aa);
        
        %can have a few options.
        if tD < 0 %easiest: sample is less than zero
            keepGoodData.goodData(a,aa)=0;
            keepGoodData.goodData_filtered(a,aa)=0;
            keepGoodData.goodDataError(a,aa) = 0;
        elseif tD - tF < 0 %does the error window include zero?
            %what is the window around it,
            keepGoodData.goodData(a,aa)=0;
            keepGoodData.goodData_filtered(a,aa)=0;
            keepGoodData.goodDataError(a,aa) = 0;
            %             elseif tE./tD*100 > 66 %is the error percent above 66%?
            %                 %added 5/18/2016 bc getting too many things that have
            %                 %values that get calculated, but the error is high cfd
            %                 %to the measured value
            %                 keepGoodData.goodData(a,aa)=0;
            %                 keepGoodData.goodDataError(a,aa) = 0;
            
        end
        clear tD tE TF aa
    end
end

% If a column ends up being a mix of zeros and NaNs, set all to zero (maybe
% not the best practice.
% Note, this with the following section could be turned into a single loop.
i = isnan(keepGoodData.goodData);
for a = 1:size(i,1)
    td = keepGoodData.goodData(a,:);
    ts = sum(td(i(a,:)~=1));
    if ts==0
        keepGoodData.goodData(a,i(a,:)==1) = 0;
    end
     clear td ki k ts
end
 clear a i

%repeat for goodData_filtered
 i = isnan(keepGoodData.goodData_filtered);
for a = 1:size(i,1)
    td = keepGoodData.goodData_filtered(a,:);
    ts = sum(td(i(a,:)~=1));
    if ts==0
        keepGoodData.goodData_filtered(a,i(a,:)==1) = 0;
    end
     clear td ki k ts
end
 clear a i

% Now ahead and *delete* the rows where all the datapoints are zero...this
% assumes that the user is familiar with the list of compounds and knows
% about compounds that could have been in the samples but were all zero.
fm = logical(keepGoodData.goodData~=0);
sfmc = sum(fm,2);
k = find(sfmc==0);

keepGoodData([k],:) = [];
     
clear fm sfmc k setStandardConcentrations


end

% put the internal functions here at the end %
% These functions were implemented and modified in considerMAVEN.m by K.
% Longnecker

    function dataOut = getErrors(xdata,ydata)
        %function dataOut = getErrors(xdata,ydata)
        %From this web site:
        %http://terpconnect.umd.edu/~toh/spectrum/LeastSquaresMatlab.txt
        %KL modifying 4/21/2014
        
        x = xdata;
        y = ydata;
        % Simple Matlab script for calculating the first-order least-square fit of y vs x,
        % including the Slope and Intercept and the predicted standard deviation of
        % the slope (SDSlope) and intercept (SDIntercept).
        
        NumPoints=length(x);
        Sxx = sum((x-mean(x)).^2);
        Syy = sum((y-mean(y)).^2);
        Sxy = sum((x-mean(x)).*(y-mean(y)));
        Slope = Sxy./Sxx;
        Intercept = mean(y)-Slope*mean(x);
        
        Sy = sqrt((Syy-Slope^2*Sxx)/(NumPoints-2));
        
        SDslope = Sy/sqrt(Sxx);
        SDintercept = Sy*sqrt(1./(NumPoints-(sum(x).^2)./sum(x.^2)));
        
        r2 = 1 - ((Syy-Slope^2*Sxx) ./Syy);
        
        %data to send out of this function (when it is a function)
        dataOut.slope = Slope;
        dataOut.intercept = Intercept;
        dataOut.SDslope = SDslope;
        dataOut.PercentSlopeError = SDslope./Slope;
        dataOut.SDintercept = SDintercept;
        dataOut.PercentInterceptError = SDintercept./Intercept;
        dataOut.r2 = r2;
        
    end %end of getErrors as a function


    function [calcError, calcConc] = useErrors(myErrorData,measuredSample)
        %function [calcError, calcConc] = useErrors(Slope,Intercept,SDslope,SDintercept,measuredSample)
        %use the errors on the line to get the errors on the samples actually
        %measured
        %KL 4/21/2014
        Intercept = myErrorData.intercept;
        Slope = myErrorData.slope;
        SDintercept = myErrorData.SDintercept;
        SDslope = myErrorData.SDslope;
        
        %calculated concentrations from my hypothetical list
        calcConc = (measuredSample - Intercept)./Slope;
        
        %apply to my list of hypothetical unknowns, split this up to make
        %it easier to keep track of where the parentheses etc. are
        fSQ = (SDintercept./(measuredSample - Intercept)).^2 + (SDslope./Slope).^2;
        calcError = calcConc .* sqrt(fSQ);
        errorPercent = calcError./calcConc*100;
        
    end %end of useErrors as a function