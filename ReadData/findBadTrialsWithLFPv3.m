% We check, for each trial, whether there is a value greater than 
% 1. threshold times the std dev. 
% 2. maxLimit
% If yes, that trial is marked as a bad trial

% Additions/changes made on 09 March 2024 by Surya S P

% Bad trials can be now be saved for each array separately in case of dual
% array and EEG recordings by entering the electrode numbers of an array as
% 'checkTheseElectrodes' variable (eg. [1:48] or [49:96]) and setting the arrayString
% variable indicating the area in which the array is placed (eg V1 or V4 or EEG).

% Option to override the marginals threshold by setting marginalFlag to 1.Default is 0. 
% Set this flag to 1 only if rejectTolerence variable is non-empty.

% Bad trials file will be saved as Badtrials[ArrayString].mat.

% If no arrayString is given and checkTheseElectrodes include all the electrodes then badTrials
% common across all the electrodes are computed like in the previous version.

% A minor bug where badElecsMarginalStats is computed is fixed. 

% Summary figure now displays the bad trials for 'CheckTheseElectrodes'
% only.

% Option is added to remove bad trials based on PSD for which checkPSDFlag should be set to 1 and 
% separate threshold (psd_threshold) and frequency range (checkFrequencyRange) should be provided.
% By default psd_threshold = 6 and checkFrequencyRange = [120 200] Hz;

function [allBadTrials,badTrials] = findBadTrialsWithLFPv3(monkeyName,expDate,protocolName,folderSourceString,gridType,checkTheseElectrodes,processAllElectrodes,threshold,maxLimit,minLimit,showElectrodes,saveDataFlag,checkPeriod,rejectTolerance,marginalsFlag,arrayString,checkPSDFlag)

if ~exist('checkTheseElectrodes','var');     checkTheseElectrodes = [33 12 80 63 44];   end
if ~exist('processAllElectrodes','var');     processAllElectrodes = 0;                  end
if ~exist('folderSourceString','var');       folderSourceString = 'G:';                 end
if ~exist('threshold','var');                threshold = 6;                             end
if ~exist('minLimit','var');                 minLimit = -2000;                          end
if ~exist('maxLimit','var');                 maxLimit = 1000;                           end
if ~exist('saveDataFlag','var');             saveDataFlag = 1;                          end
if ~exist('checkPeriod','var');              checkPeriod = [-0.7 0.8];                  end
if ~exist('rejectTolerance','var');          rejectTolerance = 1;                       end
if ~exist('showElectrodes','var');           showElectrodes = [];                       end
if ~exist('marginalsFlag','var');            marginalsFlag = 0;                         end
if ~exist('arrayString','var');              arrayString = [];                          end
if ~exist('checkPSDFlag','var');             checkPSDFlag = 0;                          end

folderName = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName);
folderSegment = fullfile(folderName,'segmentedData');

load(fullfile(folderSegment,'LFP','lfpInfo.mat'));

if processAllElectrodes % compute bad trials for all the saved electrodes
    numElectrodes = length(analogChannelsStored);
else % compute bad trials for only the electrodes mentioned
    numElectrodes = length(checkTheseElectrodes);
end

allBadTrials = cell(1,numElectrodes);
nameElec = cell(1,numElectrodes);

for i=1:numElectrodes
    if processAllElectrodes
        electrodeNum=analogChannelsStored(i); % changed from checkTheseElectrodes
        % to calculate bad trials for each electrode irrespective of the
        % electrodes to be checked
    else
        electrodeNum=checkTheseElectrodes(i);
    end
    load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeNum) '.mat']));
    
    disp(['Processing electrode: ' num2str(electrodeNum)]);
    nameElec{i} = ['elec' num2str(electrodeNum)];
    disp(nameElec{i});
    
    analogDataSegment = analogData;
    % determine indices corresponding to the check period
    checkPeriodIndices = timeVals>=checkPeriod(1) & timeVals<=checkPeriod(2);
    
    analogData = analogData(:,checkPeriodIndices);
    % subtract dc
    analogData = analogData - repmat(mean(analogData,2),1,size(analogData,2));
    
    numTrials = size(analogData,1); %#ok<*NODEF>
    meanData = mean(analogData,2)';
    stdData  = std(analogData,[],2)';
    maxData  = max(analogData,[],2)';
    minData  = min(analogData,[],2)';
    
    clear tmpBadTrials tmpBadTrials2
    tmpBadTrials = unique([find(maxData > meanData + threshold * stdData) find(minData < meanData - threshold * stdData)]);
    tmpBadTrials2 = unique(find(maxData > maxLimit));
    tmpBadTrials3 = unique(find(minData < minLimit));
    if checkPSDFlag  % borrowed from findBadTrialsWithEEG
        checkFrequencyRange = [120 200];
        psd_threshold = 6;
        
        Fs = 1/(timeVals(2) - timeVals(1)); %Hz
        params.tapers   = [3 5];
        params.pad      = -1;
        params.Fs       = Fs;
        params.fpass    = [0 200];
        params.trialave = 0;

        clear powerVsFreq;
        [powerVsFreq,freqVals] = mtspectrumc(analogData',params);
        checkFreqPos = intersect(find(freqVals>checkFrequencyRange(1)),find(freqVals<checkFrequencyRange(2)));
        powerVsFreq = powerVsFreq(checkFreqPos,:)';

        clear meanTrialData stdTrialData tDplus
        meanTrialData = mean(powerVsFreq,1,'omitnan');  % mean trial trace
        stdTrialData = std(powerVsFreq,[],1,'omitnan'); % std across trials

        tDplus = (meanTrialData + (psd_threshold)*stdTrialData);    % upper boundary/criterion
        tDMinus = (meanTrialData - (psd_threshold)*stdTrialData);    % lower boundary/criterion
        clear tBoolTrialsPlus tmpBadTrials4
        tBoolTrialsPlus = sum((powerVsFreq > ones(numTrials,1)*tDplus),2);
        tBoolTrialsMinus = sum((powerVsFreq < ones(numTrials,1)*tDMinus),2);
        tmpBadTrials4 = (find(tBoolTrialsPlus>0))';
        tmpBadTrials5 = (find(tBoolTrialsMinus>0))';
    else
        checkFrequencyRange = [];
        psd_threshold = [];
        tmpBadTrials4 = [];
        tmpBadTrials5 = [];
    end
    allBadTrials{electrodeNum} = unique([tmpBadTrials tmpBadTrials2 tmpBadTrials3 tmpBadTrials4 tmpBadTrials5]);
end

numElectrodes = length(checkTheseElectrodes); % check the list for these electrodes only to generate the overall badTrials list
j = checkTheseElectrodes(1);
badTrials=allBadTrials{j};
for i=1:numElectrodes
    j = checkTheseElectrodes(i);
    badTrials=intersect(badTrials,allBadTrials{j}); % in the previous case we took the union
end

%**************************************************************************
% [vinay] SOME ADDITIONAL CRITERIA for detecting bad repeats and electrodes
%--------------------------------------------------------------------------
% [Vinay] - decide as per a tolerance for the percent of electrodes showing
% a particular stimulus as bad
if exist('rejectTolerance','var')
    badTrials = [];
    for n=1:numTrials
        trialCount=0;
        for i=1:numElectrodes
            j = checkTheseElectrodes(i);
            if ~isempty(find(allBadTrials{j} == n, 1))
                trialCount=trialCount+1;
            end
        end
        trialPercent = trialCount/numElectrodes;
        if trialPercent>=rejectTolerance
            badTrials = cat(2,badTrials,n);
        end
    end
end
%-----

% [vinay] - statistically determine the repeats which are flagged as bad
% across a significant number of electrodes and reject them. Similarly tag
% those electrodes as bad which show a significantly higher number of
% badTrials as compared to the mean statistics across electrodes.

% look at the overall distribution of badTrials across electrodes
% and reject the repeats which are flagged as bad for >= (mean+std*threshold) 
% number of electrodes. (across trials)
% Similarly tag those electrodes as bad which show number of badTrials 
% >= (mean+std*threshold) (across electrodes)
allBadTrialsMatrix = zeros(length(allBadTrials),numTrials);
for i=1:length(allBadTrials)
    allBadTrialsMatrix(i,allBadTrials{i}) = 1;
end
marginalStim = sum(allBadTrialsMatrix(checkTheseElectrodes,:),1); % for each stimulus - no. of electrodes showing this stimulus as bad
marginalElectrodes = sum(allBadTrialsMatrix,2); % for each electrode - no. of stimuli flagged as bad

if processAllElectrodes && strcmpi(gridType,'Microelectrode')
    if strcmpi(monkeyName,'tutu')
        electrodesForMarginals = 1:81;
    else
        electrodesForMarginals = 1:96;
    end
else
    electrodesForMarginals = checkTheseElectrodes;
end

thresholdMarginal = threshold; % thresholdMarginal i.e. values> u+thresholdMarginal*sigma will be tagged as bad 
badTrialsMarginalStats = find(marginalStim>=(mean(marginalStim)+std(marginalStim)*thresholdMarginal));
badElecsMarginalStatsPos = find(marginalElectrodes(electrodesForMarginals)>=(mean(marginalElectrodes(electrodesForMarginals))+std(marginalElectrodes(electrodesForMarginals))*thresholdMarginal));

if marginalsFlag
    marginalsThreshold = mean(marginalStim)+std(marginalStim)*thresholdMarginal;
    rejectToleranceThreshold = length(checkTheseElectrodes)*rejectTolerance;
    if marginalsThreshold<rejectToleranceThreshold
        badTrials = intersect(badTrials,badTrialsMarginalStats);
        if iscolumn(badTrials)
            badTrials = badTrials';
        end
    else
        badTrials = unique(cat(2,badTrials,badTrialsMarginalStats));
    end
else
    badTrials = unique(cat(2,badTrials,badTrialsMarginalStats));
end
badElecs = electrodesForMarginals(badElecsMarginalStatsPos)'; %#ok<FNDSB> 

disp(['total Trials: ' num2str(numTrials) ', bad trials indices: ' num2str(badTrials)]);
disp(['total Elecs: ' num2str(length(marginalElectrodes)) ', bad elecs indices: ' num2str(badElecs')]);
%--------------------------------

for i=1:numElectrodes
    j = checkTheseElectrodes(i);
    if length(allBadTrials{j}) ~= length(badTrials)
        disp(['Bad trials for electrode ' num2str(checkTheseElectrodes(i)) ': ' num2str(length(allBadTrials{j}))]);
    else
        disp(['Bad trials for electrode ' num2str(checkTheseElectrodes(i)) ': common bad trials only (' num2str(length(badTrials)) ')']);
    end
end

if saveDataFlag
    disp(['Saving ' num2str(length(badTrials)) ' bad trials']);
    save(fullfile(folderSegment,['badTrials' arrayString '.mat']),'badTrials','checkTheseElectrodes','threshold','maxLimit','minLimit','checkPeriod','allBadTrials','nameElec','rejectTolerance','badElecs','badTrialsMarginalStats','psd_threshold','checkFrequencyRange');
else
    disp('Bad trials will not be saved..');
end

load(fullfile(folderSegment,'LFP','lfpInfo.mat'));

lengthShowElectrodes = length(showElectrodes);
if ~isempty(showElectrodes)
    for i=1:lengthShowElectrodes
        figure;
        subplot(2,1,1);
        channelNum = showElectrodes(i);

        clear signal analogData analogDataSegment
        load(fullfile(folderSegment,'LFP',['elec' num2str(channelNum) '.mat']));
        analogDataSegment = analogData;
        if numTrials<4000
            plot(timeVals,analogDataSegment(setdiff(1:numTrials,badTrials),:),'color','k');
            hold on;
        else
            disp('More than 4000 trials...');
        end
        if ~isempty(badTrials)
            plot(timeVals,analogDataSegment(badTrials,:),'color','g');
        end
        title(['electrode ' num2str(channelNum)]);
        axis tight;

        subplot(2,1,2);
        plot(timeVals,analogDataSegment(setdiff(1:numTrials,badTrials),:),'color','k');
        hold on;
        j = channelNum;
        if ~isempty(allBadTrials{j})
            plot(timeVals,analogDataSegment(allBadTrials{j},:),'color','r');
        end
        axis tight;
    end
end

%**************************************************************************
% summary plot
%--------------------------------------------------------------------------
% allBadTrialsMatrix = zeros(length(allBadTrials),numTrials);
% for i=1:length(allBadTrials)
%     allBadTrialsMatrix(i,allBadTrials{i}) = 1;
% end

summaryFig = figure('name',[monkeyName expDate protocolName],'numbertitle','off');
h0 = subplot('position',[0.8 0.8 0.18 0.18]); set(h0,'visible','off');
text(0.05, 0.7, ['thresholds (uV): [' num2str(minLimit) ' ' num2str(maxLimit) ']'],'fontsize',12,'unit','normalized','parent',h0);
text(0.05, 0.4, ['checkPeriod (s): [' num2str(checkPeriod(1)) ' ' num2str(checkPeriod(2)) ']'],'fontsize',12,'unit','normalized','parent',h0);
text(0.05, 0.1, ['rejectTolerance : ' num2str(rejectTolerance)],'fontsize',12,'unit','normalized','parent',h0);

h1 = getPlotHandles(1,1,[0.07 0.07 0.7 0.7]);
subplot(h1);
imagesc(1:numTrials,length(allBadTrials):-1:1,flipud(allBadTrialsMatrix),'parent',h1);
set(gca,'YDir','normal','ylim',[checkTheseElectrodes(1) checkTheseElectrodes(end)]); colormap(gray);
xlabel('# trial num','fontsize',15,'fontweight','bold');
ylabel('# electrode num','fontsize',15,'fontweight','bold');

h2 = getPlotHandles(1,1,[0.07 0.8 0.7 0.17]);
h3 = getPlotHandles(1,1,[0.8 0.07 0.18 0.7]);
subplot(h2); cla; set(h2,'nextplot','add');
stem(h2,1:numTrials,sum(allBadTrialsMatrix(checkTheseElectrodes,:),1)); axis('tight');
ylabel('#count');
if ~isempty(badTrials)
    stem(h2,badTrials,sum(allBadTrialsMatrix(checkTheseElectrodes,badTrials),1),'color','r');
end
subplot(h3); cla; set(h3,'nextplot','add');
stem(h3,1:length(allBadTrials),sum(allBadTrialsMatrix,2)); axis('tight'); ylabel('#count');
if ~isempty(badElecs)
    stem(h3,badElecs,sum(allBadTrialsMatrix(badElecs,:),2),'color','r');
end
xlim(h3,[checkTheseElectrodes(1) checkTheseElectrodes(end)]);
view([90 -90]);

saveas(summaryFig,fullfile(folderSegment,[monkeyName expDate protocolName 'summmaryBadTrials' arrayString '.fig']),'fig');
%saveas(summaryFig,[monkeyName expDate protocolName 'summmaryBadTrials' arrayString '.fig'],'fig');

%**************************************************************************
% fill the badTrials summary sheet
%**************************************************************************




end
