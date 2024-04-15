% This program combines and saves the bad trials across different arrays like V1, V4, EEG as well as
% saves the common bad trials for each array separately if not saved already.

% Description of some of the input variables:

% 1) arrayStringList is a cellarray containing the name of arrays. Eg: {V1,V4,EEG}
% 2) checkTheseElectrodes is a cellarray with each cell containing the electrode numbers of each 
% array specified in the same order as in the arrayStringList. Eg: {1:48,49:96,97:111}
% 3) showElectrodes is a cellarray with each cell containing the electrodes from respective arrays 
% for which the good and bad trial LFP traces will be displayed. It should either have the same dimension 
% as arrayStringList or should be an empty cell if no electrodes to be shown in any array
% 4) MinLimit, MaxLimit, threshold and rejectTolerence can each be passed as an array of same 
% length as arrayStringList if they are different for each array or as a scalar if they are same 
% for all the arrays.
% 5) CheckPeriod is a cellarray - a provison made to run the bad trials for
% multiple checkPeriods. CheckPeriod should either be a single cell or cell
% of length same as arrayStringList.

% Surya S P 12 March 2024
function [allBadTrials,badTrials] = combineBadTrialsWithLFP(monkeyName,expDate,protocolName,folderSourceString,gridType,saveDataFlag,arrayStringList,checkTheseElectrodes,processAllElectrodes,threshold,maxLimit,minLimit,showElectrodes,checkPeriod,rejectTolerance,marginalsFlag,checkPSDFlag)

if ~exist('folderSourceString','var');       folderSourceString = 'G:';                 end
if ~exist('gridType','var');                 gridType = 'Microelectrode';               end
if ~exist('saveDataFlag','var');             saveDataFlag = 1;                          end
if ~exist('arrayStringList','var');          arrayStringList = {'V1','V4'};             end
if ~exist('checkTheseElectrodes','var');     checkTheseElectrodes = {1:48,49:96};       end
if ~exist('threshold','var');                threshold = 6;                             end
if ~exist('minLimit','var');                 minLimit = -2000;                          end
if ~exist('maxLimit','var');                 maxLimit = 1000;                           end
if ~exist('rejectTolerance','var');          rejectTolerance = 1;                       end
if ~exist('processAllElectrodes','var');     processAllElectrodes = 0;                  end
if ~exist('checkPeriod','var');              checkPeriod = {[-0.7 0.8]};                end
if ~exist('showElectrodes','var');           showElectrodes = {};                       end
if ~exist('marginalsFlag','var');            marginalsFlag = 0;                         end
if ~exist('checkPSDFlag','var');             checkPSDFlag = 0;                          end
folderSegment = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData');

if isscalar(minLimit)
    minLimit = repmat(minLimit,1,length(arrayStringList));
end
if isscalar(maxLimit)
    maxLimit = repmat(maxLimit,1,length(arrayStringList));
end
if isscalar(threshold)
    threshold = repmat(threshold,1,length(arrayStringList));
end
if isscalar(rejectTolerance)
    rejectTolerance = repmat(rejectTolerance,1,length(arrayStringList));
end
if isempty(showElectrodes)
    showElectrodes = repmat({[]},1,length(arrayStringList));
end
if isscalar(checkPeriod)
    checkPeriod = repmat(checkPeriod,1,length(arrayStringList));
end

%%%%%%%%%%% Loading bad trials of individual arrays %%%%%%%%%%%%%%

for i=1:length(arrayStringList)
    fileName = fullfile(folderSegment,['badTrials' arrayStringList{i} '.mat']);

    if ~exist(fileName,'file')
        disp(['Saving the bad trial data for ' arrayStringList{i}]);

        findBadTrialsWithLFPv3(monkeyName,expDate,protocolName,folderSourceString,gridType,checkTheseElectrodes{i},processAllElectrodes,threshold(i),maxLimit(i),minLimit(i),showElectrodes{i},saveDataFlag,checkPeriod{i},rejectTolerance(i),marginalsFlag,arrayStringList{i},checkPSDFlag);
    end

    badTrialsData(i) = load(fileName); %#ok<*AGROW> 
end

%%%%%%%%% Combining the bad trials across arrays %%%%%%%%%%%%
clear threshold maxLimit minLimit checkPeriod rejectTolerance checkTheseElectrodes
badTrialsTMP = [];
allBadTrialsTMP = [];
allBadTrials = [];
badTrialsMarginalStatsTMP = [];
checkTheseElectrodes = [];
badElecs = [];
nameElec = {};
checkFrequencyRange = {};
psd_threshold = [];
for i=1:length(arrayStringList)
    badTrialsTMP = cat(2,badTrialsTMP,badTrialsData(i).badTrials);
    allBadTrialsTMP{i} = badTrialsData(i).allBadTrials(badTrialsData(i).checkTheseElectrodes);
    badTrialsMarginalStatsTMP = cat(2,badTrialsMarginalStatsTMP,badTrialsData(i).badTrialsMarginalStats);
    threshold(i) = badTrialsData(i).threshold;
    maxLimit(i) = badTrialsData(i).maxLimit;
    minLimit(i) = badTrialsData(i).minLimit;
    checkPeriod(i,:) = badTrialsData(i).checkPeriod;
    rejectTolerance(i) = badTrialsData(i).rejectTolerance;
    nameElec = cat(2,nameElec,badTrialsData(i).nameElec);
    badElecs = cat(1,badElecs,badTrialsData(i).badElecs);
    checkTheseElectrodes = cat(2,checkTheseElectrodes,badTrialsData(i).checkTheseElectrodes);
    checkFrequencyRange = cat(2,checkFrequencyRange,badTrialsData(i).checkFrequencyRange);
    psd_threshold = cat(2,psd_threshold,badTrialsData(i).psd_threshold);
end
badTrials = unique(badTrialsTMP);
badTrialsMarginalStats = unique(badTrialsMarginalStatsTMP);
badElecs = sort(unique(badElecs));
allBadTrials = cat(2,allBadTrials,allBadTrialsTMP{:});

if saveDataFlag
    disp(['Saving ' num2str(length(badTrials)) ' combined bad trials']);
    save(fullfile(folderSegment,'badTrials.mat'),'badTrials','checkTheseElectrodes','threshold','maxLimit','minLimit','checkPeriod','allBadTrials','nameElec','rejectTolerance','badElecs','badTrialsMarginalStats','arrayStringList','checkFrequencyRange','psd_threshold');
else
    disp('Bad trials will not be saved..');
end

%**************************************************************************
% summary plot
%--------------------------------------------------------------------------

load(fullfile(folderSegment,'LFP',[nameElec{1} '.mat'])); %#ok<LOAD> 
numTrials = size(analogData,1);
for j=1:length(arrayStringList)
    badTrialsMatrix{j} = zeros(max(checkTheseElectrodes),numTrials);
    for i=1:length(allBadTrialsTMP{j})
        badTrialsMatrix{j}(badTrialsData(j).checkTheseElectrodes(i),allBadTrialsTMP{j}{i}) = 1;
    end
end
allBadTrialsMatrix = zeros(length(allBadTrials),numTrials);
for i=1:length(allBadTrials)
    allBadTrialsMatrix(i,allBadTrials{i}) = 1;
end
for i=1:length(arrayStringList)
    arrayList{i} = [arrayStringList{i} ';  '];
    elecNumString{i} = [' [' num2str(badTrialsData(i).checkTheseElectrodes(1)) ' ' num2str(badTrialsData(i).checkTheseElectrodes(end)) ']; '];
    thresholdStringList{i} = [' [' num2str(minLimit(i)) ' ' num2str(maxLimit(i)) ']; '];
    rejectToleranceStringList{i} = [num2str(rejectTolerance(i)) ';  '];
    checkPeriodStringList{i} = ['[' num2str(checkPeriod(i,1)) ' ' num2str(checkPeriod(i,2)) ']; '];
end
summaryFig = figure('name',[monkeyName expDate protocolName],'numbertitle','off');
fontSize = 8;
h0 = subplot('position',[0.8 0.8 0.18 0.18]); set(h0,'visible','off');
text(0,0.9,['arrays : ' arrayList{:}],'fontsize',fontSize,'unit','normalized','parent',h0);
text(0,0.75,['electrode Numbers :' elecNumString{:}],'fontsize',fontSize,'unit','normalized','parent',h0)
text(0, 0.6, ['thresholds (uV):' thresholdStringList{:}],'fontsize',fontSize,'unit','normalized','parent',h0); %[' num2str(minLimit(1)) ' ' num2str(maxLimit(1)) ']' ' [' num2str(minLimit(2)) ' ' num2str(maxLimit(2)) ']'
text(0, 0.45, ['checkPeriod (s): ' checkPeriodStringList{:}],'fontsize',fontSize,'unit','normalized','parent',h0); 
text(0, 0.3, ['rejectTolerance : ' rejectToleranceStringList{:} ],'fontsize',fontSize,'unit','normalized','parent',h0); %num2str(rejectTolerance)

h1 = getPlotHandles(length(arrayStringList),1,[0.07 0.07 0.7 0.7],0,0);
for i=1:length(arrayStringList)
    startElec = badTrialsData(i).checkTheseElectrodes(1);
    endElec = badTrialsData(i).checkTheseElectrodes(end);
    yticks = startElec:10:endElec;
    subplot(h1(length(arrayStringList)+1-i));
    
imagesc(1:numTrials,endElec:-1: startElec,flipud(badTrialsMatrix{i}(startElec:endElec,:)),'parent',h1((length(arrayStringList)+1-i)));
set(gca,'YDir','normal','Ytick',yticks,'YtickLabels',yticks,'box','off'); colormap(gray);
    if i~=1
        set(gca,'XColor','none')
    end
ylabel(h1(length(arrayStringList)+1-i),arrayStringList{i},'fontsize',10,'fontweight','bold');
end
xlabel(h1(length(arrayStringList)),'# trial num','fontsize',15,'fontweight','bold');
annotation('textbox',[0.045,0.25,0.1,0.1],'String','# electrode num','fontsize',15,'FontWeight','bold','Rotation',90,'LineStyle','none','unit','normalized');


h2 = getPlotHandles(1,1,[0.07 0.8 0.7 0.17]);
h3 = getPlotHandles(length(arrayStringList),1,[0.8 0.07 0.18 0.7],0,0);
subplot(h2); cla; set(h2,'nextplot','add');
stem(h2,1:numTrials,sum(allBadTrialsMatrix,1)); axis('tight');
ylabel('#count');
if ~isempty(badTrials)
    stem(h2,badTrials,sum(allBadTrialsMatrix(:,badTrials),1),'color','r');
end
for i=1:length(arrayStringList)
    startElec = badTrialsData(i).checkTheseElectrodes(1);
    endElec = badTrialsData(i).checkTheseElectrodes(end);
    yticks = startElec:10:endElec;
    subplot(h3(length(arrayStringList)+1-i)); cla; set(h3(i),'nextplot','add');
    stem(h3(length(arrayStringList)+1-i),startElec:endElec,sum(badTrialsMatrix{i}(startElec:endElec,:),2)); axis('tight'); ylabel('#count');
    set(gca,'xtick',yticks,'XTickLabel',yticks);
    if ~isempty(badTrialsData(i).badElecs)
       hold(h3(length(arrayStringList)+1-i),'on')
       stem(h3(length(arrayStringList)+1-i),badTrialsData(i).badElecs,sum(badTrialsMatrix{i}(badTrialsData(i).badElecs,:),2),'color','r'); 
    end
    
    ylims{i} = get(gca,'ylim');
    
    if i~=1
        set(gca,'YColor','none')
    end
    set(gca,'box','off')
    view([90 -90]);
end
ylimsNew = [min(cell2mat(ylims)) max(cell2mat(ylims))];
for i=1:length(arrayStringList)
    set(h3(length(arrayStringList)+1-i),'ylim',ylimsNew)
end
saveas(summaryFig,fullfile(folderSegment,[monkeyName expDate protocolName 'summmaryBadTrials.fig']),'fig');
end
