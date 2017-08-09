
%Displays raw and processed P300 data

function [CCdata, NNdata, Spec, ff, NumSequences, AIc, AIn,  TargWait, GAZE, FaceData] = ...
    ShowMeTheP300(Name,Session,Run,RerefVal,Art,DoCSP,capver)

tic
Figures_On = 0;


currdir = cd;
runloc = currdir(1);

SequencestoRemove = 0;  %For 20 flashes, there are 10 sequences.

%
% if BuildClassifier == 0
%     if Batch == 0
%         CSession = input('What Session? (input as string)   ');
%         CRun = input('What Run? (input as string)   ');
%     elseif Batch == 1;
%         CRun = cell2mat(CRun);
%     end
%     %     NameBase = ['C:\Documents and Settings\amg5106\My Documents\'...
%     %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' Name CSession '\'...
%     %         Name 'S' CSession 'R' CRun '_'];
%
%     NameBase = ['C:\Documents and Settings\amg5106\My Documents\'...
%         'Year4\BCI2000 3.0.5\data\' Name CSession '\'...
%         Name 'S' CSession 'R' CRun '_'];
% end
%

% EOGloc  = 14:16;
% EEGloc = 1:13;
% EEGloc = 1:5;
% EOGloc = [];
% EEGloc = 1:14;
% EOGloc = 15:16;
% EOGloc = [];
EEGloc = 1:8;
EOGloc = [1];

%% Load Data

Data = [];
StimulusCode = [];
StimulusType = [];
InputText = [];
SourceTime = [];
StimulusTime = [];
SelectedStimulus = []; % For audio speller
NewRunName = [];
TTSpell = [];
GazeX = [];
GazeY = [];
SequencePhase = [];
Run


for rr = 1:length(Run)
    %     [a b c d] = load_bcidat(['C:\Documents and Settings\amg5106\My Documents\'...
    %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' ...
    %         Name Session '\' Name 'S' Session 'R' Run{rr} '.dat'],'-calibrated');
    if runloc == 'C'
        try
            [a b c d] = load_bcidat(['data\' ...
                Name Session '\' Name 'S' Session 'R' Run{rr} '.dat'],'-calibrated');
            %disp(['Time to load data ' num2str(toc)]);
        catch
            [a b c d] = load_bcidat(['P:\ALS Proj Data\' ...
                Name Session '\' Name 'S' Session 'R' Run{rr} '.dat'],'-calibrated');
            %disp(['Time to load data ' num2str(toc)]);
        end
    else
        addpath(genpath('Tools'));
        [a b c d] = load_bcidat(['/gpfs/work/a/amg5106/ALS_P300/Data/' ...
            Name Session '/' Name 'S' Session 'R' Run{rr} '.dat'],'-calibrated');
    end
    
    
    full_code = c.SamplingRate.NumericValue*c.StimulusDuration.NumericValue/1000;
    %disp(['full_code = ' num2str(full_code)]);
    
    
   
    
    fs = c.SamplingRate.NumericValue;
    %Trim Data and StimulusCodes
    trimMax = round((c.PreRunDuration.NumericValue-.5)*fs);
    Data = [Data;  a(trimMax:end,:)];
    StimulusCode = [StimulusCode; b.StimulusCode(trimMax:end)];
    StimulusType = [StimulusType; b.StimulusType(trimMax:end)];
    SourceTime = [SourceTime; b.SourceTime(trimMax:end)];
    StimulusTime = [StimulusTime; b.StimulusTime(trimMax:end)];
     if isfield(c,'ToBeCopied') %audio speller
            TTSpell = [TTSpell c.ToBeCopied.NumericValue];
            TTSpell_cell{rr} = c.ToBeCopied.NumericValue;
            SelectedStimulus = [SelectedStimulus; b.SelectedStimulus(trimMax:end)];
        else
            TTSpell = [TTSpell cell2mat(c.TextToSpell.Value)];
            TTSpell_cell{rr} = cell2mat(c.TextToSpell.Value);
     end
    if isfield(b,'EyetrackerLeftEyeGazeX')
        GazeX = [GazeX; b.EyetrackerLeftEyeGazeX(trimMax:end)];
        GazeY = [GazeY; b.EyetrackerLeftEyeGazeY(trimMax:end)];
    end
    SequencePhase = [SequencePhase; b.PhaseInSequence(trimMax:end)];
    
    NewRunName = strcat(NewRunName,Run{rr});
    if isfield(c,'ToBeCopied') %audio speller
        Speller = 4;
    else
    if runloc == 'C'
        fid = fopen(['data\' Name Session '\'...
            Name 'S' Session 'R' num2str(str2num(Run{rr})) '+_mylogfile.txt'],'r');
        if fid==-1
            fid = fopen(['P:\ALS Proj Data\' Name Session '\'...
                Name 'S' Session 'R' num2str(str2num(Run{rr})) '+_mylogfile.txt'],'r');
        end
    else
        fid = fopen(['/gpfs/work/a/amg5106/ALS_P300/Data/' Name Session '/'...
            Name 'S' Session 'R' num2str(str2num(Run{rr})) '+_mylogfile.txt'],'r');
    end
    itt=textscan(fid,'%s','delimiter','\n');
    InputText = [InputText; itt{1}];
    Speller = 1; %multiflash grid
    fclose(fid);
    end
    
    %Associate Target Codes with Targets
    if Speller == 4
        tdtemp = c.Stimuli.Value(1,:);
        TargSymb(rr,:) = tdtemp
    else
        tdtemp = c.TargetDefinitions.Value;
        if size(tdtemp,2)==1 %Multi Menu
            TargSymb(rr,:) = tdtemp{1}(:,2)'; %-- NOTT CORRECT YET
        else
            TargSymb(rr,:) = tdtemp(:,2)';
        end
        SymbSymb = {'Angry%20','Bag%20','Bed%20','Bored%20','Bathroom%20','Carer%20','Cdrink%20',...
            'Clothes%20','Cold%20','Doctor%20','Family%20','Food%20','Glasses%20','Hdrink%20',...
            'Hear%20','Help%20','Hot%20','Spouse%20','Idk%20','Light%20','Medicine%20','Yes%20',...
            'No%20','Nurse%20','Pain%20','Pill%20','Slippers%20','Telephone%20','Tired%20',...
            'Toilet%20','Tv%20','Wheelchair%20','Walker%20','Worried%20','Watch '};
    end
    
    %Number of channels
    NumChans{rr} = size(a,2);
    
    %Number of sequences per trial.  Two times this number is the number of
    %times each stimulus is flashed each trial.
    NumSequences{rr} = c.NumberOfSequences.NumericValue;
    
    %Number of StimulusCodes
    NumStimCodes = max(StimulusCode);
    
    %Duration of the Stimulus State -- This is how long the stimulus is on
    StateDuration = c.StimulusDuration.NumericValue*fs/1000;
    
    %Number of trials
    NumTrials{rr} = fix((sum(b.StimulusCode>=1)/StateDuration)/(NumSequences{rr}*NumStimCodes));
    
    %Load Spatial Filter
    SpatFiltUsed{rr} = c.SpatialFilter.NumericValue;
    
    %Load Classifier
    ClassifierUsed{rr} = c.Classifier.NumericValue;
    ClassifierUsed{rr}(:,4) = cellfun(@str2num,c.Classifier.Value(:,4));
    
end
Data = Data';
range = (round(.45*fs):round(1.2*fs)); %Use this data range for classifiation
% range = (round(.2*fs):round(.6*fs)); %Use this data for checking on the
% computer -- the previous one takes too much memory :(

%  range = (1:round(fs));



%For the new g.Nautilus Caps (with flexible cables), the channel order is
% Cz, Fz, P3, Pz, P4, PO7, PO8, Oz
%The older caps (with the ribbon cables) were
% Fz, Cz, P3, Pz, P4, PO7, PO8, Oz

if capver == 2
    Data = [Data(2,:); Data([1 3:8],:)];
end
% figure
% xx(1)=subplot(121); plot(Data(1:4,:)');
% xx(2)=subplot(122); plot(Datan(1:4,:)');
% linkaxes(xx)



%disp(['Data prep ' num2str(toc)]);
%% Check Variables
%Check if free spelling
FreeSpell = sum(StimulusType)==0&sum(StimulusCode)~=0;


%Are channels consistant across runs?
if sum(diff(cell2mat(NumChans))) ~= 0
    disp('Different number of channels in each run!');
end
NumChans = NumChans{1};

% %Do the channel labels match the number of channels in Data?
% if length(EOGloc)+length(EEGloc)~=size(Data,1)
%     disp('Channel labels incorrect'); return
% end

%Are the spatial filters the same same size?
try %Are classifiers the same size
    s2mC = cell2mat(SpatFiltUsed);
    
    if sum(sum(diff(reshape(s2mC,size(s2mC,1),size(s2mC,2)/rr,rr),[],3))) ~= 0
        disp('Spatial Filters have different values in each run!');
    end
    SpatFiltUsed = SpatFiltUsed{1};
catch err
    disp('Spatial Filters  are different sizes in each run!');
end
%If the same size, do they contain the same values?


try %Are classifiers the same size
    c2mC = cell2mat(ClassifierUsed);
    %If the same size, do they contain the same values?
    if sum(sum(diff(reshape(c2mC,size(c2mC,1),size(c2mC,2)/rr,rr),[],3))) ~= 0
        disp('Classifiers have different values in each run!');
    end
    ClassifierUsed = ClassifierUsed{1};
catch err
    disp('Classifiers are different sizes in each run!');
end

%Are number of sequences the same
if sum(cell2mat(NumSequences))/rr ~= NumSequences{1}
    disp('Not the same number of sequences in each run, code will not run correctly');
end
NumSequences = NumSequences{1};

%Number of rows, columns, stimuli
if Speller==4
    NR = 1; NC = 1;
else
NR = c.NumMatrixRows.NumericValue;
NC = c.NumMatrixColumns.NumericValue;
end

if size(tdtemp,2)==1 %Multi Menu
    NS = size(tdtemp{1},1);
else
    NS = size(tdtemp,1);
end
%Check state duration again
for scc = 1:NumStimCodes
    sdd = find(StimulusCode == scc);
    kk(scc)=2;
    while((sdd(kk(scc))-sdd(kk(scc)-1))==1)
        kk(scc) = kk(scc)+1;
    end
end
if mean(kk) == kk(1)
    StateDuration2 = mean(kk)-1;
else
    disp('We have a problem'); return
end

if StateDuration ~= StateDuration2
    disp('We have a problem'); return
end

%If using the CNEamp - remove the first 2 channels
if isfield(c,'ComPort')==1
    Data = Data(3:end,:);
    cneamp=1;
end

%disp(['Check Vars ' num2str(toc)]);
%% Check Timing of Data
%Check for Dropped Samples (in the case of BCI2000, the data is repeated in
%the next block if the processing exceeds the roundtrip time)
Timerr = diff(SourceTime(1:c.SampleBlockSize.NumericValue:end));
StTimerr = diff(StimulusTime(1:c.SampleBlockSize.NumericValue:end));

BadTime = Timerr<31 | Timerr>32;
BadTimeF = repmat(double(BadTime),1,c.SampleBlockSize.NumericValue)';
BadTimeF = BadTimeF(:);
% if Figures_On == 1
%     figure
%     plot(StTimerr); hold on;
%     plot(Timerr,'r'); ylim([0 250]);
%     figure
%     plot(Data(1,:),'r'); hold on; plot(BadTimeF*20); plot(StimulusCode,'g');
% end
%figure
%plot(Data(:,1)*1000); hold on; plot(b.SourceTime,'r');


TargWait = double(diff(StimulusType)==1); %Locations where the stimulus type changes
% TargWait(StimulusType==-1) = 0; %Remove detected regions which have been removed
TargWaitVal = diff(find(TargWait))'; %This is the time that
%preceded the current target stimulus (starts with Cloc(2))
TargWaitVal = [1 TargWaitVal];
TargWait(TargWait==1) = TargWaitVal;
%Shift forward one so that the Timing of the TargWait and StimulusType/Code
%align
TargWait = circshift(TargWait,[1 0]);

%% Artifact Correction

if Art == 1 %Remove data displaying artifacts
    try
        ArtInt = (Data(EOGloc(1),:)-Data(EOGloc(2),:)>75 | ...
            Data(EOGloc(1),:)-Data(EOGloc(2),:)<-75 | ...
            Data(EOGloc(1),:)-Data(EOGloc(3),:)>75 | ...
            Data(EOGloc(1),:)-Data(EOGloc(3),:)<-75 | ...
             Data(EOGloc(2),:)-Data(EOGloc(3),:)>75 | ...
            Data(EOGloc(2),:)-Data(EOGloc(3),:)<-75);
        AI = ArtInt;
        ArtInt = find(ArtInt);
    catch
        disp('Cannot remove data because no EOG channels were recorded');
        return
    end
    
    GoodData = ones(size(Data,2),1);
    
    %Changed from 2 seconds after 5/29/17
    %Define a range 1 seconds before and .5 seconds after all of the found
    %artifacts.  This big of a range is used because we want to ensure that
    %we dont keep stimulus codes that could be affected in an upcoming
    %eye blink, and that we dont have the .5 seconds following
    %an eye blink as the beginning of a stimulus code.
    arange = fix(-1*fs:.5*fs);
    
    BadData = repmat(ArtInt,length(arange),1)+repmat((arange)',1,length(ArtInt));
    BadData = unique(BadData);
    GoodData = ~ismember(1:length(Data),BadData);
    
    if Figures_On == 1
        figure('Position',[50 50 1200 700])
        nn(1)=subplot(3,1,[1 2])
        plot(Data(EOGloc(1),:)-Data(EOGloc(2),:)); hold on;
        plot(Data(EOGloc(2),:)-Data(EOGloc(3),:),'g');
        plot(Data(EEGloc(10),:),'k');
        plot(100*(1-GoodData),'r');
    end
    
    
    GoodData = logical(GoodData);
    
    
    %Not trimming bad data -- too complicated to change all the remaining code
    %to work with non-regular stimulus codes (much is based on having exactly
    %16 codes of 10 repetitions each. I TAKE IT BACK, I AM :)
    %
    
    %Set these codes to .1, not 0 because i still need to register them later
    %to construct StimCode.  I just wont access them for P300 building.
    StimulusCode = double(StimulusCode);
    StimulusType = double(StimulusType);
    %Find inidices within the bad data segment for which StimulusCode is up.
    %Set bad inds to 0 right now, will change to .1 in a few lines.
    Bad_ind = find(StimulusCode>0&~GoodData');
    StimulusCode(Bad_ind) = 0;
    StimulusType(Bad_ind) = 0;
    
    %Get rid of partial stim codes
    UD_sc = diff(double(StimulusCode));
    U_sc = find(UD_sc>0);
    D_sc = find(UD_sc<0);
    
    if Figures_On == 1
    nn(2) = subplot(3,1,3);
        plot(StimulusCode);
        hold on
        plot(GoodData,'r');
    end
    
    not_full = find((D_sc-U_sc)~=full_code);
    for nf = 1:length(not_full)
        Bad_ind = [Bad_ind; (U_sc(not_full(nf))+1:D_sc(not_full(nf)))'];
    end
    
    StimulusCode(Bad_ind) = -1;
    StimulusType(Bad_ind) = -1;
    TargWait(Bad_ind) = 0;
    
    if Figures_On == 1
    plot(StimulusCode,'g');
    plot(TargWait./100,'k');
    linkaxes(nn,'x');
    end
    
    
    ArtWeights = eye(size(Data,1));
    ArtWeights2 = eye(length(EEGloc),size(Data,1));
    NumChans = length(EEGloc);
    
elseif Art == 2 || Art == 4 %Artifact regression
    %     if BuildClassifier == 1
    %Build Artifact rejector (output chans x input chans)
    %Find the eog & eeg channels
    
    %     EOGc = regexp(c.ChannelNames.Value,'EOG');
    %     Cloc = cellfun(@isempty,EOGc);
    %     EEGloc = find(Cloc~=0);
    %     EOGloc = find(Cloc==0);
    %
    %         %use the regression algorithm
    %     EOGD = Data(1:5*fs,:);
    
    
    %         figure
    %         plot(Data(:,EOGloc(1)));
    %         EOGDloc = input('Data points to use for EOG regression: ');
    try
        if length(EOGloc)==3
            ArtInt = (Data(EOGloc(1),:)-Data(EOGloc(2),:)>75 | ...
                Data(EOGloc(1),:)-Data(EOGloc(2),:)<-75 | ...
                Data(EOGloc(1),:)-Data(EOGloc(3),:)>75 | ...
                Data(EOGloc(1),:)-Data(EOGloc(3),:)<-75 | ...
                Data(EOGloc(2),:)-Data(EOGloc(3),:)>75 | ...
                Data(EOGloc(2),:)-Data(EOGloc(3),:)<-75);
            AI = ArtInt;
            ArtInt = find(ArtInt);
        else
            ArtInt = (Data(EOGloc(1),:)>75 | ...
                Data(EOGloc(1),:)<-75);
            AI = ArtInt;
            ArtInt = find(ArtInt);
            
        end
    catch
        disp('Cannot remove data because no EOG channels were recorded');
        return
    end
    
    %Here, we use this range because we are only interested in capturing
    %the time of the eye blink
    Arange = fix(-.1*fs:.5*fs);
    
    ArtInt = repmat(ArtInt',1,length(Arange))+repmat(Arange,length(ArtInt),1);
    ArtInt = ArtInt';
    ArtInt = unique(ArtInt(:));
    ArtInt = ArtInt(ArtInt>0 & ArtInt<size(Data,2));
    ArtLoc = logical(zeros(size(Data,2),1));
    ArtLoc(ArtInt) = true;
    
    %         if Figures_On == 1
    %             if length(EOGloc)==3
    %                 figure
    %                 plot(Data(EOGloc(1),:)-Data(EOGloc(2),:))
    %                 hold on
    %                 plot(Data(EOGloc(1),:)-Data(EOGloc(3),:))
    %                 plot(100*ArtLoc,'r')
    %             else
    %                 figure
    %                 plot(Data(EOGloc(1),:)); hold on;
    %                 plot(100*ArtLoc,'r');
    %             end
    %         end
    
    EOGD = Data(:,ArtLoc);
    
    if ~isempty(EOGD)
        %the function covm adds an additional column of ones in front of the data
        %and is necessary for regress_eog.m
        if length(EOGloc)==3
            [R] = regress_eog(covm(EOGD','E'),EEGloc, ...
                sparse([EOGloc(1),EOGloc(3),EOGloc(2),EOGloc(1)],[1,1,2,2],[1,-1,1,-1]));
        elseif length(EOGloc)==2
            [R] = regress_eog(covm(EOGD','E'),EEGloc, ...
                sparse([EOGloc(1),EOGloc(2)],[1,1],[1,-1]));
        else
            [R] = regress_eog(covm(EOGD','E'),EEGloc, EOGloc);
        end
        %Create full matrix for online artifact reduction
        %I believe this is the way they say to do it (pad Data with a channel
        %of ones -- this introduces a bias to the output channel) (see DD2 below).
        %However, this padding is not something I want to do online, and since
        %it is only a bias, we can remove the first column of ArtWeights.
        ArtWeights = full(R.r0)';
        ArtWeights2 = ArtWeights(EEGloc,:);
        NumChans = length(EEGloc);
        %         ArtWeights = full(R.r1)';
        %         ArtWeights2 = ArtWeights(:,2:end);
        %DD2 = [ones(size(Data,1),1),Data] * ArtWeights';
    else
        ArtWeights2 = eye(length(EEGloc),size(Data,1));
        NumChans = length(EEGloc);
    end
    
    %     elseif BuildClassifier == 0
    %     end
    
else
    ArtWeights = eye(size(Data,1));
    ArtWeights2 = eye(length(EEGloc),size(Data,1));
    NumChans = length(EEGloc);
    AI = zeros(1,size(Data,2));
end


%Using the correction coefficients, transform entire training run to reduce
%artifact
Dw = ArtWeights2*Data;
% if Figures_On == 1
%     if BuildClassifier == 1
%Plot the F and O raw and artifacted, as well as EOG channels
if Figures_On == 1
    
    figure
    tt(1)=subplot(311); plot(Data(1,:)); hold on; plot(Dw(1,:),'r');
    tt(2)=subplot(312); plot(Data(10,:)); hold on;
    plot(Dw(10,:),'r');
    tt(3)=subplot(313); plot(Data(EOGloc,:)'); hold on
    linkaxes(tt,'x');
    
end
%         clear Dw;
%     end
% end
%
% %Plot the spectra of the raw vs rerefed and artifacted channels
%     figure
%     kk=1;
%     for ich =[EEGloc(1) EEGloc(2) EEGloc(end-1-length(EOGloc)) EEGloc(end-length(EOGloc))]
%         subplot(2,2,kk);
%         [Sp1,ff1]=pwelch(Data(:,ich),2*fs,fs,2*fs,fs); hold on;
%         plot(ff1,10*log10(Sp1));
%         [Sp2,ff2]=pwelch(DD(:,ich),2*fs,fs,2*fs,fs);
%         plot(ff2,10*log10(Sp2),'r');
%         kk=kk+1;
%     end
%
% DD = DD(:,EEGloc);
% NumChans = length(EEGloc);


%% Find remaining (>70 mV) artifacts in central
%  channels (4,5,6,9,10,11,14,15,16,18,19)
if Art == 4
    %If not using the contents of this cell, set GoodData to ones;
    GoodData2 = ones(size(Dw,2),1);
    
    
    BadData2 = [];
    for i = [4 5 6 9 10 11 14 15 16 18 19]
        %Find large artifacts
        a_thresh = std(Dw(i,:));
        BadData2 = [BadData2 find(Dw(i,:)>5*a_thresh|Dw(i,:)<-5*a_thresh)];
    end
    %Define a range 1 second before and 2 seconds after all of the found
    %artifacts
    BadData2 = repmat(BadData2,3*fs,1)+repmat((-1*fs+1:2*fs)',1,length(BadData2));
    BadData2 = unique(BadData2);
    GoodData2 = ~ismember(1:length(Data),BadData2);
    
    if Figures_On == 1
        figure
        for i = 1:size(Dw,1)
            if size(Dw,1) == 19
                eleclocs = [7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 29];
                tr(i) = subplot(6,5,eleclocs(i)); plot(Dw(i,:));
                hold on;
                plot(100*(1-GoodData2),'r');
            else
                tr(i) = subplot(6,5,i); plot(Dw(i,:));
                hold on;
                plot(100*(1-GoodData2),'r');
            end
        end
        linkaxes(tr);
    end
    
    
    GoodData2 = logical(GoodData2);
    
    
    %Not trimming bad data -- too complicated to change all the remaining code
    %to work with non-regular stimulus codes (much is based on having exactly
    %16 codes of 10 repetitions each. I TAKE IT BACK, I AM :)
    %
    
    %Set these codes to .1, not 0 because i still need to register them later
    %to construct StimCode.  I just wont access them for P300 building.
    StimulusCode = double(StimulusCode);
    StimulusType = double(StimulusType);
    %Find inidices within the bad data segment for which StimulusCode is up.
    %Set bad inds to 0 right now, will change to .1 in a few lines.
    Bad_ind2 = find(StimulusCode>0&~GoodData2');
    StimulusCode(Bad_ind2) = 0;
    StimulusType(Bad_ind2) = 0;
    
    %Get rid of partial stim codes
    UD_sc = diff(double(StimulusCode));
    U_sc = find(UD_sc>0);
    D_sc = find(UD_sc<0);
    
    if Figures_On == 1
        figure
        plot(StimulusCode);
        hold on
        plot(GoodData2,'r');
    end
    
    not_full = find((D_sc-U_sc)~=full_code);
    for nf = 1:length(not_full)
        Bad_ind2 = [Bad_ind2; (U_sc(not_full(nf))+1:D_sc(not_full(nf)))'];
    end
    
    StimulusCode(Bad_ind2) = -1;
    StimulusType(Bad_ind2) = -1;
    TargWait(Bad_ind2) = 0;
    
    if Figures_On == 1
    plot(StimulusCode,'g');
    plot(TargWait./100,'k');
    end
    
    
    
end

TargWait = TargWait(TargWait~=0)';
%disp(['Artifacting ' num2str(toc)]);

%% Rereference Data
switch RerefVal
    case 0
        RefFilt = eye(length(EEGloc));
        %         DataF = Data*SpatFilt;
    case 1
        RefFilt = (-1/length(EEGloc))*(ones(length(EEGloc))-eye(length(EEGloc)));
        RefFilt = RefFilt + eye(length(EEGloc));
        %         SpatFilt = blkdiag(SpatFilt, eye(length(EOGloc)));
        %         DataF = Data*SpatFilt;
    case 2
        if size(Data,1)==16
            RefFilt = [-1/4 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0 0 0;
                0 -1/4 0 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0;
                0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 0 -1/4 0;
                0 0 0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 -1/4];
        elseif size(Data,1)==22
            %This is for the 22 channel setup
            RefFilt = [-1/4 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0 0 0 0 0 0 0 0;
                0 -1/4 0 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0 0 0 0 0 0;
                0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0 0 0;
                0 0 0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0;
                0 0 0 0 0 0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 0 -1/4 0;
                0 0 0 0 0 0 0 0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 -1/4];
        end
        NumChans_old = NumChans;
        NumChans = size(RefFilt,1);
    case 3
        if size(Data,1)==16
            RefFilt = [0 0 0 0 -1 0 0 0 0 0 0 0 1 0;
                0 0 0 0 -1 0 0 0 0 0 0 0 0 1];
        elseif size(Data,1)==22
            RefFilt = [0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 1 0;
                0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1];
        end
        NumChans_old = NumChans;
        NumChans = 2;
    case 4 %reference to C3
        if size(Data,1)==16
        elseif size(Data,1)==22
            RefFilt = eye(length(EEGloc));
            tmp = find(strcmp(c.ChannelNames.Value,'C3'));
            RefFilt(:,tmp) = RefFilt(:,tmp)-1;
        end

end
%disp(['Reref ' num2str(toc)]);


%% Do CSP
if DoCSP == 1
    %     if BuildClassifier == 1
    disp('CSP not configured to work with data that has been rejected!')
    Data_C = (RefFilt*ArtWeights2*Data)';
    [W] = CSP_P300(Data_C,StimulusCode,StimulusType,...
        NumChans, NumTrials, NumStimCodes,range,StateDuration);
    %         [W2] = CSP_P300_ver2(Data_C,StimulusCode,StimulusType,...
    %             NumChans, NumTrials, NumStimCodes,range,StateDuration);
    
    %     elseif BuildClassifier == 0
    %     end
elseif DoCSP == 0
    %     if BuildClassifier == 1
    W = eye(size(RefFilt,1));
    W2 = eye(size(RefFilt,1));
    %     else
    %     end
end

%disp(['CSP ' num2str(toc)]);


%     ChansToKeep = input('Which Channels to keep?'
if DoCSP == 1
    if RerefVal == 1 %This is because the
        %last CSP channel unser CAR rereferencing is zeros.
        if NumChans <= 4
            ChansToKeep = 1:NumChans-1;
        else
            ChansToKeep = [1 2 NumChans-2 NumChans-1];
        end
    else
        if NumChans <= 3
            ChansToKeep = 1:NumChans;
        else
            ChansToKeep = [1 2 NumChans-1 NumChans];
        end
    end
    
else
    ChansToKeep = 1:size(W,1);
end

NumChans = length(ChansToKeep);
W = W(ChansToKeep,:);





SpatFilt = W*RefFilt*ArtWeights2;

DD = SpatFilt*Data;
if isempty(EOGloc);
    EOGDD = NaN(0,size(Data,2));
     EOGDD2 = NaN(0,size(Data,2));
else
EOGDD = Data(EOGloc,:);
EOGDD2 = [Data(EOGloc(2),:)-Data(EOGloc(1),:);...
    Data(EOGloc(2),:)-Data(EOGloc(3),:);
    Data(EOGloc(1),:)-Data(EOGloc(3),:)];
end



%% Compute Spectra
for cc = 1:size(DD,1)
    [Spec(:,cc), ff] = pwelch(DD(cc,:),fs*8,0,fs*8,fs);
    
    
    %Same as Chronux Multitaper Below
    % [E,V] = dpss(length(DD(cc,:)),2);
    % h = spectrum.mtm(E,V);    % Specify DPSS and concentrations
    % % when creating the MTM spectral estimator.
    % S3 = psd(h,DD(cc,:),'Fs',fs);
    % Spec3(cc,:) = S3.Data; ff3 = S3.Frequencies;
    
    
end

% params.tapers = [3 5];
% params.Fs = fs;
% [Spec2,ff2] = mtspectrumc(DD',params);



% figure
% tt(1) = subplot(211); plot(ff,Spec);
% linkaxes(tt,'x');



%disp(['Compute spectra ' num2str(toc)]);


if FreeSpell == 1
%     Clocs = find(StimulusCode>0);
%     Clocs = Clocs(1:StateDuration:end);
%     Nlocs = find(StimulusCode>0);
%     Nlocs = Nlocs(1:StateDuration:end);
    CCdata = [];%zeros(length(range),size(DD,1),length(Clocs));
    NNdata = [];%zeros(length(range),size(DD,1),length(Nlocs));
    AIc = []; AIn = []; GAZE = []; FaceData = [];
    return
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find code stimuli associations (contained in SC and SR), as well as the
%stimuli sequences for each trial.



if Speller == 1  %Use the myfile.txt associated with CB speller runs to
    %find stimulus codes and target associations
    
    for i = 1:length(InputText)
        RunEnd(i) = ~isempty(strfind(InputText{i},'******************'));
    end
    Runind = find(RunEnd==1);
    SCindF = []; SRindF = []; INITmindF = []; NEXTmindF = []; FaceOF = [];
    for rr = 1:length(Runind);
        if rr == 1;
            RunRange = 1:Runind(1);
        else
            RunRange = Runind(rr-1):Runind(rr);
        end
        clear Begin SCind SRind INITmind NEXTmind Faceind
        for ll = RunRange
            Begin(ll) = ~isempty(strfind(InputText{ll},'*****START'));
            SCind(ll) = ~isempty(strfind(InputText{ll},'SC'));
            SRind(ll) = ~isempty(strfind(InputText{ll},'SR'));
            INITmind(ll) = ~isempty(strfind(InputText{ll},'INIT_mSequence'));
            NEXTmind(ll) = ~isempty(strfind(InputText{ll},'NEXT_mSequence'));
            Faceind(ll) = ~isempty(strfind(InputText{ll},'Face Overlays'));
        end
        Begind = find(Begin==1);
        SCind = find(SCind==1)'; FSCi = find(SCind < Begind);
        SCindF = [SCindF; SCind(FSCi(end):FSCi(end)+NumTrials{rr}-1)];
        SRind = find(SRind==1)'; FSRi = find(SRind < Begind);
        SRindF = [SRindF; SRind(FSRi(end):FSRi(end)+NumTrials{rr}-1)];
        INITmind = find(INITmind==1)';
        INITmindF = [INITmindF; INITmind(1:NumTrials{rr})];
        NEXTmind = find(NEXTmind==1)';  omitN = NumSequences:...
            NumSequences:NumTrials{rr}*NumSequences;
        inrange = omitN<length(NEXTmind);
        omitN = omitN(inrange);
        NEXTmind(omitN) = [];
        NEXTmindF = [NEXTmindF; NEXTmind];
        FaceO = find(Faceind==1);
        FaceO = InputText(FaceO+1);
        FaceOF = [FaceOF; FaceO];
        
    end
    
    
    for i = 1:sum(cell2mat(NumTrials))
        tempSC = InputText{SCindF(i)}(4:end);
        tempSC = regexp(tempSC,' ','split')'; tempSC(cellfun(@isempty,tempSC)) = [];
        SC{i} = tempSC;
        tempSR = InputText{SRindF(i)}(4:end);
        tempSR = regexp(tempSR,' ','split')'; tempSR(cellfun(@isempty,tempSR)) = [];
        SR{i} = tempSR;
        tempINITm = InputText{INITmindF(i)}(16:end);
        tempINITm = regexp(tempINITm,' ','split')'; tempINITm(cellfun(@isempty,tempINITm)) = [];
        INITm{i} = tempINITm;
        for j = 1:NumSequences-1
            tempNEXTm = InputText{NEXTmindF((NumSequences-1)*(i-1)+j)}(16:end);
            tempNEXTm = regexp(tempNEXTm,' ','split')'; tempNEXTm(cellfun(@isempty,tempNEXTm)) = [];
            NEXTm{j,i} = tempNEXTm;
        end
    end
    mSeq = cat(1,INITm,NEXTm);
    %%%%%%%%%%%
    
elseif Speller == 2 %RSVP
    for i = 1:sum(cell2mat(NumTrials))
        SC{i} = 1:32;
        SR{i} = 33:64;
    end
    
elseif Speller == 3
    for i = 1:sum(cell2mat(NumTrials))
        SC{i} = repmat(5:12,1,NR);
        SR{i} = repmat(1:4,NC,1);
        SR{i} = SR{i}(:)';
    end
    
elseif Speller == 4
    for i = 1:sum(cell2mat(NumTrials))
        SC{i} = repmat(7:12,1,NR);
        SR{i} = repmat(1:6,NC,1);
        SR{i} = SR{i}(:)';
    end
    GoalLetterF = c.Stimuli.Value(1,TTSpell(:));
    GoalTrialF = 1:sum(cell2mat(NumTrials));
    ChosenLetterF = SelectedStimulus(SelectedStimulus>0);
    ChosenLetterF = ChosenLetterF(1:fs*c.PostSequenceDuration.NumericValue:end);
end

NumTrials = sum(cell2mat(NumTrials));

%Check that the markers saved in the log file are the same as those in the
%data file
tStimC = double(nonzeros(StimulusCode)); tStimC = tStimC(1:StateDuration:end);
for i = 1:NumTrials
    for j = 1:NumSequences
        StimCode{j,i} = num2cell(tStimC(NumStimCodes*NumSequences*(i-1)+NumStimCodes*(j-1)+1:...
            NumStimCodes*NumSequences*(i-1)+NumStimCodes*(j)));
    end
end

if Speller == 1
    %StimCode should contain the same values as mSeq
    ies=0;
    for i = 1:NumTrials
        for j = 1:NumSequences
            ies = isequal(cell2mat(StimCode{j,i}),str2double(mSeq{j,i}))+ies;
        end
    end
    if ies ~= NumTrials*NumSequences
        disp('Sequences inconsistant (prob because of data rejection');
        disp([num2str((NumTrials*NumSequences)-ies) ' sequences mismatched']);
    end
end


%If removing data(probably wont do) need to also remove the same portion of
%StimCode





NumSequences = NumSequences-SequencestoRemove;
StimCode = StimCode(1:NumSequences,:);


%disp(['Acess stim codes ' num2str(toc)]);


%% EyeGaze

if Speller==4 %Can include this later -- for now just ignore eye tracking for audio speller
    GAZE = [];
else
    if ~isempty(GazeX) && Art==0
        if runloc == 'C'
        end
        [GAZE.acc, GAZE.var, GAZE.invalid, GAZE.targetbox] = EyeTracking(...
            GazeX,GazeY, Data(EOGloc,:), SC,SR, StimulusCode, ...
            SequencePhase, StimulusType, StateDuration, ...
            NumSequences+SequencestoRemove, Speller,...
            fs, c, 0, 1);
        cd(currdir);
    else
        GAZE = [];
    end
end
    %disp(['Finished with eye tracking ' num2str(toc)]);

%% Organize Transformed Data into Choice/Not Choice

if FreeSpell == 0
    Clocs = find(StimulusCode>0 & StimulusType == 1);
    Clocs = Clocs(1:StateDuration:end);
    Nlocs = find(StimulusCode>0 & StimulusType == 0);
    Nlocs = Nlocs(1:StateDuration:end);
    
    %Check if Nlocs or Clocs will go out of range (run finished early)
    outofit = Nlocs+max(range)<size(DD,2);
    Nlocs = Nlocs(outofit);
    outofit = Clocs+max(range)<size(DD,2);
    Clocs = Clocs(outofit);
    
    
    CCdata = zeros(length(range),size(DD,1)+size(EOGDD,1)+size(EOGDD2,1),length(Clocs));
    NNdata = zeros(length(range),size(DD,1)+size(EOGDD,1)+size(EOGDD2,1),length(Nlocs));
    AIc = zeros(length(range),length(Clocs));
    AIn = zeros(length(range),length(Nlocs));
    
    
    for cl = 1:length(Clocs)
        CCdata(:,:,cl) = cat(2,DD(:,Clocs(cl)+range)', EOGDD(:,Clocs(cl)+range)', EOGDD2(:,Clocs(cl)+range)');
        AIc(:,cl) = AI(Clocs(cl)+range);
    end
    for nl = 1:length(Nlocs)
        NNdata(:,:,nl) = cat(2,DD(:,Nlocs(nl)+range)', EOGDD(:,Nlocs(nl)+range)', EOGDD2(:,Nlocs(nl)+range)');
        AIn(:,nl) = AI(Nlocs(nl)+range);
    end
    
    
    trlrem = squeeze(sum(sum(CCdata>100 | CCdata<-100)))>0;
    CCdata(:,:,trlrem) = [];
    trlrem = squeeze(sum(sum(NNdata>100 | NNdata<-100)))>0;
    NNdata(:,:,trlrem) = [];
    %figure(10) %Plot averages
    % for ch = 1:NumChans
    % subplot(5,5,elocs(ch)); plot(range/fs,mean(CCdata(:,ch,:),3))
    %     hold on;
    %     plot(range/fs,mean(NNdata(:,ch,:),3),'r')
    % end
    
    % figure %Plot averages
    % for ch = 1:NumChans
    %     subplot(4,5,ch); errorbar(mean(CCdataA(:,ch,:),3),std(CCdataA(:,ch,:),[],3)/sqrt(size(CCdataA,3)))
    %     hold on;
    %     errorbar(mean(NNdataA(:,ch,:),3),std(NNdataA(:,ch,:),[],3)/sqrt(size(NNdataA,3)),'r')
    % end
    
    % figure %Plot averages
    % for ch = 1:NumChans
    %     subplot(4,5,ch); errorbar(mean(CCdataAF(:,ch,:),3),std(CCdataAF(:,ch,:),[],3)/sqrt(size(CCdataAF,3)))
    %     hold on;
    %     errorbar(mean(NNdataAF(:,ch,:),3),std(NNdataAF(:,ch,:),[],3)/sqrt(size(NNdataAF,3)),'r')
    % end
    %
    % figure %Plot averages
    % for ch = 1:NumChans
    %     subplot(4,5,ch); errorbar(mean(CCdataAFC(:,ch,:),3),std(CCdataAFC(:,ch,:),[],3)/sqrt(size(CCdataAFC,3)))
    %     hold on;
    %     errorbar(mean(NNdataAFC(:,ch,:),3),std(NNdataAFC(:,ch,:),[],3)/sqrt(size(NNdataAFC,3)),'r')
    % end
    %
    % figure %Plot averages
    % for ch = 1:NumChans
    %     subplot(4,5,ch); errorbar(mean(CCdataAFCT(:,ch,:),3),std(CCdataAFCT(:,ch,:),[],3)/sqrt(size(CCdataAFCT,3)))
    %     hold on;
    %     errorbar(mean(NNdataAFCT(:,ch,:),3),std(NNdataAFCT(:,ch,:),[],3)/sqrt(size(NNdataAFCT,3)),'r')
    % end
    
    % %For NIH Grant-used no reref, no art, yes csp -- plotted by running through
    % %each subject, plotting, saving, and adding to the previous plot.
    % %saved as NIHgrantP300_st_an_sm_aj.fig
    % %S1 - AndrewP3_14001 R03/04 -- csp chan 1
    % %S2 - SumithraP3001 R03/04 -- csp chan 1
    % %S3 - AnjumGridCB001 R01 -- csp chan 1
    % %S4 - SteveGridCB002 R03 -- csp chan 1
    % ch=1;
    % figure
    % hold on
    % plot(range/fs,mean(CCdataAFC(:,ch,:),3)); hold on;
    % plot(range/fs,mean(NNdataAFC(:,ch,:),3),'r'); xlim([0 range(end)/fs])
    % xlabel('Time (s)'); ylabel('Amplitude \muV');
    
    
    
    
    % end
else %FreeSpell = 1 - should have exited already
    
end




% %For combining data as the output
% CC.chans = {'Fp1','C3','C4','Pz','O1','lEOG','cEOG','rEOG'};
% CC.data = Data([ 1 9 11 15 18 20:22],:);
% if Art == 2 || Art == 4
%     CC.data_cor = Dw([1 9 11 15 18],:);
% end
% CC.blink_ind = AI;
% CC.stimcode = StimulusCode;
% CC.stimtype = StimulusType;
% if Art == 4
%     CC.art_ind = GoodData2;
% end







%StimulusType are repeated
%
% [aa bb] = butter(5,[.5 30]/128);
% aaa = filtfilt(aa,bb,Data);
%
% figure
% plot(Data(:,1)); hold on; plot(aaa(:,1),'r');
%
% figure
% pwelch(Data(:,1),1024,0,512,256); hold on;
% pwelch(aaa(:,1),1024,0,512,256);
%
% figure
% plot(StimulusCode);
% hold on;
% plot(StimulusType,'r');


%     %Reshape data so it can be multiplied by Wp
%     Data2=reshape(AvgEEGData,NumStimCodes*NumTrials,NumChans,size(AvgEEGData,2));
%     Z = zeros(size(Data2,1),size(Wp,1),size(Data2,3));
%     %Transform data with spatial filter Wp
%     for tt = 1:size(Data2,1)
%         Z(tt,:,:) = Wp*squeeze(Data2(tt,:,:));
%     end
%     %Reshape data
%     Z2 = reshape(Z,NumStimCodes*NumTrials*NumChans,size(Data2,3));


% P300 sanity check (2)





% %% Save/Load Spatial Filter
% if BuildClassifier == 1
%
%     %Which CSP channels to Keep?
%
%
%     %     ChansToKeep = input('Which Channels to keep?'
%     if DoCSP == 1
%         if RerefVal == 1 %This is because the
%             %last CSP channel unser CAR rereferencing is zeros.
%             if NumChans <= 4
%                 ChansToKeep = 1:NumChans-1;
%             else
%                 ChansToKeep = [1 2 NumChans-2 NumChans-1];
%             end
%         else
%             if NumChans <= 3
%                 ChansToKeep = 1:NumChans;
%             else
%                 ChansToKeep = [1 2 NumChans-1 NumChans];
%             end
%         end
%
%     else
%         ChansToKeep = 1:size(W,1);
%     end
%
%     NumChans = length(ChansToKeep);
%     W = W(ChansToKeep,:);
%
%
%
%     SpatFilt = W*RefFilt*ArtWeights2;
% %     SpatFilt2 = W2*RefFilt*ArtWeights2;
%     %Combine Spatial Filter and Artifact Matrix into Spatial Filter for online
%     %use.  Here i am first doing artifact rejection, then spatial
%     %filtering.  There is a slight difference with the order.
%
%     %Save the Spatial Filter
% %     dlmwrite(['C:\Documents and Settings\amg5106\My Documents\'...
% %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' Name Session '\'...
% %         Name 'S' Session 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
%      dlmwrite(['C:\Documents and Settings\amg5106\My Documents\'...
%         'Year4\BCI2000 3.0.5\data\' Name Session '\'...
%         Name 'S' Session 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
%
%     DD = SpatFilt*Data;
% elseif BuildClassifier == 0
%     SpatFilt = dlmread([NameBase 'SpatialFilter.txt'],'\t');
%
%
% %     disp('DONT FORGET TO CHANGE THIS BACK')
% %     SpatFilt = eye(size(Data,1));
%
%
%     DD = SpatFilt*Data;
%     NumChans = size(DD,1);
%
%
%
%
% end
%
% %% Process and Organize Data
%
% %%%%
% %Associate Target Codes with Targets

% TargSymb = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N',...
%     'O','P','Q','R','S','T','U','V','W','X','Y','Z','_','0','1','2',...
%     '3','4'};

% SymbSymb = {'Yes','Rest','Wash','Toil','Call','NO','Doc','Nurse','Hus','Fam',...
%     'Vis','Wheel','Walk','Pillow','Blanket','Clothes','Med','Pray','Clerg',...
%     'Drink','Eat','Slip','Pap','Tv','Music','Book','Glass','Breat','Cold',...
%     'Hot','Time','Date'};
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Organize Data into choice and not-choice targets
%
%
% if FreeSpell==0
%     targetstim = double(StimulusCode(logical(StimulusType)));
%     targetstim = targetstim(1:StateDuration:length(targetstim));
%     targetstim = [targetstim(1:(NumSequences+SequencestoRemove)*2:length(targetstim))...
%         targetstim(2:(NumSequences+SequencestoRemove)*2:length(targetstim))];
% else
%     targetstim = repmat([1 2],NumTrials,1); %Put fake target stims in... just
%     %so the code below runs and we group the data for classification.
%     %We will ignore this when classifying later
% end
%
%
%


FaceData.C = cell(1,24);
FaceData.N = cell(1,24);
if Speller == 1 && (Art==0 || Art==2)
    for i = 1:NumTrials
        for j = 1:24 %Loop through face stimuli
            Faceloc = find(str2num(FaceOF{i})==j);
            for k = Faceloc
                tmpSC = str2num(cell2mat((SC{i}(k)))); tmpSR = str2num(cell2mat(SR{i}(k)));
                TargC = find(StimulusCode == tmpSC); TargC = TargC(1:StateDuration:end); TargC = TargC((i-1)*NumSequences+1:i*NumSequences);
                TargR = find(StimulusCode == tmpSR); TargR = TargR(1:StateDuration:end); TargR = TargR((i-1)*NumSequences+1:i*NumSequences);
                tl = c.TextToSpell.Value{:}(i);
                C_ind = find(strcmp(tl,TargSymb));
                for l = 1:NumSequences
                    
                    if k == C_ind
                        FaceData.C{j} = cat(3,FaceData.C{j},DD(:,TargC(l)+range(1)-1:TargC(l)+range(end)-1));
                        FaceData.C{j} = cat(3,FaceData.C{j},DD(:,TargR(l)+range(1)-1:TargR(l)+range(end)-1));
                    else
                        FaceData.N{j} = cat(3,FaceData.N{j},DD(:,TargC(l)+range(1)-1:TargC(l)+range(end)-1));
                        FaceData.N{j} = cat(3,FaceData.N{j},DD(:,TargR(l)+range(1)-1:TargR(l)+range(end)-1));
                    end
                end
            end
        end
    end
else
end

% %Initialize counters for each stimulus code
% Tind = ones(1,NumStimCodes);
% for i = 1:NumTrials
%     for j = 1:NumStimCodes %Loop through stimulus Codes
%         TargLoc = find(StimulusCode == j); TargLoc = TargLoc(1:StateDuration:end);
% 
%         %Now take the locations of the first NumSeq of each of these vectors, and
%         %increment a counter by NumSeq
%         TrialTargLoc{i,j} = TargLoc(Tind(j):Tind(j)+NumSequences-1);
%         Tind(j) = Tind(j)+NumSequences+SequencestoRemove;
% 
%         for ch = 1:NumChans %Loop through channels
%             TrialTargData{i,j}{ch} = [];
% %             figure
%             for s = 1:NumSequences %Loop through sequences per trial, extract data 0-.6 seconds after stimulus
% 
%                 TrialTargData{i,j}{ch} = [TrialTargData{i,j}{ch}; DD(ch,TrialTargLoc{i,j}(s)+range(1)-1:TrialTargLoc{i,j}(s)+range(end)-1)];
% 
% %            plot(StimulusCode); hold on;
% %            plot(TrialTargLoc{i,j}(s)+range(1)-1:TrialTargLoc{i,j}(s)+range(end)-1,...
% %                 DD(TrialTargLoc{i,j}(s)+range(1)-1:TrialTargLoc{i,j}(s)+range(end)-1,ch)'); hold on;
%             end
%         end
%     end
% end
%
% disp(['Time for organizing data ' num2str(toc)]);
%
% % %% Old method (slower)
% % tic
% % oEEGData = [];
% % oChanLabel = [];
% % oClassLabel = [];
% % oAvgEEGData = [];
% % oAvgChanLabel = [];
% % oAvgClassLabel = [];
% % for ch = 1:NumChans %Loop through channels
% %     odataloc{ch} = [];
% %     for i = 1:size(StimCode,2) %Loop through Trials
% %
% %         for j = targetstim(i,:)
% %             oEEGData = [oEEGData; TrialTargData{i,j}{ch}];
% %             oAvgEEGData = [oAvgEEGData; mean(TrialTargData{i,j}{ch},1)];
% %             oAvgChanLabel = [oAvgChanLabel; ch];
% %             oAvgClassLabel = [oAvgClassLabel; 'C'];
% %             oChanLabel = [oChanLabel; repmat(ch,NumSequences,1)];
% %             oClassLabel = [oClassLabel; repmat('C',NumSequences,1)];
% %             odataloc{ch} = [odataloc{ch}; TrialTargLoc{i,j}];
% %         end
% %         for k = find(~ismember(1:NumStimCodes,targetstim(i,:))==1)
% %             oEEGData = [oEEGData; TrialTargData{i,k}{ch}];
% %             oAvgEEGData = [oAvgEEGData; mean(TrialTargData{i,k}{ch},1)];
% %             oAvgChanLabel = [oAvgChanLabel; ch];
% %             oAvgClassLabel = [oAvgClassLabel; 'N'];
% %             oChanLabel = [oChanLabel; repmat(ch,NumSequences,1)];
% %             oClassLabel = [oClassLabel; repmat('N',NumSequences,1)];
% %             odataloc{ch} = [odataloc{ch}; TrialTargLoc{i,k}];
% %         end
% %     end
% % end
% % disp(['Time to rear data (1)' num2str(toc)]);
% tic
%
% %Changed from 'N' and 'C' based to 2 and 1 based.
%
% % EEGData = zeros(NumSequences*NumTrials*NumStimCodes*NumChans,length(range));
% % ChanLabel = zeros(size(EEGData,1),1);
% % ClassLabel = repmat(char(0),size(EEGData,1),1);
% % AvgEEGData = zeros(NumTrials*NumStimCodes*max(ChanLabel),length(range));
% % AvgChanLabel = zeros(size(AvgEEGData,1),1);
% % AvgClassLabel = repmat(char(0),size(AvgEEGData,1),1);
% % for ch = 1:NumChans %Loop through channels
% %     dataloc{ch} = [];
% %     for i = 1:NumTrials %Loop through Trials
% %
% %         for j = 1:NumStimCodes
% %             index = (ch-1)*(NumTrials*NumStimCodes*NumSequences)+(i-1)*...
% %                 (NumStimCodes*NumSequences)+(j-1)*NumSequences+1:(ch-1)*...
% %                 (NumTrials*NumStimCodes*NumSequences)+(i-1)*(NumStimCodes*...
% %                 NumSequences)+(j)*NumSequences;
% %             avg_index = (ch-1)*(NumTrials*NumStimCodes)+(i-1)*...
% %                 (NumStimCodes)+(j-1)+1:(ch-1)*(NumTrials*NumStimCodes)+...
% %                 (i-1)*(NumStimCodes)+(j);
% %             EEGData(index,:) = TrialTargData{i,j}{ch};
% %             AvgEEGData(avg_index,:) = mean(TrialTargData{i,j}{ch},1);
% %
% %
% %             if ismember(j,targetstim(i,:)) %Choice targets = 1
% %              ChanLabel(index) = repmat(ch,NumSequences,1);
% %             ClassLabel(index) = repmat('C',NumSequences,1);
% %             AvgChanLabel(avg_index) = ch;
% %             AvgClassLabel(avg_index) = 'C';
% %             else
% %                ChanLabel(index) = repmat(ch,NumSequences,1);
% %             ClassLabel(index) = repmat('N',NumSequences,1);
% %             AvgChanLabel(avg_index) = ch;
% %             AvgClassLabel(avg_index) = 'N';
% %             end
% %
% %
% %             dataloc{ch}(index) = TrialTargLoc{i,j};
% %         end
% %
% %     end
% % end
%
%
% EEGData = zeros(NumSequences*NumTrials*NumStimCodes*NumChans,length(range));
% ChanLabel = zeros(size(EEGData,1),1);
% ClassLabel = zeros(size(EEGData,1),1);
% AvgEEGData = zeros(NumTrials*NumStimCodes*max(ChanLabel),length(range));
% AvgChanLabel = zeros(size(AvgEEGData,1),1);
% AvgClassLabel = zeros(size(AvgEEGData,1),1);
% for ch = 1:NumChans %Loop through channels
%     dataloc{ch} = [];
%     for i = 1:NumTrials %Loop through Trials
%
%         for j = 1:NumStimCodes
%             index = (ch-1)*(NumTrials*NumStimCodes*NumSequences)+(i-1)*...
%                 (NumStimCodes*NumSequences)+(j-1)*NumSequences+1:(ch-1)*...
%                 (NumTrials*NumStimCodes*NumSequences)+(i-1)*(NumStimCodes*...
%                 NumSequences)+(j)*NumSequences;
%             avg_index = (ch-1)*(NumTrials*NumStimCodes)+(i-1)*...
%                 (NumStimCodes)+(j-1)+1:(ch-1)*(NumTrials*NumStimCodes)+...
%                 (i-1)*(NumStimCodes)+(j);
%
%             EEGData(index,:) = TrialTargData{i,j}{ch};
%             AvgEEGData(avg_index,:) = mean(TrialTargData{i,j}{ch},1);
% %             if ch==1 && i==1 && j==1
% %             disp('Detrending P300 epochs!');
% %             end
% %             EEGData(index,:) = detrend(TrialTargData{i,j}{ch},0);
% %             AvgEEGData(avg_index,:) = detrend(mean(TrialTargData{i,j}{ch},1),0);
% %
%             if ismember(j,targetstim(i,:)) %Choice targets = 1
%              ChanLabel(index) = repmat(ch,NumSequences,1);
%             ClassLabel(index) = ones(NumSequences,1);
%             AvgChanLabel(avg_index) = ch;
%             AvgClassLabel(avg_index) = 1;
%             else
%                ChanLabel(index) = repmat(ch,NumSequences,1);
%             ClassLabel(index) = 2*ones(NumSequences,1);
%             AvgChanLabel(avg_index) = ch;
%             AvgClassLabel(avg_index) = 2;
%             end
%
%
%             dataloc{ch}(index) = TrialTargLoc{i,j};
%         end
%
%     end
% end
%
% if FreeSpell == 0
%     AvgClassLabel = AvgClassLabel';
%     AvgChanLabel = AvgChanLabel';
% elseif FreeSpell == 1
%     clear ClassLabel AvgClassLabel
% end
% disp(['Time to rearrange data (2) ' num2str(toc)]);
%
%
% % seqq = 1;
% % sccc = 10;
% % doop = [TrialTargData{1,sccc}{1}(seqq,:)' TrialTargData{1,sccc}{2}(seqq,:)'...
% %     TrialTargData{1,sccc}{3}(seqq,:)' TrialTargData{1,sccc}{4}(seqq,:)'];
% %
% %
% % figure
% % mmm(1) = subplot(211); plot(doop)
% % mmm(2) = subplot(212); plot(sc11)
% % linkaxes(mmm)
% if Figures_On == 1
% if FreeSpell == 0
%     figure
%     for ch = 1:NumChans
%         subplot(4,5,ch);
%         plot(mean(AvgEEGData(AvgClassLabel==1&AvgChanLabel==ch,:)));
%         hold on;
%         plot(mean(AvgEEGData(AvgClassLabel==2&AvgChanLabel==ch,:)),'r')
%     end
% elseif FreeSpell == 1
%     figure
%     for ch = 1:NumChans
%         subplot(4,5,ch);
%         plot(mean(AvgEEGData(AvgChanLabel==ch,:)));
%     end
% end
% end
%
%
% %% Data Reduction - for LDA Classifiers
% tic
% if ReducT == 1
%     %Decimate Data - use decimated data in classifier
% %
% %     clear dEEGData dAvgEEGData
% %     for dd = 1:size(EEGData,1)
% %         dEEGData(dd,:) = decimate(EEGData(dd,range),8);
% %     end
% %     for dd = 1:size(AvgEEGData,1)
% %         dAvgEEGData(dd,:) = decimate(AvgEEGData(dd,range),8);
% %     end
% %     LDAtimes = round(decimate(range,8));
% %
% %
% %     inpch = [];
% %     if Avv == 0
% %         TrainD = [];
% %         for ch = 1:max(ChanLabel)
% %             TrainD = [TrainD double(dEEGData(ChanLabel==ch,:))];
% %             inpch = [inpch ch*ones(1,size(dEEGData,2))];
% %         end
% %     else
% %         AvgTrainD = [];
% %         for ch=1:max(AvgChanLabel)
% %             AvgTrainD = [AvgTrainD double(dAvgEEGData(AvgChanLabel==ch,:))];
% %             inpch = [inpch ch*ones(1,size(dAvgEEGData,2))];
% %         end
% %     end
%
% elseif ReducT == 2
%     %     %Find timepoints where averages are significantly different - pick from
%     %     %these to be the features of the LDA.
%     %     for ch = 1:max(ChanLabel)
%     %         %             figure
%     %         Cch = ChanLabel==ch;
%     %         Cn = ClassLabel=='N';
%     %         Cc = ClassLabel=='C';
%     %         nloc = Cch&Cn;
%     %         cloc = Cch&Cc;
%     %         %             plot(mean(EEGData(nloc,range)),'r');
%     %         %             hold on;
%     %         %             plot(mean(EEGData(cloc,range)),'b');
%     %         for i = 1:length(range)
%     %             [tstat(ch,i) pval(ch,i)] =ttest(EEGData(cloc,range(i)),mean(EEGData(nloc,range(i))));
%     %         end
%     %         [~, LDAtimes(ch)]=min(pval(ch,:));
%     %         dEEGData(Cch,:) = EEGData(Cch,LDAtimes(ch));
%     %
%     %
%     %         %             figure
%     %         Cch = AvgChanLabel==ch;
%     %         Cn = AvgClassLabel=='N';
%     %         Cc = AvgClassLabel=='C';
%     %         nloc = Cch&Cn;
%     %         cloc = Cch&Cc;
%     %         %             plot(mean(AvgEEGData(nloc,range)),'r');
%     %         %             hold on;
%     %         %             plot(mean(AvgEEGData(cloc,range)),'b');
%     %         for i = 1:length(range)
%     %             [Atstat(ch,i) Apval(ch,i)] =ttest(AvgEEGData(cloc,range(i)),mean(AvgEEGData(nloc,range(i))));
%     %         end
%     %         [~, AvgLDAtimes(ch)]=min(Apval(ch,:));
%     %         dAvgEEGData(Cch,:) = AvgEEGData(Cch,AvgLDAtimes(ch));
%     %
%     %     end
%     %
%     %     inpch = [];
%     %     if Avv == 0
%     %         TrainD = [];
%     %         for ch = 1:max(ChanLabel)
%     %             TrainD = [TrainD double(dEEGData(ChanLabel==ch,:))];
%     %             inpch = [inpch ch*ones(1,size(dEEGData,2))];
%     %         end
%     %     elseif Avv == 1
%     %         AvgTrainD = [];
%     %         for ch=1:max(AvgChanLabel)
%     %             AvgTrainD = [AvgTrainD double(dAvgEEGData(AvgChanLabel==ch,:))];
%     %             inpch = [inpch ch*ones(1,size(dAvgEEGData,2))];
%     %         end
%     %     end
%
% elseif ReducT == 3 %%Low pass filter, downsample to 20Hz
%
%     %lowpass filter at 20 Hz
%     LPv = 20;
%     [bb2 aa2] = butter(5,LPv/(fs/2),'low');
%
%     inpch = [];
%         x_filtered = filtfilt(bb2,aa2,EEGData');
%         x_down = downsample(x_filtered,fix(fs/(LPv)))';
%         LDAtimes = downsample(range,fix(fs/(LPv)));
%         TrainD = [];
%         for ch = 1:max(ChanLabel)
%             TrainD = [TrainD double(x_down(ChanLabel==ch,:))];
%             inpch = [inpch ch*ones(1,size(x_down,2))];
%         end
%
%
%
%         x_filtered = filtfilt(bb2,aa2,AvgEEGData');
%         x_down = downsample(x_filtered,fix(fs/(LPv)))';
%
% %                 figure
% %                 plot(AvgEEGData(142,:));
% %                 hold on;
% %                 plot(x_filtered(:,142),'g')
% %                 plot(LDAtimes,x_down(142,:),'r');
%         AvgTrainD = [];
%         for ch=1:max(AvgChanLabel)
%             AvgTrainD = [AvgTrainD double(x_down(AvgChanLabel==ch,:))];
%             inpch = [inpch ch*ones(1,size(x_down,2))];
%         end
%
%
%
% elseif ReducT == 0
%
% %     inpch = [];
% %     if Avv == 0
% %         TrainD = [];
% %         for ch = 1:max(ChanLabel)
% %             TrainD = [TrainD double(EEGData(ChanLabel==ch,:))];
% %             inpch = [inpch ch*ones(1,size(EEGData,2))];
% %         end
% %     else
% %         AvgTrainD = [];
% %         for ch=1:max(AvgChanLabel)
% %             AvgTrainD = [AvgTrainD double(AvgEEGData(AvgChanLabel==ch,:))];
% %             inpch = [inpch ch*ones(1,size(AvgEEGData,2))];
% %         end
% %     end
% %
% %     LDAtimes = range;
%
% end
%
% disp(['Time to data reduce ' num2str(toc)]);
%
% %% P300 Sanity Check (3)
% if Figures_On == 1
% if FreeSpell == 0
%     figure
%     for ch = 1:NumChans
%         subplot(4,5,ch);
%         plot(LDAtimes,AvgTrainD(AvgClassLabel(AvgChanLabel==1)==2,...
%             (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes)),'r'); hold on;
%         plot(LDAtimes,AvgTrainD(AvgClassLabel(AvgChanLabel==1)==1,...
%             (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes)),'b');
%     end
% elseif FreeSpell == 1
% end
% end
% %% Classification
% if BuildClassifier == 1
%
%
%         %Open the parameter file for writing
%
%     fpm = fopen(['C:\Documents and Settings\amg5106\My Documents\Year4\'...
%         'BCI2000 3.0.5\parms\ALSP3\FINAL_CrossSpeller_test_22_3EOG_2amp.prm'],'r');
%     fpm2=textscan(fpm,'%s','delimiter','\n');
%     fclose(fpm);
%     for ll = 1:length(fpm2{1})
%         LinInd(ll) = ~isempty(strfind(fpm2{1}{ll},'Filtering:Linear'));
%         SNInd(ll) = ~isempty(strfind(fpm2{1}{ll},'SubjectName'));
%         SEInd(ll) = ~isempty(strfind(fpm2{1}{ll},'SubjectSession'));
%         SFInd(ll) = ~isempty(strfind(fpm2{1}{ll},'SpatialFilter matrix'));
%     end
%     LinInd = find(LinInd==1);
%     SNInd = find(SNInd==1);
%     SEInd = find(SEInd==1);
%     SFInd = find(SFInd==1);
%     SBegInd = strfind(fpm2{1}{SNInd},'='); SEndInd = strfind(fpm2{1}{SNInd},'Name %');
%     SEBegInd = strfind(fpm2{1}{SEInd},'='); SEEndInd = strfind(fpm2{1}{SEInd},'% %');
%     SFBegInd = strfind(fpm2{1}{SFInd},'='); SFEndInd = strfind(fpm2{1}{SFInd},'//');
%     BegInd = strfind(fpm2{1}{LinInd},'}'); EndInd = strfind(fpm2{1}{LinInd},'//');
%     BegInd2 = strfind(fpm2{1}{LinInd},'= '); EndInd2 = strfind(fpm2{1}{LinInd},' {');
%
%
%
%     %      NumCrossVal = 10;
%     %
%     %         for ccc = 1:NumCrossVal
%     %             test_ind = (ccc-1)*size(AvgTrainD,1)/NumCrossVal+1:...
%     %                 ccc*size(AvgTrainD,1)/NumCrossVal;
%     %             train_ind = find(~ismember(1:size(AvgTrainD),test_ind));
%     %             cv_train = AvgTrainD(train_ind,:);
%     %             cv_test = AvgTrainD(test_ind,:);
%     %             cv_cl = AvgClassLabel(train_ind);
%     %
%     %
%     %             if Cfier == 1
%     %                 %Stepwise Linear Regression
%     %                 if Avv == 0
%     %                     [B, SE, PVAL, INMODEL] = stepwisefit(TrainD,ClassLabel(ChanLabel==1));
%     %                 else
%     %                     [B, SE, PVAL, INMODEL] = stepwisefit(cv_train,cv_cl);
%     %                 end
%     %                 SigLocs = find(INMODEL==1);
%     %                 Classweight = -B;
%     %
%     %             elseif Cfier == 2
%     %                 %LDA
%     %                 %Classifier Built from individual trials (dont think this is what I
%     %                 %want to do)
%     %                 if Avv == 0
%     %                     [class,err,post,logp,coeff] = classify(TrainD,TrainD,ClassLabel(ChanLabel==1));
%     %                 elseif Avv == 1
%     %                     [class,err,post,logp,coeff] = classify(AvgTrainD,AvgTrainD,AvgClassLabel(AvgChanLabel==1));
%     %                 end
%     %                 Classweight = coeff(1,2).linear;
%     %                 SigLocs = find(Classweight~=0);
%     %
%     %
%     %             elseif Cfier == 3
%     %                 % Steve's LDA
%     %                 %This produces a crazy gamma matrix similar to the regular LDA, which
%     %                 %doesnt work.  Figure out how to make regular LDA output a useful gamma
%     %                 %matrix and then come back to this
%     %
%     %                 if Avv == 0
%     %                     cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
%     %                         'ALS Project\Part2 - P300\P300Classifier\LDA Code\Redo'])
%     %                     [Gam,z] = runLDA('Andrew',TrainD,ClassLabel(ChanLabel==1));
%     %                     cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
%     %                         'ALS Project\Part2 - P300\P300Classifier'])
%     %                 elseif Avv == 1
%     %                     cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
%     %                         'ALS Project\Part2 - P300\P300Classifier\LDA Code\Redo'])
%     %                     [Gam,z] = runLDA('Andrew',AvgTrainD,AvgClassLabel(AvgChanLabel==1));
%     %                     cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
%     %                         'ALS Project\Part2 - P300\P300Classifier'])
%     %                 end
%     %                 Classweight = Gam;
%     %                 SigLocs = find(Gam~=0);
%     %             end
%
%     if Cfier == 1
%         %Stepwise Linear Regression
%         if Avv == 0
%             [B, SE, PVAL, INMODEL] = stepwisefit(TrainD,...
%                 ClassLabel(ChanLabel==1));%,'penter',.005,'premove',.01);
%         else
%             [B, SE, PVAL, INMODEL] = stepwisefit(AvgTrainD,...
%                 AvgClassLabel(AvgChanLabel==1));%,'penter',.005,'premove',.01);
%         end
%         Classweight = -B(INMODEL==1);
%         SigLocs = find(INMODEL==1);
%
% %         [~, keepind] = sort(PVAL);
% %         if length(keepind)>10
% %             keepind = keepind(1:10);
% %         end
% %         SigLocs = keepind;
% %         TrainD2 = TrainD(:,SigLocs);
% %         [B, SE, PVAL, INMODEL, STATS] = stepwisefit(TrainD2,ClassLabel(ChanLabel==1));
% %         Classweight = -B(INMODEL==1);
% %         SigLocs = SigLocs(INMODEL==1);
%
%
%     elseif Cfier == 2
%         %LDA
%         %Classifier Built from individual trials (dont think this is what I
%         %want to do)
%         if Avv == 0
%             [class,err,post,logp,coeff] = classify(TrainD,TrainD,ClassLabel(ChanLabel==1));
%         elseif Avv == 1
%             [class,err,post,logp,coeff] = classify(AvgTrainD,AvgTrainD,AvgClassLabel(AvgChanLabel==1));
%         end
%         Classweight = coeff(1,2).linear;
%         SigLocs = find(Classweight~=0);
%
%
%     elseif Cfier == 3
%         % Steve's LDA
%         %This produces a crazy gamma matrix similar to the regular LDA, which
%         %doesnt work.  Figure out how to make regular LDA output a useful gamma
%         %matrix and then come back to this
%
%         if Avv == 0
%             cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
%                 'ALS Project\Part2 - P300\P300Classifier\LDA Code\Redo'])
%             [Gam,z,SigLocs] = runLDA2([],TrainD,ClassLabel(ChanLabel==1),10);
%             cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
%                 'ALS Project\Part2 - P300\P300Classifier'])
%         elseif Avv == 1
%             cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
%                 'ALS Project\Part2 - P300\P300Classifier\LDA Code\Redo'])
%             [Gam,z,SigLocs] = runLDA2([],AvgTrainD,AvgClassLabel(AvgChanLabel==1),10);
%             cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
%                 'ALS Project\Part2 - P300\P300Classifier'])
%         end
%         %Sometimes the polarity is reversed on the output.  We always want
%         %the first class label (p300) to be greater
%         if mean(z(1:length(z)/2))-mean(z(length(z)/2+1:end))>0
%             Classweight(SigLocs) = Gam;
%         else
%             Classweight(SigLocs) = -Gam;
%         end
%     end
%
%      %Construct Classifier Matrix
%     CM = zeros(length(SigLocs),4);
%     if isempty(CM)
%         disp('No significant features!');
%         CM = NaN; choicesums = NaN; notchoicesums = NaN; targetstim = NaN;
%         WORDD = NaN; allsums = NaN;
%         skiptest = 1;
%         return
%     else
%         CM(:,1) = inpch(SigLocs); %Input Channel
%         Ctms = mod(SigLocs,length(LDAtimes));  Ctms(Ctms==0)=length(LDAtimes);
%         CM(:,2) = LDAtimes(Ctms); %Input element (sample #)
%         CM(:,3) = 1; %Output channel
%         CM(:,4) = Classweight; %Weight
%        % Changed from CM(:,4) = Classweight(SigLocs);  to reflect changes
%        % to SW classifier that limits features to a specified number.
%
%     ClassTypes = {'SW','L','SS'};
%     dlmwrite(['C:\Documents and Settings\amg5106\My Documents\'...
%         'Year4\BCI2000 3.0.5\data\' Name Session '\'...
%         Name 'S' Session 'R' NewRunName '_' ClassTypes{Cfier} 'Classifier.txt'],CM,'delimiter','\t')
%     skiptest = 0;
%
%
%     %save in parameter file
%     prmcls = num2str(reshape(CM',1,size(CM,1)*size(CM,2)));
%     %Change Subject Name
%     fpm2{1}{SNInd} = [fpm2{1}{SNInd}(1:SBegInd+1) Name...
%         fpm2{1}{SNInd}(SEndInd-1:end)];
%     %Change Session
%     fpm2{1}{SEInd} = [fpm2{1}{SEInd}(1:SEBegInd+1) Session...
%         fpm2{1}{SEInd}(SEEndInd-1:end)];
%     %Change Spatial Filter
%     prmsf = [num2str(size(SpatFilt))  ' ' num2str(reshape(SpatFilt',1,...
%         size(SpatFilt,1)*size(SpatFilt,2)))];
%     fpm2{1}{SFInd} = [fpm2{1}{SFInd}(1:SFBegInd+1) prmsf...
%         fpm2{1}{SFInd}(SFEndInd-1:end)];
%     %Change Classifier
%     fpm2{1}{LinInd} = [fpm2{1}{LinInd}(1:BegInd2+1) num2str(size(CM,1)) ...
%         fpm2{1}{LinInd}(EndInd2:BegInd+1) prmcls...
%         fpm2{1}{LinInd}(EndInd-1:end)];
%
%     %Save new parameter file
%     dlmcell(['C:\Documents and Settings\amg5106\My Documents\Year4\'...
%         'BCI2000 3.0.5\data\' Name Session '\ParamFile'...
%         Name 'S' Session 'R' NewRunName '_' ClassTypes{Cfier} '.prm'],fpm2{1});
%
%
%     end
%
%
% else %Testrun - do not build classifier, prompt for a previously built one
%      ClassTypes = {'SW','L','SS'};
%     %Load Classifier Data
%     disp(['Loading classifier: ' NameBase ClassTypes{Cfier} 'Classifier.txt']);
%     CM = dlmread([NameBase ClassTypes{Cfier} 'Classifier.txt'],'\t');
%
%     %Unnecessary for P300
% %     NormalizerGain = CM(end-1,1);
% %     NormalizerOffset = CM(end,1);
% %     CM = CM(1:end-2,:);
%
% % Name = 'AndrewGrid-gUSB';
% % CSession = '002';
% % CRun = '02';
% % CM = dlmread(['C:\Documents and Settings\amg5106\My Documents\'...
% %         'Year2\Research\BCI2000\BCI2000sourcetree\data\AndrewGrid-gUSB002\'...
% %         'AndrewGrid-gUSBS002R01_SWClassifier.txt'],'\t');
%
% end
%
%
%
%      %Looking at the expected classification for trial averages
%     for ii = 1:size(AvgTrainD,1)
%         cursA(ii) = 0;
%         for jj = 1:size(CM,1)
%             ffr = (CM(jj,1)-1)*length(LDAtimes)+1:CM(jj,1)*length(LDAtimes);
%             ffrx = LDAtimes==CM(jj,2);
%             cursA(ii) = cursA(ii) + CM(jj,4)*AvgTrainD(ii,ffr(ffrx));
%         end
%     end
%     %What is optimal decision boundary?
%     pcurs = cursA(AvgClassLabel(AvgChanLabel==1)==1);
%     ncurs = cursA(AvgClassLabel(AvgChanLabel==1)==2);
%     k=1;
%     rangeC = min(ncurs):.01:max(pcurs);
%     for thr = rangeC %Cycle through possible thresholds
%         pp = sum(pcurs>thr)/length(pcurs);
%         pn = sum(ncurs<thr)/length(ncurs);
%         psums(k) = pp+pn; k=k+1;
%     end
%     [ma mi] = max(psums)
%     %Set the gain and the offset of the Normalizer
%     NormalizerOffset = rangeC(mi);
%     NormalizerGain = 1/var(cursA);
%
%     %Looking at the expected classification for individual segments
%     for ii = 1:size(TrainD,1)
%         curs(ii) = 0;
%         for jj = 1:size(CM,1)
%             ffr = (CM(jj,1)-1)*length(LDAtimes)+1:CM(jj,1)*length(LDAtimes);
%             ffrx = LDAtimes==CM(jj,2);
%             curs(ii) = curs(ii) + CM(jj,4)*TrainD(ii,ffr(ffrx));
%         end
%     end
%     for tt = 1:NumTrials*NumStimCodes
%         trlindex = (tt-1)*NumSequences+1:(tt)*NumSequences;
%         cursA2(tt) = mean(curs(trlindex));
%     end
%
%     if Figures_On == 1
%     figure
%     subplot(211);
%     hist(cursA(AvgClassLabel(AvgChanLabel==1)==1),10);
%     hold on;
%     hist(cursA(AvgClassLabel(AvgChanLabel==1)==2),10);
%     h = findobj(gca,'Type','patch');
%     set(h(1),'FaceColor','r','EdgeColor','w','FaceAlpha',.5);
%     %weird that this labeling of h(1) and h(2) is backwards
%     set(h(2),'FaceColor','b','EdgeColor','w','FaceAlpha',.5);
%     line([NormalizerOffset NormalizerOffset],[0 20],'Color','k','LineWidth',3)
%     subplot(212);
%     hist(curs(ClassLabel(ChanLabel==1)==1),100);
%     hold on;
%     hist(curs(ClassLabel(ChanLabel==1)==2),100);
%     h = findobj(gca,'Type','patch');
%     set(h(1),'FaceColor','r','EdgeColor','w','FaceAlpha',.5);
%     set(h(2),'FaceColor','b','EdgeColor','w','FaceAlpha',.5);
%     end
%
%
%
%     %Instead of the line on the histogram plot, do an ROC curve
%     %If we consider "truth" as being a left trial
%     k=1;
%     for thr = rangeC
%         tp(k) = sum(pcurs>thr); fn(k) = sum(pcurs<thr);
%         tn(k) = sum(ncurs<thr); fp(k) = sum(ncurs>thr); k = k+1;
%     end
%     tpr = tp./(tp+fn);
%     fpr = fp./(fp+tn);
%     auckinda = tpr.*(1-fpr);
%     [~, emaxi] = max(auckinda);
%     rangeC(emaxi);
%     lineslope = 1;
%     rangeL = 0:.01:1;
%     deg45line = lineslope*(rangeL-fpr(emaxi))+tpr(emaxi);
%     if Figures_On == 1
%     figure
%     plot(fpr,tpr); hold on
%     plot(rangeL,deg45line,'r'); xlim([0 1]); ylim([0 1]);
%     end
%     %Do the two methods produce the same result?
%     if mi==emaxi
%         disp('Yay!');
%     else
%         disp('Different :(');
%     end
%
%
% %     CM(size(CM,1)+1,:) = [NormalizerGain 0 0 0];
% %     CM(size(CM,1)+1,:) = [NormalizerOffset 0 0 0];
%
% %% View Results
%
% if FreeSpell == 0 %Find the classified choice and not choice data
%     %Absolutely no difference in whether the data is classified using the
%     %Averaged or non averaged data. The result of choice and notchoice are
%     %the same.  FOr this reason we remove the functionality of switching
%     %between Avv=0 and Avv=1
%
%
%
%
%     %     if Avv == 0
%     %         ChoiceOut = []; NotChoiceOut = [];
%     %         % ChoiceOut = coeff(1,2).const*ones(sum(ClassLabel(ChanLabel==1)=='C'),1);
%     %         % NotChoiceOut = coeff(1,2).const*ones(sum(ClassLabel(ChanLabel==1)=='N'),1);
%     %         ChoiceL = ClassLabel(1:NumStimCodes*NumSequences*NumTrials)=='C';
%     %         NChoiceL = ClassLabel(1:NumStimCodes*NumSequences*NumTrials)=='N';
%     %         %Plot the classification results with the generated classifier
%     %         for i = 1:size(ClassifierMatrix,1)
%     %             tform = EEGData((ClassifierMatrix(i,1)-1)*(NumStimCodes*NumSequences*NumTrials)+1:...
%     %                 ClassifierMatrix(i,1)*(NumStimCodes*NumSequences*NumTrials),ClassifierMatrix(i,2))...
%     %                 *ClassifierMatrix(i,4);
%     %             ChoiceOut = [ChoiceOut tform(ChoiceL)];
%     %             NotChoiceOut = [NotChoiceOut tform(NChoiceL)];
%     %         end
%     %         %Consider for each trial what the sum of the choice and not choice
%     %         %sequences are
%     %         for i = 1: NumTrials
%     %             choicesums{i} = sum(reshape(sum(ChoiceOut(i*size(targetstim,2)*NumSequences...
%     %             -(size(targetstim,2)*NumSequences-1):i*size(targetstim,2)*NumSequences,:),2),...
%     %             NumSequences,size(targetstim,2)),1);
%     %
%     %             notchoicesums{i} = sum(reshape(sum(NotChoiceOut(i*...
%     %                 (NumStimCodes-size(targetstim,2))*NumSequences...
%     %             -((NumStimCodes-size(targetstim,2))*NumSequences-1):i*...
%     %             (NumStimCodes-size(targetstim,2))*NumSequences,:),2),...
%     %             NumSequences,NumStimCodes-size(targetstim,2)),1);
%     %
%     %             figure
%     %             hist(notchoicesums{i});hold on;
%     %             hist(choicesums{i});
%     %             h = findobj(gca,'Type','patch');
%     %             set(h(1),'FaceColor','r');
%     %         end
%     %
%     %         figure
%     %         plot(range/fs,EEGData(ClassLabel=='N'&ChanLabel==Clbl,:)','r');
%     %         hold on;
%     %         plot(range/fs,EEGData(ClassLabel=='C'&ChanLabel==Clbl,:)','b');
%     %
%     %         figure
%     %         errorbar(mean(EEGData(ClassLabel=='C'&ChanLabel==Clbl,range),1),...
%     %             std(EEGData(ClassLabel=='C'&ChanLabel==Clbl,range),[],1)/...
%     %             sqrt(size(EEGData(ClassLabel=='C'&ChanLabel==Clbl,range),1))); hold on;
%     %         errorbar(mean(EEGData(ClassLabel=='N'&ChanLabel==Clbl,range),1),...
%     %             std(EEGData(ClassLabel=='N'&ChanLabel==Clbl,range),[],1)/...
%     %             sqrt(size(EEGData(ClassLabel=='N'&ChanLabel==Clbl,range),1)),'r');
%     %         xlim([0 length(range)])
%     %    else
%     AllOut = [];
%     %Plot the classification results with the generated classifier
%     for i = 1:size(CM,1)
%         tform = AvgEEGData((CM(i,1)-1)*(NumStimCodes*NumTrials)+1:...
%             CM(i,1)*(NumStimCodes*NumTrials),CM(i,2)-range(1)+1)....
%             *CM(i,4);
%
%         AllOut = [AllOut tform];
%     end
%     for i = 1: NumTrials
%
%         allsums{i} = sum(AllOut((i-1)*NumStimCodes+1:i*NumStimCodes,:),2)';
%
%
%         cclabel = targetstim(i,:);
%         nnlabel = find(~ismember(1:NumStimCodes,targetstim(i,:))==1);
%         choicesums{i} = allsums{i}(1,cclabel);
%         notchoicesums{i} = allsums{i}(1,nnlabel);
%
%
%         %         figure
%         %         hist(notchoicesums{i});hold on;
%         %         hist(choicesums{i});
%         %         h = findobj(gca,'Type','patch');
%         %         set(h(1),'FaceColor','r');
%         if Figures_On == 1
%             warning('off','all')
%             figure
%             scatter(cclabel,choicesums{i},'r');
%             hold on
%             scatter(nnlabel,notchoicesums{i},'b');
%             for scc = 1:NumStimCodes
%                 if Speller==1
%                     inR = strcmp(num2str(scc),SR{1,i});
%                     if sum(inR)>0
%                         incodesR = find(strcmp(num2str(scc),SR{1,i}));
%                         for icR = 1:length(incodesR)
%                             text(double(scc),icR*2,TargSymb{incodesR(icR)})
%                         end
%                     else
%                         incodesC = find(strcmp(num2str(scc),SC{1,i}));
%                         for icC = 1:length(incodesC)
%                             text(double(scc),icC*2,TargSymb{incodesC(icC)})
%                         end
%                     end
%                 elseif Speller == 2 || Speller == 3
%                     inR = scc==SR{i};
%                     if sum(inR)>0
%                         incodesR = find(scc==SR{i});
%                         for icR = 1:length(incodesR)
%                             text(double(scc),icR*2,TargSymb{incodesR(icR)})
%                         end
%                     else
%                         incodesC = find(scc==SC{i});
%                         for icC = 1:length(incodesC)
%                             text(double(scc),icC*2,TargSymb{incodesC(icC)})
%                         end
%                     end
%                 end
%             end
%             ylim([min([choicesums{i} notchoicesums{i} 2]) ...
%                 max([choicesums{i} notchoicesums{i} ...
%                 2*max([length(incodesR) length(incodesC)])])]);
%             warning('on','all')
%         end
%     end
%
%
%     if Figures_On == 1
%         figure
%         for Clbl = 1:NumChans
%             subplot(4,5,Clbl);
%             plot(range/fs,AvgEEGData(AvgClassLabel==2&AvgChanLabel==Clbl,:)','r');
%             hold on;
%             plot(range/fs,AvgEEGData(AvgClassLabel==1&AvgChanLabel==Clbl,:)','b');
%         end
%
%         figure
%         for Clbl = 1:NumChans
%             subplot(4,5,Clbl); errorbar(mean(AvgEEGData(AvgClassLabel==1&AvgChanLabel==Clbl,:),1),...
%                 std(AvgEEGData(AvgClassLabel==1&AvgChanLabel==Clbl,:),[],1)/...
%                 sqrt(size(AvgEEGData(AvgClassLabel==1&AvgChanLabel==Clbl,:),1))); hold on;
%             errorbar(mean(AvgEEGData(AvgClassLabel==2&AvgChanLabel==Clbl,:),1),...
%                 std(AvgEEGData(AvgClassLabel==2&AvgChanLabel==Clbl,:),[],1)/...
%                 sqrt(size(AvgEEGData(AvgClassLabel==2&AvgChanLabel==Clbl,:),1)),'r');
%             xlim([0 length(range)]);
%         end
%     end
%     %     figure
%     %     hist(sum(NotChoiceOut,2),20);hold on;
%     %     hist(sum(ChoiceOut,2),20);
%     %     h = findobj(gca,'Type','patch');
%     %     set(h(1),'FaceColor','r');
%     %
% else % In Free Spelling mode
%
%     %     if Avv == 0
%     %         ChoiceOut = [];
%     %         for i = 1:size(ClassifierMatrix,1)
%     %             tform = EEGData((ClassifierMatrix(i,1)-1)*(NumStimCodes*NumSequences*NumTrials)+1:...
%     %                 ClassifierMatrix(i,1)*(NumStimCodes*NumSequences*NumTrials),ClassifierMatrix(i,2))...
%     %                 *ClassifierMatrix(i,4);
%     %             ChoiceOut = [ChoiceOut tform];
%     %         end
%     %         for i = 1: NumTrials
%     %             choicesums{i} = sum(reshape(sum(ChoiceOut(i*NumStimCodes*NumSequences...
%     %                 -(NumStimCodes*NumSequences-1):i*NumStimCodes*NumSequences,:),2),...
%     %                 NumSequences,NumStimCodes),1);
%     %             notchoicesums{i} = [];
%     %         end
%     %     else
%     AllOut = [];
%     for i = 1:size(CM,1)
%         tform = AvgEEGData((CM(i,1)-1)*(NumStimCodes*NumTrials)+1:...
%             CM(i,1)*(NumStimCodes*NumTrials),CM(i,2))....
%             *CM(i,4);
%         AllOut = [AllOut tform];
%     end
%     for i = 1: NumTrials
%         choicesums{i} = sum(AllOut(i*NumStimCodes...
%             -(NumStimCodes-1):i*NumStimCodes,:),2)';
%         notchoicesums{i} = [];
%     end
% end
% % end
%
%
%
% %% Trial Choice
% for i = 1:NumTrials
%     %Sort the individual stimulus code sums
%     [TrialVal(i,:) TrialChoice(i,:)] = sort([choicesums{i} notchoicesums{i}],'descend');
%
%     %Find the actual stimulus codes associated with these
%     trlcode = [targetstim(i,:) find(~ismember(1:NumStimCodes,targetstim(i,:))==1)];
%
%     ll=1;
%     clear trlchc chcsums
%     for jj = 1:size(TrialChoice,2) %Loop through code 1
%         for kk = jj+1:size(TrialChoice,2) %Loop through code 2
%             %For each pair of codes, find the pairing of actual codes
%             trlchc(ll,:) = trlcode(TrialChoice(i,[jj kk]));
%             %Find the sum of these codes
%             chcsums(ll,:) = sum(TrialVal(i,[jj,kk]));
%             ll = ll+1;
%         end
%     end
%     [~, chcsumsort] = sort(chcsums,1,'descend');
%
%     matchfind = 0;
%     while matchfind == 0
%         for jj = chcsumsort'
%             jj;
%
%             if Speller == 1
%                 if ismember(trlchc(jj,1),cellfun(@str2num,SR{i}))
%                     mtch = find(cellfun(@str2num,SR{i})==trlchc(jj,1)&...
%                         cellfun(@str2num,SC{i})==trlchc(jj,2));
%                 else
%                     mtch = find(cellfun(@str2num,SC{i})==trlchc(jj,1)&...
%                         cellfun(@str2num,SR{i})==trlchc(jj,2));
%                 end
%             elseif Speller == 2 || Speller == 3
%                 if ismember(trlchc(jj,1),SR{i})
%                     mtch = find(SR{i}==trlchc(jj,1)&...
%                         SC{i}==trlchc(jj,2));
%                 else
%                     mtch = find(SC{i}==trlchc(jj,1)&...
%                         SR{i}==trlchc(jj,2));
%                 end
%             end
%
%             if ~isempty(mtch)
%                 matchfind = 1;
%                 StimChoice(i,:) = trlchc(jj,:);
%                 TargChoice(i) = mtch;
%
%                 break;
%             end
%         end
%     end
%
%
% end
%
% WORDD = TargSymb(TargChoice)
% WORDD = sum(cell2mat(WORDD) == TTSpell);
% SYMBB = SymbSymb(TargChoice)

fprintf(['Plot ' num2str(toc) '\n\n']);

end
