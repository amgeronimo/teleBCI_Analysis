
%Displays raw and processed P300 data

function [Spec, BSpec, ff] = ...
    ShowMeTheP300(Name,Session,Run,RerefVal,Art,DoCSP)

tic
Figures_On = 0;


currdir = cd;
runloc = currdir(1);

SequencestoRemove = 0;  %For 20 flashes, there are 10 sequences.



% EOGloc  = 14:16;
% EEGloc = 1:13;
% EEGloc = 1:5;
% EOGloc = [];
% EEGloc = 1:14;
% EOGloc = 15:16;
% EOGloc = [];
EEGloc = 1:8;
EOGloc = [];

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
 
    %disp(['Shift Data ' num2str(toc)]);
    
    
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




%% Artifact Rejection  ---removed this

DD = Data;
BD = [];


%% Bandpass filter between 1 and 60, then compute spectra
bpr = [1 60];
[bb,aa] = butter(4,bpr/(fs/2));
ff = [];

if size(DD,2)<512
    Spec = [];
else
    DD2 = filtfilt(bb,aa,DD')';
    clear Spec
    for cc = 1:size(DD2,1)
        [Spec(:,cc), ff] = pwelch(DD2(cc,:),fs*2,fs/2,fs*2,fs);
    end
end

if size(BD,2)<512
    BSpec = [];
else
    BD2 = filtfilt(bb,aa,BD')';
    clear BSpec
    for cc = 1:size(BD2,1)
        [BSpec(:,cc), ff] = pwelch(BD2(cc,:),fs*2,fs/2,fs*2,fs);
    end
end



%     %Same as Chronux Multitaper Below
%     [E,V] = dpss(length(DD2(cc,:)),2);
%     h = spectrum.mtm(E,V);    % Specify DPSS and concentrations
%     % when creating the MTM spectral estimator.
%     S3 = psd(h,DD2(cc,:),'Fs',fs);
%     Spec3(cc,:) = S3.Data; ff3 = S3.Frequencies;
    
    
% params.tapers = [3 5];
% params.Fs = fs;
% [Spec2,ff2] = mtspectrumc(DD2',params);
% [BSpec2,fb2] = mtspectrumc(BD2',params);


% ch = 10;
% figure
% plot((1:size(DD,2))/256,DD(ch,:)');
% hold on
% plot((1:size(DD,2))/256,DD2(ch,:)');


% figure;
% kk = 1;
% for i = [1 10 19]
%     xx = subplot(3,1,kk); plot(ff,Spec(:,i)); hold on
%     plot(ff,BSpec(:,i));
%     kk=kk+1;
% end
% linkaxes(xx)









end
