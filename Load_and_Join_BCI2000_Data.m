function [Data, State, Var, InputText, NewRunName, TTSpell,Art,c] = ...
    Load_and_Join_BCI2000_Data(Name,Sessions,Runs,Var,Art)
%% Load Data
% currdir = cd;
% runloc = currdir(1);
Data = [];
State.StimulusCode = [];
State.StimulusType = [];
InputText = [];
State.SourceTime = [];
State.StimulusTime = [];
State.SelectedStimulus = []; % For audio speller
NewRunName = [];
TTSpell = [];
State.GazeX = [];
State.GazeY = [];
State.SequencePhase = [];
Var.AllTrials = [];
Var.AppPos = [];


kk = 1; max_sess = 0;
for ss = 1:length(Sessions)
    Session = Sessions{ss};
    curr_sess = str2double(Session);
    if curr_sess>max_sess
        max_sess = curr_sess;
        Var.max_sess_i = ss;
    end
    for rr = 1:length(Runs)
        %     [a b c d] = load_bcidat(['C:\Documents and Settings\amg5106\My Documents\'...
        %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' ...
        %         Name Session '\' Name 'S' Session 'R' Run{rr} '.dat'],'-calibrated');
        try
            [a b c d] = load_bcidat(['data\' ...
                Name Session '\' Name 'S' Session 'R' Runs{rr} '.dat'],'-calibrated');
            %disp(['Time to load data ' num2str(toc)]);
        catch
            [a b c d] = load_bcidat(['P:\ALS Proj Data\' ...
                Name Session '\' Name 'S' Session 'R' Runs{rr} '.dat'],'-calibrated');
            %disp(['Time to load data ' num2str(toc)]);
        end
        
     
        
        
        if isfield(c,'SineChannelX') %overwrite headset type if Signal Generator
            Var.capt = '_siggen';
        end
        
        
        if size(a,2)  == 14 %Emotiv Headset
            Var.capt = 'emotiv';
            Var.ChannelNames = cell(14,1);
            Var.EEGloc = 1:14; Var.EOGloc = [];
        elseif size(a,2)  == 22 %standard EEG cap (2 gusb)
            Var.capt = '22_3EOG_2amp';
            Var.ChannelNames = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8',...
                'P7','P3','Pz','P4','P8','O1','O2'};
            Var.EEGloc = 1:19; Var.EOGloc = 20:22;
        elseif size(a,2)  == 16 % 16 chan gNautilus
            Var.capt = 'gNautilus16';
            Var.ChannelNames = {'Fp1','Fp2','F3','Fz','F4','T7','C3','Cz','C4','T8',...
                'P3','Pz','P4','PO7','PO8','Oz'};
            Var.EEGloc = 1:16; Var.EOGloc = 1;
        elseif size(a,2)  == 8 % 8 chan gNautilus or Signal Generator
            if isfield(c,'DeviceIDs')
                Var.capt = 'gnautilus_8wet';
                Var.EEGloc = 1:8;
                if sum(strcmp(Name,{'T02','T05'}))
                    Var.ChannelNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};
                    Var.EOGloc = 1;
                else
                    disp('DEFAULTING TO NEWER FLEXCAP CONFIG - FZ = CH2')
                    Var.ChannelNames = {'Cz','Fz','P3','Pz','P4','PO7','PO8','Oz'};
                    Var.EOGloc =2;
                end
            else
                Var.capt = 'siggen';
                Var.ChannelNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};
                Var.EEGloc = 1:8; Var.EOGloc = 1;
            end
            
        elseif size(a,2)  == 3
            Var.capt = '3ch';
            Var.ChannelNames = cell(3,1);
        else
            disp('Undefined electrode configuration!')
            return
        end
        
        
        if strcmp(Var.capt,'siggen')
            Art = 0;
        else
        end
        
        disp([num2str(c.PreRunDuration.NumericValue) ' ' num2str(c.PostRunDuration.NumericValue) ...
            ' ' num2str(c.PreSequenceDuration.NumericValue) ' ' num2str(c.PostSequenceDuration.NumericValue)]);
        
        
        
        
        Var.fs{rr} = c.SamplingRate.NumericValue;
        %Trim Data and StimulusCodes
        trimMax = round((c.PreRunDuration.NumericValue-.5)*Var.fs{rr});
        Data = [Data;  a(trimMax:end,:)];
        State.StimulusCode = [State.StimulusCode; b.StimulusCode(trimMax:end)];
        State.StimulusType = [State.StimulusType; b.StimulusType(trimMax:end)];
        State.SourceTime = [State.SourceTime; b.SourceTime(trimMax:end)];
        State.StimulusTime = [State.StimulusTime; b.StimulusTime(trimMax:end)];
        if isfield(c,'ToBeCopied') %audio speller
            TTSpell = [TTSpell c.ToBeCopied.NumericValue];
            TTSpell_cell{kk} = c.ToBeCopied.NumericValue;
            State.SelectedStimulus = [State.SelectedStimulus; b.SelectedStimulus(trimMax:end)];
        else
            TTSpell = [TTSpell cell2mat(c.TextToSpell.Value)];
            TTSpell_cell{kk} = cell2mat(c.TextToSpell.Value);
        end
        if isfield(b,'EyetrackerLeftEyeGazeX')
            State.GazeX = [State.GazeX; b.EyetrackerLeftEyeGazeX(trimMax:end)];
            State.GazeY = [State.GazeY; b.EyetrackerLeftEyeGazeY(trimMax:end)];
        else
            State.GazeX = [State.GazeX; NaN(size(a(trimMax:end,:),1),1)];
            State.GazeY = [State.GazeY; NaN(size(a(trimMax:end,:),1),1)];
        end
        State.SequencePhase = [State.SequencePhase; b.PhaseInSequence(trimMax:end)];
        NewRunName = strcat(NewRunName,Runs{rr});
        
        if isfield(c,'ToBeCopied') %audio speller
            Var.Speller = 4;
        else
            %Checkerboard Speller
            if max(State.StimulusCode) > (c.NumMatrixRows.NumericValue+c.NumMatrixColumns.NumericValue)
                fid = fopen(['data\' Name Session '\'...
                    Name 'S' Session 'R' num2str(str2num(Runs{rr})) '+_mylogfile.txt'],'r');
                if fid==-1
                    fid = fopen(['P:\ALS Proj Data\' Name  '\'...
                        Name 'S' Session 'R' num2str(str2num(Runs{rr})) '+_mylogfile.txt'],'r');
                end
                
                itt=textscan(fid,'%s','delimiter','\n');
                InputText = [InputText; itt{1}];
                Var.Speller = 1; %multiflash grid
                fclose(fid);
            end
            %RC Speller
            if max(State.StimulusCode) == (c.NumMatrixRows.NumericValue+c.NumMatrixColumns.NumericValue)
                Var.Speller = 3;
            end
            %Also check for an RSVP speller, which can have any number of stimulus
            %codes, but will always have NumMatrixRow = 1. max(StimulusCode) will
            %equal 2*NumColumns
            if c.NumMatrixRows.NumericValue == 1 %This is the RSVP speller
                Var.Speller = 2;
            end
        end
        Speller = Var.Speller;
        if ~exist('Speller','var')
            disp('Speller format not detected'); return;
        end
        
        %Associate Target Codes with Targets
        if Var.Speller == 4
            tdtemp = c.Stimuli.Value(1,:);
            TargSymb(kk,:) = tdtemp
        else
            tdtemp = c.TargetDefinitions.Value;
            if size(tdtemp,2)==1 %Multi Menu
                TargSymb(kk,:) = tdtemp{1}(:,2)'; %-- NOTT CORRECT YET
            else
                TargSymb(kk,:) = tdtemp(:,2)';
            end
            SymbSymb = {'Angry%20','Bag%20','Bed%20','Bored%20','Bathroom%20','Carer%20','Cdrink%20',...
                'Clothes%20','Cold%20','Doctor%20','Family%20','Food%20','Glasses%20','Hdrink%20',...
                'Hear%20','Help%20','Hot%20','Spouse%20','Idk%20','Light%20','Medicine%20','Yes%20',...
                'No%20','Nurse%20','Pain%20','Pill%20','Slippers%20','Telephone%20','Tired%20',...
                'Toilet%20','Tv%20','Wheelchair%20','Walker%20','Worried%20','Watch '};
        end
        
        %Number of channels
        Var.NumChans{kk} = size(a,2);
        
        %Number of sequences per trial.  Two times this number is the number of
        %times each stimulus is flashed each trial.
        Var.NumSequences{kk} = c.NumberOfSequences.NumericValue;
        
        %Number of StimulusCodes
        Var.NumStimCodes{kk} = max(b.StimulusCode);
        
        %Duration of the Stimulus State -- This is how long the stimulus is on
        Var.StateDuration{kk} = c.StimulusDuration.NumericValue*Var.fs{kk}/1000;
        
        %Number of trials
        Var.NumTrials{kk} = fix((sum(b.StimulusCode>=1)/Var.StateDuration{kk})/(Var.NumSequences{kk}*Var.NumStimCodes{kk}));
        
        %Assign each trial to a run number
        Var.AllTrials = [Var.AllTrials; kk*ones(Var.NumTrials{kk},1)];
        
        %Load Spatial Filter
        Var.SpatFiltUsed{kk} = c.SpatialFilter.NumericValue;
        
        %Load Classifier
        Var.ClassifierUsed{kk} = c.Classifier.NumericValue;
        Var.ClassifierUsed{kk}(:,4) = cellfun(@str2num,c.Classifier.Value(:,4));
        
        Var.FreeSpell{kk} = sum(b.StimulusType)==0&sum(b.StimulusCode)~=0;
        
        if isfield(b,'EyetrackerLeftEyeGazeX')
            Var.AppPos = [Var.AppPos; repmat([c.WindowLeft.NumericValue c.WindowTop.NumericValue ...
                c.WindowWidth.NumericValue c.WindowHeight.NumericValue],Var.NumTrials{kk},1)];
        end
        
        kk = kk+1;
    end
end
    Data = Data';

%% Check Variables

%Is sampling rate consistant across runs?
if sum(diff(cell2mat(Var.fs))) ~= 0
    disp('Different number of channels in each run!');
end
Var.fs = Var.fs{1};

%Is the feedback type all copy spelling ?
if sum(diff(cell2mat(Var.FreeSpell))) ~= 0
    disp('Different number of channels in each run!');
end
Var.FreeSpell = Var.FreeSpell{1};
if Var.FreeSpell==0
    disp('Run(s) are not copy spelling')
end

%Are channels consistant across runs?
if sum(diff(cell2mat(Var.NumChans))) ~= 0
    disp('Different number of channels in each run!');
end
Var.NumChans = Var.NumChans{1};

%Is the max stimulus code consistant?
if sum(diff(cell2mat(Var.NumStimCodes))) ~= 0
    disp('Different max Stim Code in each run!');
end
Var.NumStimCodes = Var.NumStimCodes{1};

%Is the state duration consistant?
if sum(diff(cell2mat(Var.StateDuration))) ~= 0
    disp('Different State Duration in each run!');
end
Var.StateDuration = Var.StateDuration{1};

%Are the spatial filters the same same size?
try %Are classifiers the same size
    s2mC = cell2mat(Var.SpatFiltUsed);
    
    if sum(sum(diff(reshape(s2mC,size(s2mC,1),size(s2mC,2)/rr,rr),[],3))) ~= 0
        disp('Spatial Filters have different values in each run!');
    end
    Var.SpatFiltUsed = Var.SpatFiltUsed{1};
catch err
    disp('Spatial Filters  are different sizes in each run!');
end
%If the same size, do they contain the same values?


try %Are classifiers the same size
    c2mC = cell2mat(Var.ClassifierUsed);
    %If the same size, do they contain the same values?
    if sum(sum(diff(reshape(c2mC,size(c2mC,1),size(c2mC,2)/rr,rr),[],3))) ~= 0
        disp('Classifiers have different values in each run!');
    end
    Var.ClassifierUsed = Var.ClassifierUsed{1};
catch err
    disp('Classifiers are different sizes in each run!');
end

% %Are number of sequences the same
% if sum(cell2mat(Var.NumSequences))/rr ~= Var.NumSequences{1}
%     disp('Not the same number of sequences in each run, code will not run correctly');
% end
% Var.NumSequences = Var.NumSequences{1};

%Number of rows, columns, stimuli (these are not yet checked for
%consistancy across multiple runs)
if Var.Speller==4
    Var.NR = 1; 
    Var.NC = 1;
    Var.AudioStim = c.Stimuli.Value;
    Var.NumStim = size(c.Stimuli.NumericValue,2);
else
    Var.NR = c.NumMatrixRows.NumericValue;
    Var.NC = c.NumMatrixColumns.NumericValue;
end

if size(tdtemp,2)==1 %Multi Menu
    Var.NS = size(tdtemp{1},1);
else
    Var.NS = size(tdtemp,1);
end

clear kk
%Check state duration again
for scc = 1:Var.NumStimCodes
    sdd = find(State.StimulusCode == scc);
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

if Var.StateDuration ~= StateDuration2
    disp('We have a problem'); return
end

% %If using the CNEamp - remove the first 2 channels
% if isfield(c,'ComPort')==1
%     Data = Data(3:end,:);
%     cneamp=1;
% end

%disp(['Check Vars ' num2str(toc)]);
%% Check Timing of Data
%Check for Dropped Samples (in the case of BCI2000, the data is repeated in
%the next block if the processing exceeds the roundtrip time)
Timerr = diff(State.SourceTime(1:c.SampleBlockSize.NumericValue:end));
StTimerr = diff(State.StimulusTime(1:c.SampleBlockSize.NumericValue:end));
thTime = c.SampleBlockSize.NumericValue/Var.fs*1000;

BadTime = Timerr<thTime-5 | Timerr>thTime+5;
if Var.Figures_On == 1
    figure('Name','Stimulus time and Source time and Stimulus Code')
    plot(StTimerr); hold on;
    plot(Timerr,'r'); ylim([-12 250]);
    plot(State.StimulusCode(1:c.SampleBlockSize.NumericValue:end),'g');
    plot(-BadTime*10,'k')
end
%figure
%plot(Data(:,1)*1000); hold on; plot(b.SourceTime,'r');


end