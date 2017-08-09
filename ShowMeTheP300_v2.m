
%Displays raw and processed P300 data
%This has been updated from BuildClassifier_teleBCI7 and modularized 8/7/17

function [CCdata, NNdata, Spec, ff, NumSequences, AIc, AIn,  TargWait, GAZE, FaceData] = ...
    ShowMeTheP300_v2(Name,Sessions,Runs,RerefVal,Art,DoCSP,capver)
tic
currdir = cd;
runloc = currdir(1);






DoHB = 0;

Var.Figures_On = 0;
Var.SequencestoRemove = 0;  %For 20 flashes, there are 10 sequences.
Var.EEGloc = 1:8;
Var.EOGloc = 1;


%Load and Join Data
[Data, State, Var, InputText, NewRunName, TTSpell, Art,c]=...
    Load_and_Join_BCI2000_Data(Name,Sessions,Runs,Var,Art);
NumSequences = Var.NumSequences;
TargWait = 0;  %Used to be in ShowMeTheP300, but removed because not sure what it was for.

Var.range = (round(.45*Var.fs):round(1.2*Var.fs)); %Use this data range for classifiation with g.Nautilus
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


%% Heartbeat Filter
[HBData] = HeartbeatFilter(Data,State,Var,DoHB);

%% Artifact Correction
[ArtWeights2, State, Var] = OcularArtifact(HBData,Var,State,Art);



%% Rereference
[RefFilt] = Rereference(HBData,Var,RerefVal);



%% Do CSP
[W] = CSPalg(HBData,State,Var,RefFilt,ArtWeights2,DoCSP);



%% Save Spatial Filter
[DD, Var] = ApplySpatialFilter(Data,Var,ArtWeights2,RefFilt,W,NewRunName,Name,Sessions,DoCSP);
disp(['Done with filtering ' num2str(toc)]);


%% Spectral Analysis
for cc = 1:size(DD,1)
    [Spec(:,cc), ff] = pwelch(DD(cc,:),Var.fs*8,0,Var.fs*8,Var.fs);
    
    
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


if Var.FreeSpell == 1
%     Clocs = find(StimulusCode>0);
%     Clocs = Clocs(1:StateDuration:end);
%     Nlocs = find(StimulusCode>0);
%     Nlocs = Nlocs(1:StateDuration:end);
    CCdata = [];%zeros(length(range),size(DD,1),length(Clocs));
    NNdata = [];%zeros(length(range),size(DD,1),length(Nlocs));
    AIc = []; AIn = []; GAZE = []; FaceData = [];
    return
end



%% Parse Log File
[Var] = ParseLogFile(Data,State,Var,InputText,TTSpell);


%% EyeGaze


if Var.Speller==4 %Can include this later -- for now just ignore eye tracking for audio speller
    GAZE = [];
else
    if sum(~isnan(State.GazeX))~=0 && Art==0
        [GAZE.acc, GAZE.var, GAZE.invalid, GAZE.targetbox] = P300_EyeTracking_inP300Classifier_v4(...
         State,Var,0,TTSpell,c);
        disp(['Finished with eye tracking ' num2str(toc)]);
    else
        disp('No eye tracking data detected');
    end
    
    waitbar(.25,wtbrr);
    eyye.gazeacc = gazeacc;
    eyye.gazevar = gazevar;
    eyye.gazeinvalid = gazeinvalid;
end






if Var.Speller==4 %Can include this later -- for now just ignore eye tracking for audio speller
    GAZE = [];
else
    if ~isempty(State.GazeX) && Art==99  %This should change back to Art ==0
        if runloc == 'C'
        end
        [GAZE.acc, GAZE.var, GAZE.invalid, GAZE.targetbox] = EyeTracking(...
            State.GazeX,State.GazeY, Data(Var.EOGloc,:), Var.SC,Var.SR, State.StimulusCode, ...
            State.SequencePhase, State.StimulusType, Var.StateDuration, ...
            Var.NumSequences+Var.SequencestoRemove, Var.Speller,...
            Var.fs, c, 0, 1);
        cd(currdir);
    else
        GAZE = [];
    end
end
    %disp(['Finished with eye tracking ' num2str(toc)]);

%% Organize Transformed Data into Choice/Not Choice

if Var.FreeSpell == 0
    Clocs = find(State.StimulusCode>0 & State.StimulusType == 1);
    Clocs = Clocs(1:Var.StateDuration:end);
    Nlocs = find(State.StimulusCode>0 & State.StimulusType == 0);
    Nlocs = Nlocs(1:Var.StateDuration:end);
    
    %Check if Nlocs or Clocs will go out of range (run finished early)
    outofit = Nlocs+max(Var.range)<size(DD,2);
    Nlocs = Nlocs(outofit);
    outofit = Clocs+max(Var.range)<size(DD,2);
    Clocs = Clocs(outofit);
    
    
    CCdata = zeros(length(Var.range),size(DD,1)+size(Var.EOGDD,1),length(Clocs));
    NNdata = zeros(length(Var.range),size(DD,1)+size(Var.EOGDD,1),length(Nlocs));
    AIc = zeros(length(Var.range),length(Clocs));
    AIn = zeros(length(Var.range),length(Nlocs));
    
    
    for cl = 1:length(Clocs)
        CCdata(:,:,cl) = cat(2,DD(:,Clocs(cl)+Var.range)', Var.EOGDD(:,Clocs(cl)+Var.range)');
        AIc(:,cl) = Var.AI(Clocs(cl)+Var.range);
    end
    for nl = 1:length(Nlocs)
        NNdata(:,:,nl) = cat(2,DD(:,Nlocs(nl)+Var.range)', Var.EOGDD(:,Nlocs(nl)+Var.range)');
        AIn(:,nl) = Var.AI(Nlocs(nl)+Var.range);
    end
    
    
    trlrem = squeeze(sum(sum(CCdata>100 | CCdata<-100)))>0;
    CCdata(:,:,trlrem) = [];
    trlrem = squeeze(sum(sum(NNdata>100 | NNdata<-100)))>0;
    NNdata(:,:,trlrem) = [];
else
end

%Not working --- need to edit the parselog in order to get FaceOF var.

% FaceData.C = cell(1,24);
% FaceData.N = cell(1,24);
% if Var.Speller == 1 && (Art==0 || Art==2)
%     for i = 1:sum(cell2mat(Var.NumTrials))
%         for j = 1:24 %Loop through face stimuli
%             Faceloc = find(str2num(FaceOF{i})==j);
%             for k = Faceloc
%                 tmpSC = str2num(cell2mat((SC{i}(k)))); tmpSR = str2num(cell2mat(SR{i}(k)));
%                 TargC = find(State.StimulusCode == tmpSC); TargC = TargC(1:StateDuration:end); TargC = TargC((i-1)*NumSequences+1:i*NumSequences);
%                 TargR = find(State.StimulusCode == tmpSR); TargR = TargR(1:StateDuration:end); TargR = TargR((i-1)*NumSequences+1:i*NumSequences);
%                 tl = c.TextToSpell.Value{:}(i);
%                 C_ind = find(strcmp(tl,TargSymb));
%                 for l = 1:NumSequences
%                     
%                     if k == C_ind
%                         FaceData.C{j} = cat(3,FaceData.C{j},DD(:,TargC(l)+range(1)-1:TargC(l)+range(end)-1));
%                         FaceData.C{j} = cat(3,FaceData.C{j},DD(:,TargR(l)+range(1)-1:TargR(l)+range(end)-1));
%                     else
%                         FaceData.N{j} = cat(3,FaceData.N{j},DD(:,TargC(l)+range(1)-1:TargC(l)+range(end)-1));
%                         FaceData.N{j} = cat(3,FaceData.N{j},DD(:,TargR(l)+range(1)-1:TargR(l)+range(end)-1));
%                     end
%                 end
%             end
%         end
%     end
% else
% end
FaceData.C = [];
FaceData.N = [];


fprintf(['Plot ' num2str(toc) '\n\n']);

end
