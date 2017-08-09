function [ArtWeights2, State, Var] = OcularArtifact(Data,Var,State,Art)
%% Artifact Correction
State.StimulusCodeUNArt = State.StimulusCode;
State.StimulusTypeUNArt = State. StimulusType;

BuildClassifier=1;

%Remove dead data
SmInt = [];
for i= 256:256:size(Data,2)-255
    if std(Data(Var.EOGloc(1),i:i+255))<1;
        SmInt  = [SmInt i];
    end
end
if ~isempty(SmInt)
    Srange = fix(-3*Var.fs:5*Var.fs);
    SmInt = repmat(SmInt,length(Srange),1)+repmat(Srange,length(SmInt),1)';
    SmInt = SmInt(:);
    SmInt = SmInt(SmInt>0 & SmInt<size(Data,2));
    SmLoc = zeros(size(Data,2),1);
    SmLoc(SmInt)=1;
    State.StimulusType(logical(SmLoc))=0; %= StimulusType(~ArtLoc);
    State.StimulusCode(logical(SmLoc))=0; %= StimulusCode(~ArtLoc);
else
    SmLoc = [];
end

%Remove superlarge data artifacts 1s before, 5 seconds after
LInt = find(abs(Data(Var.EOGloc(1),:))>1e3);
if ~isempty(LInt)
    Lrange = fix(-1*Var.fs:5*Var.fs);
    LInt = repmat(LInt,length(Lrange),1)+repmat(Lrange,length(LInt),1)';
    LInt = LInt(:);
    LInt = LInt(LInt>0 & LInt<size(Data,2));
    LLoc = zeros(size(Data,2),1);
    LLoc(LInt)=1;
    State.StimulusType(logical(LLoc))=0; %= StimulusType(~ArtLoc);
    State.StimulusCode(logical(LLoc))=0; %= StimulusCode(~ArtLoc);
else
    LLoc = [];
end

%Ocular artifacts
if isempty(Var.EOGloc)
    disp('Cannot reject artifacts because no EOG channels were recorded');
    Art = 0;
end
if Art == 1 %Remove data displaying artifacts
    disp('Removing ocular artifacts....');
    
    ArtInt = find(Data(Var.EOGloc(1),:)>50&circshift(Data(Var.EOGloc(1),:),[0 1])<50);
    ArtInt = [ArtInt find(Data(Var.EOGloc(1),:)<-50&circshift(Data(Var.EOGloc(1),:),[0 1])>-50)];
    
    Arange = fix(-1*Var.fs:.5*Var.fs);  %This was changed 5/29/17 from -.1 to -1
    %because (a) we have to account for the .45 second lag in EEG data
    %to account for the fact that the large peak was corrupting the EEG
    %data for preceding stimulus codes that were only set to zero .1
    %seconds before.
    ArtInt = repmat(ArtInt,length(Arange),1)+repmat(Arange,length(ArtInt),1)';
    %     ArtInt = ArtInt';
    ArtInt = ArtInt(:);
    ArtInt = ArtInt(ArtInt>0 & ArtInt<size(Data,2));
    ArtLoc = zeros(size(Data,2),1);
    ArtLoc(ArtInt)=1;
    Var.AI = ArtLoc;
    State.StimulusType(logical(ArtLoc))=0; %= StimulusType(~ArtLoc);
    State.StimulusCode(logical(ArtLoc))=0; %= StimulusCode(~ArtLoc);
    
    if Var.Figures_On == 1
        figure('Name','EOG channel 1, locations of artifact, and stimulus code')
        plot(Data(Var.EOGloc(1),:)); hold on;
        plot(50*ArtLoc,'r');
        plot(50*SmLoc,'c');
        plot(50*LLoc,'m');
        %     figure
        %     tt(1) = subplot(211); plot(Data(EEGloc,:)');
        %     tt(2) = subplot(212); plot(Data(EOGloc,:)');
        %     linkaxes(tt)
    end
    
    %Find incomplete trials
    sdtemp = find(diff(double(State.StimulusCode))~=0);
    nftemp = diff(sdtemp)<Var.StateDuration;
    nfstart = sdtemp(nftemp)+1;
    nfend = sdtemp(find(nftemp)+1);
    
    for ff = 1:length(nfstart)
        State.StimulusCode(nfstart(ff):nfend(ff)) = 0;
        State.StimulusType(nfstart(ff):nfend(ff)) = 0;
    end
    
    
    %     ismember(diff(double(StimulusCode))>0
    
    
    if Var.Figures_On == 1
        plot(State.StimulusCode,'k','LineWidth',2);
        plot(State.StimulusType,'g','LineWidth',2);
    end
    
    ArtWeights = eye(size(Data,1));
    ArtWeights2 = eye(length(Var.EEGloc),size(Data,1));
    NumChans = length(Var.EEGloc);
    disp(['Done with Art ' num2str(toc)]);
elseif Art == 2 %Artifact regression
    if BuildClassifier == 1
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
            if length(Var.EOGloc)==3
                ArtInt = find(Data(Var.EOGloc(1),:)-Data(Var.EOGloc(2),:)>75 | ...
                    Data(Var.EOGloc(1),:)-Data(Var.EOGloc(2),:)<-75 | ...
                    Data(Var.EOGloc(1),:)-Data(Var.EOGloc(3),:)>75 | ...
                    Data(Var.EOGloc(1),:)-Data(Var.EOGloc(3),:)<-75);
            else
                ArtInt = find(Data(Var.EOGloc(1),:)>75 | ...
                    Data(Var.EOGloc(1),:)<-75);
            end
        catch
            disp('Cannot remove data because no EOG channels were recorded');
            return
        end
        
        Arange = fix(-1*Var.fs:.5*Var.fs);
        
        ArtInt = repmat(ArtInt',1,length(Arange))+repmat(Arange,length(ArtInt),1);
        ArtInt = ArtInt';
        ArtInt = unique(ArtInt(:));
        ArtInt = ArtInt(ArtInt>0 & ArtInt<size(Data,2));
        ArtLoc = logical(zeros(size(Data,2),1));
        ArtLoc(ArtInt) = true;
        Var.AI = ArtLoc;
        if Var.Figures_On == 1
            if length(Var.EOGloc)==3
                figure('Name','Differential EOG and location of artifacts')
                plot(Data(Var.EOGloc(1),:)-Data(Var.EOGloc(2),:))
                hold on
                plot(Data(Var.EOGloc(1),:)-Data(Var.EOGloc(3),:))
                plot(100*ArtLoc,'r')
            else
                figure('Name','EOG and location of artifacts')
                plot(Data(Var.EOGloc(1),:)); hold on;
                plot(100*ArtLoc,'r');
            end
        end
        
        EOGD = Data(:,ArtLoc);
        
        if ~isempty(EOGD)
            %the function covm adds an additional column of ones in front of the data
            %and is necessary for regress_eog.m
            if length(Var.EOGloc)==3
                [R] = regress_eog(covm(EOGD','E'),Var.EEGloc, ...
                    sparse([Var.EOGloc(1),Var.EOGloc(3),Var.EOGloc(2),Var.EOGloc(1)],[1,1,2,2],[1,-1,1,-1]));
            elseif length(Var.EOGloc)==2
                [R] = regress_eog(covm(EOGD','E'),Var.EEGloc, ...
                    sparse([Var.EOGloc(1),Var.EOGloc(2)],[1,1],[1,-1]));
            else
                [R] = regress_eog(covm(EOGD','E'),Var.EEGloc, Var.EOGloc);
            end
            %Create full matrix for online artifact reduction
            %I believe this is the way they say to do it (pad Data with a channel
            %of ones -- this introduces a bias to the output channel) (see DD2 below).
            %However, this padding is not something I want to do online, and since
            %it is only a bias, we can remove the first column of ArtWeights.
            ArtWeights = full(R.r0)';
            ArtWeights2 = ArtWeights(Var.EEGloc,:);
            %         ArtWeights = full(R.r1)';
            %         ArtWeights2 = ArtWeights(:,2:end);
            %DD2 = [ones(size(Data,1),1),Data] * ArtWeights';
        else
            ArtWeights2 = eye(length(Var.EEGloc),size(Data,1));
        end
        
        Var.NumChans = length(Var.EEGloc);
        
        
    elseif BuildClassifier == 0
    end
    disp(['Done with Art ' num2str(toc)]);
else
    ArtWeights = eye(size(Data,1));
    ArtWeights2 = eye(length(Var.EEGloc),size(Data,1));
    Var.NumChans = length(Var.EEGloc);
    Var.AI = zeros(size(Data,2),1);
    %     StimulusCodeUNArt = StimulusCode;
    %     StimulusTypeUNArt = StimulusType;
    disp('No artifacting done');
end
%Using the correction coefficients, transform entire training run to reduce
%artifact
%
if Var.Figures_On == 1
    if BuildClassifier == 1
        Dw = ArtWeights2*Data;
        %
        %
        %Plot the F and O raw and artifacted, as well as EOG channels
        figure('Name','Original and Corrected Ch1 data, Ch2 data, and EOG channels')
        tt(1)=subplot(311); plot(Data(1,:)); hold on; plot(Dw(1,:),'r');
        tt(2)=subplot(312); plot(Data(2,:)); hold on;
        plot(Dw(2,:),'r');
        tt(3)=subplot(313); plot(Data(Var.EOGloc,:)'); hold on
        linkaxes(tt,'x');
        
        clear Dw;
    end
end


end
