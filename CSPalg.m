function [W] = CSPalg(Data,State,Var,RefFilt,ArtWeights2,DoCSP)


if DoCSP == 1
    
    Data_C = (RefFilt*ArtWeights2*Data)';
    [W] = CSP_P300(Data_C,State.StimulusCodeUNArt,State.StimulusType,...
        Var.NumChans, Var.NumTrials, Var.NumStimCodes,Var.range,Var.StateDuration);
    [W2] = CSP_P300_ver2(Data_C,State.StimulusCodeUNArt,State.StimulusType,...
        Var.NumChans, Var.NumTrials, Var.NumStimCodes,Var.range,Var.StateDuration);
    disp(['Done with CSP ' num2str(toc)]);
    
elseif DoCSP == 0
    
    W = eye(size(RefFilt,1));
    W2 = eye(size(RefFilt,1));
    
    disp('CSP not done');
end

