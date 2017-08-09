function [RefFilt] = Rereference(Data,Var,RerefVal)

%% Rereference Data
switch RerefVal
    case 0
        RefFilt = eye(length(Var.EEGloc));
        %         DataF = Data*SpatFilt;
    case 1
        RefFilt = (-1/length(Var.EEGloc))*(ones(length(Var.EEGloc))-eye(length(Var.EEGloc)));
        RefFilt = RefFilt + eye(length(Var.EEGloc));
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
        Var.NumChans = size(RefFilt,1);
    case 3
        if size(Data,1)==16
            RefFilt = [0 0 0 0 -1 0 0 0 0 0 0 0 1 0;
                0 0 0 0 -1 0 0 0 0 0 0 0 0 1];
        elseif size(Data,1)==22
            RefFilt = [0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 1 0;
                0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1];
        end
        NumChans_old = NumChans;
        Var.NumChans = 2;
    case 4 %reference to C3
        if size(Data,1)==16
        elseif size(Data,1)==22
            RefFilt = eye(length(Var.EEGloc));
            tmp = find(strcmp(c.ChannelNames.Value,'C3'));
            RefFilt(:,tmp) = RefFilt(:,tmp)-1;
        end

end
%disp(['Reref ' num2str(toc)]);
end