function [DD, Var] = ApplySpatialFilter(Data,Var,ArtWeights2,RefFilt,W,NewRunName,Name,Sessions,DoCSP)

    %     ChansToKeep = input('Which Channels to keep?'
    if DoCSP == 1
        if RerefVal == 1 %This is because the
            %last CSP channel unser CAR rereferencing is zeros.
            if Var.NumChans <= 4
                ChansToKeep = 1:Var.NumChans-1;
            else
                ChansToKeep = [1 2 Var.NumChans-2 Var.NumChans-1];
            end
        else
            if Var.NumChans <= 3
                ChansToKeep = 1:Var.NumChans;
            else
                ChansToKeep = [1 2 Var.NumChans-1 Var.NumChans];
            end
        end
        
    else
        ChansToKeep = 1:size(W,1);
    end
    
    Var.NumChans = length(ChansToKeep);
    W = W(ChansToKeep,:);
    
    
    Data(isnan(Data))=0;
    SpatFilt = W*RefFilt*ArtWeights2;
    %     SpatFilt2 = W2*RefFilt*ArtWeights2;
    %Combine Spatial Filter and Artifact Matrix into Spatial Filter for online
    %use.  Here i am first doing artifact rejection, then spatial
    %filtering.  There is a slight difference with the order.
    
    %Save the Spatial Filter
    %     dlmwrite(['C:\Documents and Settings\amg5106\My Documents\'...
    %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' Name Session '\'...
    %         Name 'S' Session 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
    try
        dlmwrite(['data\' Name Sessions{Var.max_sess_i} '\'...
            Name 'S' Sessions{Var.max_sess_i} 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
    catch
        dlmwrite(['P:\ALS Proj Data\' Name Sessions{Var.max_sess_i} '\'...
            Name 'S' Sessions{Var.max_sess_i} 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
    end
    DD = SpatFilt*Data;
    
    
    if isempty(Var.EOGloc);
        Var.EOGDD = NaN(0,size(Data,2));
    else
        Var.EOGDD = Data(Var.EOGloc,:);
    end
