function [traj signal] = getTrajFromInd(ind,tracks,Me)

traj =  cell(1,size(ind,1));

Nsignal =  length(Me{1}(1).MeanIntensity);

signal = zeros([size(ind) Nsignal]);%cell(1,size(ind,1));

for i=1:size(ind,1)
    
    tmp = zeros(length(tracks),2);
    tmpSignal = zeros(length(tracks),1);

    for f=1:length(tracks)

            p = tracks(f).carth;
            
            s = Me{f};
            
            if(ind(i,f) >0)

            tmp(f,1) = p(ind(i,f),1);
            tmp(f,2) = p(ind(i,f),2);
            
            for k=1:Nsignal
            signal(i,f,k) = s(ind(i,f)).MeanIntensity(k); 
            end
                        
            %signal(i,f) = quantile( double( [s(ind(i,f)).PixelValues ]), 0.75);
            
            
            %tmpSignal(f) = s(ind(i,f)).MeanIntensity;
            
            end

    end
    
    traj{i} = tmp;
    %signal{i} = tmpSignal;

end