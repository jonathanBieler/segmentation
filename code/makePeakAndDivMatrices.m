
peakMatrix = zeros(size(signal)); 


for k=1:size(peakMatrix,3)

    for i=1:length(longTraces)

        idx = longTraces(i);
        tmp = signal(idx,:,k);

        maxp = getPeaks(tmp,expe,peakMethod,doDraw);

        for j=1:size(maxp,1)
            peakMatrix(idx,maxp(j,1),k)=1;
        end

    end

    %remove peaks that are too close in time    
    for i=1:length(longTraces)

        idx = longTraces(i);
        tmp = signal(idx,:,k);

        pp = find(peakMatrix(idx,:,k));     

        for j=1:length(pp)-1

            p2pTime = expe.t(pp(j+1))-expe.t(pp(j));        
            if( p2pTime < minTimeBetweenPeaks)           

                    %keep the best one
                    score1 = tmp(pp(j));                
                    score2 = tmp(pp(j+1));

                    if(score1 >= score2)                   
                        peakMatrix(idx,pp(j+1),k) = 0;
                    else
                        peakMatrix(idx,pp(j),k) = 0;
                    end

            end
        end
    end


    %peaks at the begining are almost always false

    peakMatrix(:,1:2,k) = 0;
    %peakMatrix = 0*peakMatrix;

    

end


imagesc(peakMatrix(longTraces,:,1) )
save peakMatrix.mat peakMatrix

%% build division matrix

minTimeBetweenDivisions = 4; 

divMatrix = zeros(size(ind)); 

k = 1;

    for i=1:length(longTraces)

        idx = longTraces(i);
        tmp = signal(idx,:,k);

        %tracking division
        dd = find([divisions.motherInd] == idx);
        for j=1:length(dd)

            ff = divisions(dd(j)).motherFrame;        
            divMatrix(idx,ff,k)=1;

        end

        dd = find([divisions.sisterInd] == idx);    
        for j=1:length(dd)

            ff = divisions(dd(j)).sisterFrame;        
            divMatrix(idx,ff,k)=1;        
        end

        %detect from trace with low false neg (hopefully)
        maxp = getDivisionsFromTrace(tmp,1,0.1);
        for j=1:size(maxp,1)

            %check if there is no division already close by          
            sel = (maxp(j,1)-4):(maxp(j,1)+4);
            sel = sel(sel>0); sel = sel(sel <= size(ind,2));

            if( sum( divMatrix(idx,sel,k)) == 0) 
                divMatrix(idx,maxp(j,1),k)=1;
            end

        end

        %correct divisions that are too close in time    
        divs = find(divMatrix(idx,:,k));    
        for j=1:length(divs)-1

            div2divTime = expe.t(divs(j+1))-expe.t(divs(j));        
            if( div2divTime < minTimeBetweenDivisions)           

                    tmp = fillTrace(tmp);
                    tmp = imnorm(tmp);

                    f = [1 1 -5 1 1];
                    pe = conv(tmp,f,'same');
                    pe(pe<0)=0;
                    pe = smooth(smooth(pe,3)); %%allow for small time differences

                    score1 = pe(divs(j));                
                    score2 = pe(divs(j+1));

                    if(score1 >= score2)                   
                        divMatrix(idx,divs(j+1),k) = 0;
                    else
                        divMatrix(idx,divs(j),k) = 0;
                    end

            end
        end
    end


clf
imagesc(divMatrix(longTraces,:,1))
z = (1-linspace(0,1,100));
cmap = [z; z; z]';
colormap(cmap)

save divMatrix.mat divMatrix

%%

areaMatrix = zeros(size(ind));

for t=1:N

    Measurements = Me{t};

    for ii=1:size(ind,1)
        if( ind(ii,t) > 0 && ind(ii,t) < length(Measurements) )

            areaMatrix(ii,t) = Measurements(ind(ii,t)).Area;
        end
    end
end

imagesc(areaMatrix)

save areaMatrix.mat areaMatrix

