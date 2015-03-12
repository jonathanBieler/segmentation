[N1,N2] = getImageDimensions(expe);

for n=traces
    
    idx = longTraces(n);
    disp(100*n / length(longTraces));
    
    if(strcmp(threshFoler,'zStackedThreshCorrectedRefined'))

        i = floor(refinedTrajX(idx,:));
        j = floor(refinedTrajY(idx,:));
        
        ifull = floor(binsize*(refinedTrajX(idx,:)));
        jfull = floor(binsize*(refinedTrajY(idx,:)));
        
    else
        i = floor(trajX(idx,:));
        j = floor(trajY(idx,:));
        
        ifull = floor(binsize*(trajX(idx,:)));
        jfull = floor(binsize*(trajY(idx,:)));
    end
  
   
    w = 2*s+1;

    %find fist non zero index
    nonZeroIdx = find(ind(idx,:) > 0);
    nonZeroIdx = nonZeroIdx(1);

    timeInd = floor(linspace(1,w,length(1:NToTrack)));
    
    %make some frames before the trace begin
    for k= max(nonZeroIdx-8,1):nonZeroIdx
        
        
        %a = imread([inputFolder num2str(k) '.png']);
        
        if( useFullSizeImages )
            a = imread(['fullSizeimg/' getImageName(expe.colorNames{colorIndex},k)]);
        else
            a = imread(['img/' getImageName(expe.colorNames{colorIndex},k)]);
        end
        
        m = imread([threshFoler '/'  num2str(k) '.png']);
        
        [seli selj] = getNeiInd(i(nonZeroIdx),j(nonZeroIdx),s,N1,N2);

        if( useFullSizeImages)
            %[seli selj] = getNeiInd(1+binsize*(i(nonZeroIdx)-1),1+binsize*(j(nonZeroIdx)-1),binsize*s,binsize*N1,binsize*N2);
            [seli selj] = getNeiInd(ifull(nonZeroIdx),jfull(nonZeroIdx),binsize*s,binsize*N1,binsize*N2);
        end
                
        sub_a = a(selj,seli);

        %save image as png for GUI
        imwrite(sub_a,['snapShots/' num2str(idx) '_' num2str(k) '.png']);

        
    end
    
    %%
    
    for k=nonZeroIdx:N
               
        %a = imread([inputFolder num2str(k) '.png']);
                
        if( useFullSizeImages )
            a = imread(['fullSizeimg/' getImageName(expe.colorNames{colorIndex},k)]);
        else
            a = imread(['img/' getImageName(expe.colorNames{colorIndex},k)]);
        end
                
        m = imread([threshFoler '/'  num2str(k) '.png']);
        
        [seli selj] = getNeiInd(i(k),j(k),s,N1,N2);
        
        if( ind(idx,k)~=0 )
         
            if( useFullSizeImages)
                %[selif seljf] = getNeiInd(1+binsize*(i(k)),1+binsize*(j(k)-1),binsize*s,binsize*N1,binsize*N2);
                [selif seljf] = getNeiInd(ifull(k),jfull(k),binsize*s+1,binsize*N1,binsize*N2);
                sub_a = a(seljf,selif);
            else
                sub_a = a(selj,seli);
            end
                                    
%             if(ind(idx,k)~=0)
%                 px = Me{k}; %erase other objects
%                 px=px(ind(idx,k)).PixelIdxList; 
%                 m=double(m);
%                 m(px)=2;
%                 m = m==2;
%             end

            if( useFullSizeImages )
               m = imresize(m, size(a),'nearest');
               m = m(seljf,selif);  
            else
                
               m = m(selj,seli);   
            end
            
            

         
            %make a nice image            
            b1 = bwmorph(m,'remove');
            sub_a(b1) = min(sub_a(:));
 
            %save image as png for GUI
            imwrite(sub_a,['snapShots/' num2str(idx) '_' num2str(k) '.png']);
            
            if(doDraw)
                clf;
                imagesc(sub_a); drawnow;
                pause(0.1)
            end
                    
        end
    end
end