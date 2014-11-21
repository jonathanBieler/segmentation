for n=1:length(longTraces)
    
    idx = longTraces(n);
    disp(100*n / length(longTraces));
    
    i = round(traj{idx}(:,1));
    j = round(traj{idx}(:,2));

    clf;

    s = 50;
    w = 2*s+1;

    %find fist non zero index
    nonZeroIdx = find(ind(idx,:) > 0);
    nonZeroIdx = nonZeroIdx(1);

    timeInd = round(linspace(1,w,length(1:NToTrack)));
    
    %make some frames before the trace begin
    for k= max(nonZeroIdx-8,1):nonZeroIdx
        
        
        a = imread([inputFolder num2str(k) '.png']);
        m = imread(['zStackedThreshCorrected/' num2str(k) '.png']);
        
        [seli selj] = getNeiInd(i(nonZeroIdx),j(nonZeroIdx),s,N1,N2);

        sub_a = a(selj,seli);

  
        %save image as png for GUI
        imwrite(sub_a,['snapShots/' num2str(idx) '_' num2str(k) '.png']);

        sub_a =  ind2rgb((sub_a)/2.0,hot);
        [imind,cm] = rgb2ind(sub_a,256);
        
        
    end
    
    for k=nonZeroIdx:NToTrack
               
        a = imread([inputFolder num2str(k) '.png']);
        m = imread(['zStackedThreshCorrected/' num2str(k) '.png']);
        
        [seli selj] = getNeiInd(i(k),j(k),s,N1,N2);

        if( ind(idx,k)~=0 )
         
            sub_a = a(selj,seli);
            
            if(ind(idx,k)~=0)
                px = Me{k}; %erase other objects
                px=px(ind(idx,k)).PixelIdxList; 
                m=double(m);
                m(px)=2;
                m = m==2;
            end
            
            m = m(selj,seli);                        
         
            %make a nice gif            
            b1 = bwmorph(m,'remove');
            sub_a(b1) = min(sub_a(:));
 
            %save image as png for GUI
            imwrite(sub_a,['snapShots/' num2str(idx) '_' num2str(k) '.png']);

            
            if(doDraw)
                imagesc(sub_a); drawnow;
                %pause(0.1)
            end
                    
        end
    end
end