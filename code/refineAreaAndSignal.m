mkdirIfNotExist('zStackedThreshCorrectedRefined')
 


for n=1:length(longTraces)
        
    %n = 3
    idx = longTraces(n);        

    clf;
    disp(100*n/length(longTraces))
   
    w=2*s+1;
    nR = 6;
    tot = N;

    d = round(tot/nR);
    out = zeros(d*w,nR*w);
    
    %i = round(traj{idx}(:,1));
    %j = round(traj{idx}(:,2));
    
    for k=1:N
        
        disp(k);
        
        %%
        if( ~exist(['zStackedThreshCorrectedRefined/' num2str(k) '.png'],'file') )
            m = imread(['zStackedThreshCorrected/' num2str(k) '.png']);
            imwrite(m<-1,['zStackedThreshCorrectedRefined/' num2str(k) '.png']);
        end

        
        if( ind(idx,k)~=0 )
            
            a = imread(['zStackedYFP/' num2str(k) '.png']);

            data = {};
            for j=1:expe.numberOfColors            
                data{j} = imread( ['img/' getImageName(expe.colorNames{j},k)] );
            end
                        
            m = imread(['zStackedThreshCorrected/' num2str(k) '.png']);
            mrefined = imread(['zStackedThreshCorrectedRefined/' num2str(k) '.png']);
            
            name = ['Measures/' num2str(k) '.mat'];
                
            pos = Me{k}(ind(idx,k)).Centroid;
            pos = round(pos);

            i = pos(1);
            j = pos(2);

            %[seli selj] = getNeiInd( i,j,s,N1,N2 );
            
            seli = (i-s):(i+s); 
            selj = (j-s):(j+s); 
            
            valid = seli > 0 & selj >0  & seli <= N2 & selj <= N1;

                
            seli = seli(valid);
            selj = selj(valid);

            sub_a = a(selj,seli);
            
            sub_data = {};
            for j=1:expe.numberOfColors  
                sub_data{j} = data{j}(selj,seli);
            end
                        
            bgkMask = m(selj,seli);
            bgkMask = imdilate(bgkMask,strel('disk',bgkSize));   
            
            if(doDrawBkg)
               clf 
               imagesc( imnorm(sub_a).*double(~bgkMask) )
               drawnow                
            end

            for j=1:expe.numberOfColors  

                if(~isempty(sub_a(~bgkMask)))
                    bkg(idx,k,j) = median( sub_data{j}(~bgkMask) );
                else
                    bkg(idx,k,j) = -1;
                end

            end
            
            px = Me{k}; %erase other objects
            px=px(ind(idx,k)).PixelIdxList; 
            m=double(m);
            m(px)=2;
            m = m==2;

            m = m(selj,seli);      
            %imagesc(m)
            %drawnow
            
            %check if object touch the border of the frame
            updateArea = 1;
            if( (sum(m(:,1)) + sum(m(:,end)) + sum(m(1,:)) + sum(m(end,:)))>0 )
           
                touchBorder(idx,k) = 1;                
                updateArea = 0;
            end
            
            m = imdilate(m,strel('disk',3));   
    
            if( superSampling > 1)
                sub_a = imresize(sub_a,superSampling);
                
                for j=1:expe.numberOfColors  
                    sub_data{j} = imresize(sub_data{j},superSampling);
                end 
                m = imresize(m,superSampling);
            end
           
            if(updateArea)
                seg = region_seg(sub_a, m, NIteration,0.6,doDraw && mod(k,1)==0); %-- Run segmentation
            else
                seg = m;   
            end
            
            if( sum(seg(:)) == 0 )
                seg = m;
            end
            
            if( dilateSizeAfterRefine >= 1 )
                seg = imdilate(seg,strel('disk',dilateSizeAfterRefine));
            end
            
            mrefined(selj,seli) = mrefined(selj,seli)  +  imresize(seg,1./superSampling);
            imwrite(mrefined,['zStackedThreshCorrectedRefined/' num2str(k) '.png']);

            %
            tmp = imclearborder(seg);

            if( sum(tmp(:)) > 0 )
                 seg = tmp;
            end

            %quantify different colors            
            for j=1:expe.numberOfColors  
                
                tmp = double(sub_data{j}(seg==1));            
                %qu  = quantile(tmp,0.99);    ql = quantile(tmp,0.01);
                %tmp = tmp(  (tmp < qu) & (tmp > ql) );

                refinedMean(idx,k,j)    = mean(tmp);
                refinedSum(idx,k,j)     = sum(tmp);
                refinedStd(idx,k,j)     = std(tmp); 
            end
            
            refinedArea(idx,k) = sum(seg(:))/superSampling^2;
            
            
            %make small image while we are at it   
            b1 = bwmorph(seg,'remove');
            sub_a(b1) = min(sub_a(:));
 
            %save image as png for GUI
            imwrite(sub_a,['snapShots/' num2str(idx) '_' num2str(k) '.png']);

            
            
        end      
    end
    
    
    clf;
    subplot_tight(2,1,1,0.1); hold on
        plot(areaMatrix(idx,:),'k');
        plot(refinedArea(idx,:),'r');
        plot(10*touchBorder(idx,:),'g');
        title('area')
    subplot_tight(2,1,2,0.1); hold on
        %plot(refinedMean(idx,:)-refinedStd(idx,:),'r--');
        %plot(refinedMean(idx,:)+refinedStd(idx,:),'r--');

        plot(refinedMean(idx,:,1),'r');        
        plot(signal(idx,:,1),'k')
        plot(bkg(idx,:,1),'r--')
        
        if(size(signal,3)>1)
            plot(signal(idx,:,2),'k')
            plot(refinedMean(idx,:,2),'g');
            plot(bkg(idx,:,2),'g--')
        end
        
        title('Signal and background')
                
    drawnow
    
end