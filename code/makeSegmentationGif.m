
for k=1:N

    
        disp(100*k/N)
        
         
        seg = [threshFolder num2str(k) '.png'];
        a = imread(['zStackedYFP/' num2str(k) '.png']);
        %a=double(a);
        %a=a.^0.8;
        b = imread(seg);
        
        %superpose both
        b = (edge(b,'log',0));
        b = imdilate(b,strel('disk',1));
        
        a(b==1) = 250;
        
        
        %a=double(a)
                              
        [a ] =  ind2rgb((a)/1.5,jet);    
        [a,cm] = rgb2ind(a,256);
        
        if(doDraw)
            imagesc(a)
            colormap(cm)
            drawnow
        end
    
        if k == 1;
            imwrite(a,cm,['segmentation.gif'],'gif', 'Loopcount',inf,'DelayTime',0.15);
        else
            imwrite(a,cm,['segmentation.gif'],'gif','WriteMode','append','DelayTime',0.15);
        end

   
end

