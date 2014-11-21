function [] = removeCytoplasm(imInd,threshFolder,ratioThresh)

%ratioThresh = 1.2; %inside signal must be at least ratioThresh*outside signal

a = imread(['zStackedYFP/' num2str(imInd) '.png']);
b = imread(['zStackedThresh/' num2str(imInd) '.png'] );

labeledImage = bwlabel(b, 8);     % Label each blob so we can make measurements of it
Measurements = regionprops(labeledImage, a, {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});   


for i=1:length(Measurements)

    m=b;
    
    %erase all objects except the current one   
    px = Measurements(i).PixelIdxList; 
    m=double(m);
    m(px)=2;
    m = m==2;
        
    %measure the signal inside the object
    signalIn = mean( a(m==1) ); 
    
    %measure the signal outside the object
    se = strel('disk',6);
    m = m-imdilate(m,se);
    m = m==-1;
    
    signalOut = mean( a(m==1) );
        
    ratio = signalIn/signalOut;
    
    %probably cytoplasm: erase the object from the image        
    if(ratio < ratioThresh)
        b(px) = 0;        
    end
    
           
end

imwrite(b,['zStackedThresh/' num2str(imInd) '.png'] );