function [nObj] = doMeasures(N,threshFolder,expe)


mkdirIfNotExist('Measures')

nObj=zeros(1,N);

for k=1:N
    
    disp(100*k/N);
    
    for i=1:expe.numberOfColors
        a{i} = imread(['img/' getImageName(expe.colorNames{i},k)]);         
    end
    
    b = imread([threshFolder num2str(k) '.png']);

    labeledImage = bwlabel(b, 8);     % Label each blob so we can make measurements of it

    Measurements  = regionprops(labeledImage, a{1}, {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});   
    
    %
                        
    for j=1:expe.numberOfColors        
        tmpMeasurements{j}  = regionprops(labeledImage, a{j}, {'PixelValues','MeanIntensity'});           
    end
    
    for i=1:length(Measurements)
    
        newPixelValues = zeros([expe.numberOfColors size(Measurements(i).PixelValues)]);
        newMeanIntensity = zeros([expe.numberOfColors size(Measurements(i).MeanIntensity)]);
    
        for j=1:expe.numberOfColors 
            newPixelValues(j,:) = tmpMeasurements{j}(i).PixelValues;
            newMeanIntensity(j,:) = tmpMeasurements{j}(i).MeanIntensity;
        end    
        
        Measurements(i).PixelValues = newPixelValues;
        Measurements(i).MeanIntensity = newMeanIntensity;
    end
    
    %
    name = ['Measures/' num2str(k) '.mat'];
            
    save(name,'Measurements')
        
    
    nObj(k)=length(Measurements);
end

