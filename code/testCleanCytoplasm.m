moviePath = '/Users/bieler/Desktop/matlab/23september_40deg/movie2/';

imInd = 80;

a = imread([moviePath 'zStackedYFP/' num2str(imInd) '.png']);

inputFolder = 'zStackedYFP/';

segPara.minSize= 1.7; 
segPara.maxSize= 3.0;
segPara.thresh = 0.91;
segPara.openSize = 5.0;

doDraw = 0;

tic
generateFilters;
segmentImage(imInd,segPara,inputFolder,filters,openFilter,doDraw);
toc

%%

ratioThresh = 1.2; %inside signal must be at least ratioThresh*outside signal

a = imread(['zStackedYFP/' num2str(imInd) '.png']);
b = imread(['zStackedThresh/' num2str(imInd) '.png'] );


labeledImage = bwlabel(b, 8);     % Label each blob so we can make measurements of it
Measurements = regionprops(labeledImage, a, {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});   

ratios = [];

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
        
    ratios(i) = signalIn/signalOut;
    
    %probably cytoplasm: erase the object from the image        
    if(ratios(i) < ratioThresh)
        b(px) = 0;        
    end
    
        
end


min(ratios)

%hist(ratios,20)

%

%superpose both
b = (edge(b,'log',0));

a(b==1) = min(a(:));

imagesc(imnorm(a).^0.5)
colormap jet