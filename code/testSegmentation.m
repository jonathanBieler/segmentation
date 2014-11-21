%train


inputFolder = 'zStackedYFP/';

load segPara.mat
load segParaFrames.mat

frames = [5 20];

coeff = trainSegmentation(inputFolder,segPara,segParaFrames,frames)

save ../segmentationCoeff.mat coeff





%%

K = [coeff.const]; K=K(1);

L = [coeff.linear];
L = L(:,1);

seg = (K + features*L) > 0;
seg = reshape(seg,512,512);

seg = imfill(seg,'holes');
seg = imopen(seg, openFilter);

imagesc(seg+trueSeg)

err = ( sum(sum( abs(seg-trueSeg) )))

%%

imagesc(imnorm(b)-0.0*double(seg))
    
%% 

frame = 24;

seg = segmentLDA(inputFolder,segPara,segParaFrames,frame,coeff);

imagesc(seg)

%% Compare with manual corrected


trueSeg = imread(['zStackedThreshCorrected/' n2s(frame) '.png']);

imagesc(seg+trueSeg)

err = ( sum(sum( abs(seg-trueSeg) )))


%% compare with auto segment

trueSeg = imread(['zStackedThresh/' n2s(frame) '.png']);

clf
imagesc(seg+trueSeg)

err = ( sum(sum( abs(seg-trueSeg) )))

%%

clf
imagesc(seg2+trueSeg)
err= ( sum(sum( abs(seg2-trueSeg) )))    

%% filtering objects


nObj=zeros(1,N);

Me = [];
labels = [];
frames = [];
objs = [];

threshFolder = 'zStackedThreshSplit/';

for k=1:N
    %disp(100*k/N);
    
    a = imread(['zStackedYFP_Data/' num2str(k) '_1.png']);
    a2 = imread(['zStackedYFP_Data/' num2str(k) '_2.png']);
    b = imread([threshFolder num2str(k) '.png']);
    bTrue = imread(['zStackedThreshCorrected/' num2str(k) '.png']);

    labeledImage = bwlabel(b, 8);     % Label each blob so we can make measurements of it

    Measurements  = regionprops(labeledImage, a, {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});   
    Measurements2 = regionprops(labeledImage, a2, {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});   
        
    for i=1:size(Measurements)
        Measurements(i).PixelValues2 = Measurements2(i).PixelValues;
        Measurements(i).MeanIntensity2 = Measurements2(i).MeanIntensity;
        
        % check if the object is present on the correct image
        if(sum( bTrue(Measurements(i).PixelIdxList) ) > 10) 
            labels = [labels; 1];    
        else
            labels = [labels; 0];
        end
        
        frames = [ frames; k];
        objs = [ objs; i];
                
    end
    
        
    Me = [Me; [[Measurements.MeanIntensity]/10; [Measurements.Area]; 30*[Measurements.Eccentricity] ]'];
    
    nObj(k)=length(Measurements);
end

clf;
plot(nObj)
ylabel('number of objects')

[class,err,post,logp,coeff]  = classify(Me,Me,labels,'quadratic');

err
sum(class == 0)
sum(labels == 0)

save ../objectFilteringCoef.mat coeff

%% compare

k = 13;
b = imread(['zStackedThreshSplit/' num2str(k) '.png']);
bTrue = imread(['zStackedThreshCorrected/' num2str(k) '.png']);

imagesc(b+bTrue)

%% corrected

b = imread(['zStackedThreshSplit/' num2str(k) '.png']);
bTrue = imread(['zStackedThreshCorrected/' num2str(k) '.png']);


name = ['Measures/' num2str(k) '.mat'];
load(name)
    

for i=1:length(Measurements)
    
    features = [[Measurements(i).MeanIntensity]/10; [Measurements(i).Area]; 30*[Measurements(i).Eccentricity] ]';

    %0 < K + x*L + x*Q*x'
    
    K = [coeff.const]; K=K(1);

    L = [coeff.linear];
    L = L(:,1);
    
    Q = [coeff.quadratic]; Q =  Q(1:3,1:3);

    if( (K + features*L  + features*Q*features' ) < -2 )
        b(Measurements(i).PixelIdxList) = 0;
    end
end


imagesc( b + bTrue )

err = ( sum(sum( abs(b-bTrue) )))


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% segment with prior

frame = 28;

seg = imread(['zStackedThresh/' n2s(frame-1) '.png']);
seg = bwmorph(seg,'shrink',Inf);

P = imfilter(double(seg),fspecial('Gaussian',200,25));
prior = imnorm(P);

seg = imread(['zStackedThresh/' n2s(frame) '.png']);

imagesc(prior + double(seg))
colorbar

%%

b = imread(['zStackedYFP/' num2str(frame) '.png']);

load segPara.mat
load segParaFrames

doDraw = 0;

para.minSize  = interp1(segParaFrames,[segPara(:).minSize  ],frame);
para.maxSize  = interp1(segParaFrames,[segPara(:).maxSize  ],frame);
para.thresh   = interp1(segParaFrames,[segPara(:).thresh   ],frame);
para.openSize = interp1(segParaFrames,[segPara(:).openSize ],frame);
para.ratioThresh = interp1(segParaFrames,[segPara(:).ratioThresh ],frame);

[filters openFilter] = generateFilters(para,doDraw);


Nf=size(filters,3);
c = zeros(size(b,1),size(b,2),Nf+1);

for i=1:Nf
    
   tmp = imfilter(b,filters(:,:,i),'replicate');
   
   tmp(tmp<0)=0;
   
   %c(:,:,i) = imnorm(medfilt2(tmp,[2 2]));
   c(:,:,i) = imnorm(tmp);
   

   
    %imagesc(c(:,:,i))
    %drawnow
   
   
end

%c = max(c,[],3);

c(:,:,end) = 2*prior;
c = mean(c,3);
c = imnorm(double(c));

th = quantile(c(:),para.thresh);

b = im2bw(c,th);

b = imfill(b,'holes');
seg = imopen(b, openFilter);

b = imread(['zStackedYFP/' num2str(frame) '.png']);

border = (edge(seg,'log',0));


b = double(b);
b(border==1) = 1.1*max(b(:));

imagesc( b )
drawnow



%%

