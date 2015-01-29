%% experiment paramaters

cd /Users/bieler/Desktop/matlab/segmentation2

expe = experimentPara();

mainDir = expe.mainDir;
imgDir = expe.imgDir;
N = expe.numberOfFrames;

% Define some stuff, move into the right folder

movie = 1;

outDir = [mainDir '/movie' num2str(movie) '/'];
cd(outDir)

addpath ../
addpath ../code
addpath ../code/regionbased_seg
set(0,'defaultlinelinewidth',2)

%% combine stacks, denoise, ...

Nz = expe.numberOfColors; %number of images to combine

deNoise = {'none','BM3D','median','localNorm'}; %denoise algo on each stack
deNoise = deNoise{3};

medianSize = 3;

weightsSegmentation = [1 1 1]; %weights for summing the different channels
compressionQuantile = 0.99;      %signal above this quantile will be cut off, set to 1 to disable
gaussianFilterSize = 50;         %typycal length of the background

doDraw = 1;

mkdirIfNotExist('zStackedYFP');
mkdirIfNotExist('zStackedYFP_Data');

for k=1:N

    disp(100*k/N);
          
    images = {};
    for i=1:Nz
        imagesPath{i} = [outDir 'img/' getImageName(expe.colorNames{i},k)];
    end
            
    [out,N1,N2] = combineStack(imagesPath,Nz,deNoise,medianSize,compressionQuantile,gaussianFilterSize,weightsSegmentation,doDraw);
    imwrite(out,['zStackedYFP/' num2str(k) '.png']);        
end

%% just display the combined images

clf;colormap jet
for k=1:N
    
        disp(100*k/N)                 
        a = imread(['zStackedYFP/' num2str(k) '.png']);        
        imagesc(a); caxis([0 250])        
        pause(0.08); drawnow;        
end

%% test segmentation on a few frames 
       
guiSeg

%% Do the segmentation

inputFolder = 'zStackedYFP/';

load segPara.mat
load segParaFrames.mat

doDraw = 0;
segMethod = 1;
%

for frame=1:N
    disp(100*frame/N);
    
    para = [];
        
    para.minSize  = interp1(segParaFrames,[segPara(:).minSize  ],frame);
    para.maxSize  = interp1(segParaFrames,[segPara(:).maxSize  ],frame);
    para.thresh   = interp1(segParaFrames,[segPara(:).thresh   ],frame);
    para.openSize = interp1(segParaFrames,[segPara(:).openSize ],frame);
    para.ratioThresh = interp1(segParaFrames,[segPara(:).ratioThresh ],frame);
    
    %generateFilters;
        
    if(segMethod ==1)
        
        [filters openFilter] = generateFilters(para,0);
        segmentImage(frame,para,inputFolder,filters,openFilter,doDraw);
    else
    
        load ../segmentationCoeff.mat
        seg = segmentLDA(inputFolder,segPara,segParaFrames,frame,coeff,N1,N2);
        clf
        imagesc(seg)
        drawnow    
    end
end

%% do the measures

threshFolder = 'zStackedThresh/';
doMeasures(N,threshFolder,expe);

% try to split some merged cells

saveFolder = 'zStackedThreshSplit/';
mkdirIfNotExist(saveFolder)

doDraw = 0;

threshold = 1.6; %low value -> split everything

splitMergeCells;

threshFolder = saveFolder;
doMeasures(N,threshFolder,expe);

%% look at difference between splited and original

doDraw = 1;
if doDraw 
    for k=1:N
        
         a = imread(['zStackedThresh/' num2str(k) '.png']);
         b = imread(['zStackedThreshSplit/' num2str(k) '.png']);

         a = double(a);
         b = double(b);

         imagesc(a+b);
         pause(0.01)
    end
end

%% correct segmentation by hand

linksGui

%% redo the measures with corrected images

threshFolder = 'zStackedThreshCorrected/';
nObj = doMeasures(N,threshFolder,expe);

clf;
plot(nObj)
ylabel('number of objects')

%% tweak tracking parameters if needed
% main parameters are frame_displacement and split_cost

edit('get_struct.m');

%% load all measures & do the final Tracking

NToTrack = N;

doLinksOnly = 0;

clc
Me = loadMeasures(N);
[tracks, signal, traj, ind, divisions,divPerframe,trajX,trajY] = doTracking(NToTrack, Me,doLinksOnly);

clf; imagesc( signal(:,:,1) ); colormap jet

%% reload things if needed

load tracks.mat 
load signal.mat 
load traj.mat
load ind.mat 
load divisions.mat 

%% plot Tracking

pauseTime = 0.1;

plotTracking;

%% select good traces, based on length

lengthOfTrace = zeros(size(ind,1),1);
lengthOfGaps  = zeros(size(ind,1),1);

indAnnotation = zeros(size(ind));

if( exist('lengthThresh.mat','file') )
    load lengthThresh.mat;
else
    lengthThresh = 0.3; %note: to change lengthThresh value you first need to delete the file if it exists: !rm lengthThresh.mat
end

for i=1:size(ind,1)
        
    %look for continous traces
    indAnnotation(i,:) = markTrace(ind(i,:));

    lengthOfGaps(i) = sum( indAnnotation(i,:) == -3);
    lengthOfTrace(i) = sum( indAnnotation(i,:) == 1);
    
end

longTraces = find( (lengthOfGaps < 1) .* (lengthOfTrace/NToTrack>lengthThresh) );
clf

A = signal(longTraces,:);
[tmp,ia,ic] = unique(A,'rows');

longTraces = sort( longTraces(ia) );

imagesc(signal(longTraces,:,1))

length(longTraces)
colormap jet

save lengthThresh.mat lengthThresh
save longTraces.mat longTraces

%% build area, sum of signal, peak and div matrices

minTimeBetweenPeaks = 5;
peakMethod ='diff';
doDraw =0;

makePeakAndDivMatrices

%% refine area and signal around each cell, and do images for guiTraces

doDrawBkg = 0;  %display there area where the background is measured 
doDraw = 0;     %display area refinement result

bgkSize = 4;    %size around the cell where the background is not quantified
superSampling = 1; %increase the resolution of the image

mkdirIfNotExist('snapShots')

NToTrack = N;
a = imread(['zStackedYFP/' num2str(1) '.png']);
N1=size(a,1);
N2=size(a,2);

touchBorder = zeros(size(areaMatrix)); %is equale to one if the cell is touching the border of the image
refinedArea = zeros(size(areaMatrix));

refinedMean = zeros([size(areaMatrix) size(signal,3)]); %mean of the pixel values 
refinedSum = zeros([size(areaMatrix) size(signal,3)]);  %sum of the pixel values 
refinedStd  = zeros([size(areaMatrix) size(signal,3)]); %standard deviation of the pixel values 

bkg  = zeros([size(areaMatrix) size(signal,3)]);

refineAreaAndSignal

save refinedMean.mat refinedMean
save refinedSum.mat refinedSum
save bkg.mat bkg
save refinedArea.mat refinedArea
save touchBorder.mat touchBorder

%% plot mean signal of trace i with std

i=1

clfh
sel = ind(longTraces(i),:) > 0;
col = {'r','g'};
for k=1:min(2,expe.numberOfColors)
    errorbar(expe.t(sel),refinedMean(longTraces(i),sel,k),refinedStd(longTraces(i),sel,k),col{k});
end

xlabel('time')

%% plot mean signal, trace i 

i=1

clfh
sel = ind(longTraces(i),:) > 0;
plot(expe.t(sel),refinedMean(longTraces(i),sel,1),'r');
%plot(expe.t(sel),refinedMean(longTraces(i),sel,2),'g');


%% plot sum of signal, trace i 

i=1

clfh
sel = ind(longTraces(i),:) > 0;
plot(expe.t(sel),refinedSum(longTraces(i),sel,1),'r');
%plot(expe.t(sel),refinedSum(longTraces(i),sel,2),'g');

%% plot area and division trace i 

i=1

clfh
sel = ind(longTraces(i),:) > 0;
plot(expe.t(sel),imnorm(areaMatrix(longTraces(i),sel,1)),'k');
plot(expe.t(sel),0.5*divMatrix(longTraces(i),sel),'r');

legend('area','divisions')

%% plot trajectories 

clfh
for i=1:size(ind,1)
    sel= ind(i,:)>0;
    plot(trajX(i,sel),trajY(i,sel),'color',rand(3,1))
end

%% make small images around each cell for the trace tool 
% use if you don't want to do the refine area thing above

mkdirIfNotExist('snapShots')

doDraw = 0;
inputFolder = 'zStackedYFP/';
  
makeImagesForTraceTool


%% delete peakMatrixFinal and divMatrixFinal (reset guiTraces)

%!rm divMatrixFinal.mat
%!rm peakMatrixFinal.mat

%% Traces tool, correct divs and peaks

guiTraces

%% Load corrected peaks and divs

load divMatrixFinal.mat
load peakMatrixFinal.mat

clf
imagesc(peakMatrix(:,:,1)-divMatrix)

%%








