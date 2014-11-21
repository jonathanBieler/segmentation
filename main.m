%% experiment paramaters

cd /Users/bieler/Desktop/matlab/segmentation2

expe = experimentPara();

% Define some stuff, move into the right folder

movie = 21;

outDir = [mainDir '/movie' num2str(movie) '/'];
cd(outDir)


addpath ../
addpath ../code
addpath ../code/regionbased_seg


%% combine stacks, denoise, ...

Nz = expe.numberOfColors; %number of images to combine

deNoise = {'none','BM3D','median','localNorm'}; %denoise algo on each stack
deNoise = deNoise{3};

medianSize = 3;

weightsSegmentation = [1 0.2 1]; 

doDraw = 1;

mkdirIfNotExist('zStackedYFP');
mkdirIfNotExist('zStackedYFP_Data');


for k=1:N

    disp(100*k/N);
          
    images = {};
    for i=1:Nz
        imagesPath{i} = [outDir 'img/' getImageName(expe.colorNames{i},k)];
    end
            
    [out,N1,N2] = combineStack(imagesPath,Nz,deNoise,medianSize,weightsSegmentation,weightsData,doDraw);
    imwrite(out,['zStackedYFP/' num2str(k) '.png']);
        
end

%% just display the combined images

clf;colormap jet
for k=1:N
    
        disp(100*k/N)                 
        a = imread(['zStackedYFP/' num2str(k) '.png']);        
        imagesc(a); caxis([0 150])        
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
        
        [filters openFilter] = generateFilters(para,doDraw);
        segmentImage(frame,para,inputFolder,filters,openFilter,0);
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

threshold = 1.5; %low value -> split everything

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
         pause(0.1)
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

%% load all measures & do the final Tracking

NToTrack = N;

doLinksOnly = 0;

Me = loadMeasures(N);
[tracks, signal, traj, ind, divisions,divPerframe,trajX,trajY] = doTracking(NToTrack, Me,doLinksOnly);

clf; imagesc( signal(:,:,1) )

%%

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
    lengthThresh = 0.5;
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

%% build area, peak and div matrices

minTimeBetweenPeaks = 5;
peakMethod ='diff';
doDraw =0;

makePeakAndDivMatrices

%% refine area and do little images

doDraw = 1;
superSampling = 1;

mkdirIfNotExist('snapShots')

NToTrack = N;
a = imread(['zStackedYFP/' num2str(1) '.png']);
N1=size(a,1);
N2=size(a,2);

touchBorder = zeros(size(areaMatrix));
refinedArea = zeros(size(areaMatrix));

refinedMean = zeros([size(areaMatrix) size(signal,3)]);
refinedStd  = zeros([size(areaMatrix) size(signal,3)]);

bkg  = zeros([size(areaMatrix) size(signal,3)]);

for n=1:length(longTraces)
        
    idx = longTraces(n);        

    clf;

    s = 40;

    w=2*s+1;
    nR = 6;
    tot = NToTrack;

    d = round(tot/nR);
    out = zeros(d*w,nR*w);
    
    %i = round(traj{idx}(:,1));
    %j = round(traj{idx}(:,2));
    
    for k=1:N

        k
        if( ind(idx,k)~=0 )
            a = imread(['zStackedYFP/' num2str(k) '.png']);

            data = {};
            for j=1:expe.numberOfColors            
                data{j} = imread( ['img/' getImageName(expe.colorNames{j},k)] );
            end
                        
            m = imread(['zStackedThreshCorrected/' num2str(k) '.png']);
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
            bgkMask = imdilate(bgkMask,strel('disk',3));   

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
            
            %m = imdilate(m,strel('disk',1));   
    
            if( superSampling > 1)
                sub_a = imresize(sub_a,superSampling);
                
                for j=1:expe.numberOfColors  
                    sub_data{j} = imresize(sub_data{j},superSampling);
                end 
                m = imresize(m,superSampling);
            end
           
            if(updateArea)
                seg = region_seg(sub_a, m, 12,0.8,doDraw && mod(k,1)==0); %-- Run segmentation
            else
                seg = m;   
            end

            %
            tmp = imclearborder(seg);

            if( sum(tmp(:)) > 0 )
                 seg = tmp;
            end

            %quantify green
            
            for j=1:expe.numberOfColors  
                
                tmp = double(sub_data{j}(seg==1));            
                qu = quantile(tmp,0.95);    ql = quantile(tmp,0.05);
                tmp = tmp(  (tmp < qu) & (tmp > ql) );

                refinedMean(idx,k,j)    = mean(tmp);            
                refinedStd(idx,k,j) = std(tmp);
            
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
    subplot_tight(2,1,2,0.1); hold on
        %plot(refinedMean(idx,:)-refinedStd(idx,:),'r--');
        %plot(refinedMean(idx,:)+refinedStd(idx,:),'r--');

        plot(refinedMean(idx,:,1),'r');
        plot(refinedMean(idx,:,2),'g');
        plot(signal(idx,:,1),'k')
        plot(signal(idx,:,2),'k')
        
        plot(bkg(idx,:,1),'r--')
        plot(bkg(idx,:,2),'g--')
    drawnow
    
end


save refinedMean.mat refinedMean
save bkg.mat bkg
save refinedArea.mat refinedArea
save touchBorder.mat touchBorder


%% plot trace i

i=1

clfh
sel = ind(longTraces(i),:) > 0;
errorbar(refinedMean(longTraces(i),sel,1),refinedStd(longTraces(i),sel,1),'r');
errorbar(refinedMean(longTraces(i),sel,2),refinedStd(longTraces(i),sel,2),'g');

%% make small images around each cell for the trace tool 
% use if you don't want to do the refine area thing above

mkdirIfNotExist('snapShots')

doDraw = 0;
inputFolder = 'zStackedYFP/';
  
makeImagesForTraceTool

%% Traces tool, correct divs and peaks

guiTraces

%% Load corrected peaks and divs

load divMatrixFinal.mat
load peakMatrixFinal.mat

clf
imagesc(peakMatrix-divMatrix)

%%

load signalr.mat signalr
load signalg.mat signalg
load ind.mat

%%

ind = ind(longTraces,:);
signalr = signalr(longTraces,:);
signalg = signalg(longTraces,:);
areaMatrix = areaMatrix(longTraces,:);
peakMatrix = peakMatrix(longTraces,:);
divMatrix  = divMatrix(longTraces,:);

%%

sel = (sum(peakMatrix,2) + sum(divMatrix,2)) > 0;

ind = ind(sel,:);
signalr = signalr(sel,:);
signalg = signalg(sel,:);
areaMatrix = areaMatrix(sel,:);
peakMatrix = peakMatrix(sel,:);
divMatrix  = divMatrix(sel,:);


%% export all to dropbox

save /Users/bieler/Dropbox/NicolasVilla/3t3FucciData/movie3/ind.mat ind
save /Users/bieler/Dropbox/NicolasVilla/3t3FucciData/movie3/signalr.mat signalr
save /Users/bieler/Dropbox/NicolasVilla/3t3FucciData/movie3/signalg.mat  signalg
save /Users/bieler/Dropbox/NicolasVilla/3t3FucciData/movie3/areaMatrix.mat  areaMatrix
save /Users/bieler/Dropbox/NicolasVilla/3t3FucciData/movie3/peakMatrix.mat  peakMatrix
save /Users/bieler/Dropbox/NicolasVilla/3t3FucciData/movie3/divMatrix.mat   divMatrix

%% cluster

sel = ind>0;

normFun = @(x) (x/1);
%normFun = @(x) imnorm(x);

normSignalr = signalr;
normSignalg = signalr;

L = 40;
smoothing = 40;
p = 0.01;

x=[]; y=[];
for i=1:size(ind,1)
    
    sel = ind(i,:)>0;
    %x = [x normFun(signalg(i,sel))];
    %y = [y normFun(signalr(i,sel))];    
    
    %    
    tmp = signalr(i,sel);          
    [mi ma] = getMovingMinMax(tmp,L,smoothing,p);    
    y = [y (tmp-mi')./ma']; 
    
    normSignalr(i,sel) = (tmp-mi')./ma';
    
    
    
%     clfh
%     plot(tmp)
%     plot(mi,'r')
%     plot(ma,'r')
%     plot( 7000*(tmp-mi')./ma','k' )
%     pause(1)
    
    tmp = signalg(i,sel);          
    [mi ma] = getMovingMinMax(tmp,L,smoothing,p);    
    x = [x (tmp-mi')./ma']; 
    
    normSignalg(i,sel) = (tmp-mi')./ma';


end
%

clfh
scatterWithDensity(x,y)
colormap jet

options = statset('Display','final');
obj = gmdistribution.fit([x' y'],3,'Options',options);

mu    = obj.mu;
sigma = obj.Sigma;

for i=1:size(mu,1)
    h = plot_gaussian_ellipsoid(mu(i,:), sigma(:,:,i));
    set(h,'color','k'); 
end

c = [0.15 0.15];
plot(c(1),c(2),'ks')


paperSize(24,18)
setFonts

print -dpdf cellCyclePhasesClassification.pdf
!open cellCyclePhasesClassification.pdf

%% decode trace using GMU

i = 1;
clf

sel = ind(i,:)>0;
x = [normFun(normSignalg(i,sel))];
y = [normFun(normSignalr(i,sel))];    

c = [0.1 0.15];

th = -atan2(y-c(2),x-c(1));

idx = cluster(obj,[x(:) y(:)]);

clfh
plot(0.1+x,'r')
plot(0.1+y,'g')
%plot(0.5+th/6,'k')

L = {'G1','G2','S'};
for j=1:4:length(idx)
   
   % text(j,0.05,L{idx(j)})    
end

axis([0 length(x)+1 0 1.1])
xlabel('Temps')

%%

i = 1;
clf

sel = ind(i,:)>0;
x = [normFun(normSignalg(i,sel))];
y = [normFun(normSignalr(i,sel))];   

plot(x,y,'k')

%%

i = 1;
clf

sel = ind(i,:)>0;
x = [normFun(normSignalg(i,sel))];
y = [normFun(normSignalr(i,sel))];   

for t=1:length(x)
    
    sel = max(1,t-10):min(length(x),t);

    clfh
    plot(x,y,'k','lineWidth',1)
    plot(x(sel),y(sel),'r')
    plot(x(t),y(t),'rs')
        
    plot(c(1),c(2),'ks')
    plot([c(1) x(t)],[c(2) y(t)],'k--')
    
    
    
    pause(0.05)
    
end

%% plot for thesis defense

for k=1:6

    Nt = 1e3;
    dt= 1e-1;
    v = zeros(Nt,1);
    x = zeros(Nt,1);

    t = linspace(0,Nt*dt,Nt);

    v(1) = 1;

    g = 1 + (k-3)/10;

    v = zeros(Nt,1);
    x = zeros(Nt,1);

    t = linspace(0,Nt*dt,Nt);

    v(1) = 1;

    for t=1:Nt-1

        x(t+1) = x(t) + dt*v(t);
        v(t+1) = v(t) - g*0.02*dt;

    end

    t = linspace(0,Nt*dt,Nt);

    clfh
    plot(t(1:20:end),x(1:20:end),'r')


    v(1) = 1;

    for t=1:Nt-1

        x(t+1) = x(t) + dt*v(t);
        v(t+1) = v(t) - 0.02*dt;

    end

    t = linspace(0,Nt*dt,Nt);


    %plot(t(1:20:end),x(1:20:end),'ks')


    axis([0 100 0 30])
    xlabel('distance')
    ylabel('hauteur')

    setFonts
    paperSize(24,18)

    print('-dpdf',['tmp' n2s(k) '.pdf'])

end


%%

i = 3;
clf

sel = ind(i,:)>0;
x = [normFun(normSignalg(i,sel))];
y = [normFun(normSignalr(i,sel))];   
z = [normFun(areaMatrix(i,sel))];   

for t=1:length(x)
    
    sel = max(1,t-10):min(length(x),t);

    clfh
    plot3(x,y,z,'k','lineWidth',1)
    plot3(x(sel),y(sel),z(sel),'r')
    plot3(x(t),y(t),z(t),'rs')
        
    %plot(c(1),c(2),'ks')
    %plot([c(1) x(t)],[c(2) y(t)],'k--')
    
    
    view(-150+t/10,14)
    axis([-0.2 0.8 -0.2 0.8 0 2000])
    axis vis3d
    
    pause(0.01)
    
end


%%

th_idx_1 = [];
th_idx_2 = [];
th_idx_3 = [];
th_idx_4 = [];

for i=1:size(signalg,1)
    
    sel = ind(i,:)>0;
    x = [normFun(signalg(i,sel))];
    y = [normFun(signalr(i,sel))];    

    th = -atan2(y-c(2),x-c(1))+pi;

    idx = cluster(obj,[x(:) y(:)]);

    th_idx_1 = [th_idx_1; th(idx==1)'];
    th_idx_2 = [th_idx_2; th(idx==2)'];
    th_idx_3 = [th_idx_3; th(idx==3)'];
    th_idx_4 = [th_idx_4; th(idx==4)'];
    
end

clfh
bins = linspace(0,2*pi,30);

n  = ksdensity(th_idx_1,bins,'width',0.3);
plot(bins,n,'b')
n  = ksdensity(th_idx_2,bins,'width',0.3);
plot(bins,n,'r')
n  = ksdensity(th_idx_3,bins,'width',0.3);
plot(bins,n,'g')
n  = ksdensity(th_idx_4,bins,'width',0.3);
plot(bins,n,'k')

xlabel('Cell cycle phase')
ylabel('density')

axis([0 2*pi 0 1.4])

paperSize(24,18)
setFonts

print -dpdf cellCyclePhases.pdf
!open cellCyclePhases.pdf

%% 


normFun = @(x) imnorm(x);

sel = ind>0;

x=[]; y=[]; z=[]; t=[];

for i=1:size(ind,1)
 
    sel = ind(i,:)>0;
    x = [x normFun(signalr(i,sel))];
    y = [y normFun(signalg(i,sel))];
    z = [z normFun(areaMatrix(i,sel))];
        
end

clf
scatter3(x,y,z,30,z,'filled')
colormap winter







