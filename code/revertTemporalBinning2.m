tb =  round( length( dir('zStackedYFP_unbinned/*.png') ) ./ length( dir('zStackedYFP/*.png') ));
No = expe.numberOfFrames;


if( tb > 1 )

%%

ind2 = zeros(size(ind,1),No);

for i=1:size(ind,1)    
    tmp = reshape(repmat(ind(i,:),temporalBinning,1),1,[]);
    ind2(i,1:length(tmp)) = tmp;        
end

%fill missing values 
idx = size(ind,2)*temporalBinning;
if( idx < No )
    ind2(:,(idx+1):end) = repmat(ind(:,end),1,No-idx);
end
imagesc(ind2)

%% Unbing matrices

[i j] = meshgrid(1:size(signal,2),1:size(signal,1));
[ip jp] = meshgrid(linspace(1,size(signal,2),No),1:size(signal,1)) ;

%ind = interp2(i,j,ind,ip,jp,'nearest');
signal_ = zeros( [size(ip) expe.numberOfColors] );

for k=1:expe.numberOfColors    
    signal_(:,:,k) = interp2(i,j,signal(:,:,k),ip,jp,'nearest');    
end

trajX_  = zeros( [size(ip) ] );
trajY_  = zeros( [size(ip) ] );
for i=1:size(trajX_,1)

    
    sel = find(ind(i,:) > 0);
    sel2 = find(ind2(i,:) > 0);

    if(length(sel) > 1)

        trajX_(i,sel2) = interp1(sel,trajX(i,sel),linspace(sel(1),sel(end),length(sel2)));
        trajY_(i,sel2) = interp1(sel,trajY(i,sel),linspace(sel(1),sel(end),length(sel2)));
    else

        trajX_(i,sel2) = trajX(i,sel);
        trajY_(i,sel2) = trajY(i,sel);
    end

end

imagesc(trajX_)

%% interpolate segmentation


idxBin = 1:No/temporalBinning;
k = reshape(repmat(idxBin,temporalBinning,1),1,[]);

k = [k ones(No-length(k),1)*k(end)]; %when temporalBinning doesn't divide No nicely

for i = 1:No
    

    if( i <= length(k) )           
        M1 = Me{  k(i)  };
        M2 = Me{ min(k(i) +1,k(end)) };
    else
        M1 = Me{ k(end) };
        M2 = Me{ k(end) };
    end
    
    idx = ind(:,k(i)); 
    a = zeros(N1,N2);
   
    for j=1:length(idx)   
        
        if(idx(j)>0)
        
            p = M1(idx(j)).PixelIdxList;
            [ip,jp] = ind2sub([N1,N2],p);

            ip = ip - M1(idx(j)).Centroid(2);
            jp = jp - M1(idx(j)).Centroid(1);

            ip = ip + trajY_(j,i);
            jp = jp + trajX_(j,i);

            ip = round(ip);
            jp = round(jp);

            ip(ip<1) = 1; ip(ip>N1) = N1;
            jp(jp<1) = 1; jp(jp>N2) = N2;

            p = sub2ind([N1,N2],ip,jp);

            a(p) = 1;       
        
        end
    end
    
    if(doDraw)
        clf
        imagesc(a)
        title([n2s(i) '/' n2s(No)])
        drawnow
    end
    
    imwrite(a,['zStackedThreshCorrected/' num2str( i  ) '.png']);
  
end

%% replace stacked images

system('rm zStackedYFP/*.png')
system('cp zStackedYFP_unbinned/*.png zStackedYFP/')

%% redo the measures with corrected images

N = No;

threshFolder = 'zStackedThreshCorrected/';
nObj = doMeasures(N,threshFolder,expe);


%% load all measures & do the final Tracking

NToTrack = N;

doLinksOnly = 0;

clc
Me = loadMeasures(N);
[tracks, signal, traj, ind, divisions,divPerframe,trajX,trajY] = doTracking(NToTrack, Me,doLinksOnly);




end
