methods = {'constant','AND','OR'};
method = methods{1};

tb =  round( length( dir('zStackedYFP_unbinned/*.png') ) ./ length( dir('zStackedYFP/*.png') ));  

if( tb > 1 )

No = expe.numberOfFrames;

idxBin = 1:No/temporalBinning;
k = reshape(repmat(idxBin,temporalBinning,1),1,[]);

mkdirIfNotExist('tmp')

for i = 1:No
    
    %copy segmentation
    if( i <= length(k) )           
        a = imread(['zStackedThreshCorrected/' num2str( k(i)  ) '.png']);
        b = imread(['zStackedThreshCorrected/' num2str( min(k(i) +1,k(end)) ) '.png']);
    else
        a = imread(['zStackedThreshCorrected/' num2str( k(end)  ) '.png']);
        b = imread(['zStackedThreshCorrected/' num2str( k(end)  ) '.png']);
    end
    
    switch method
        
        case 'constant'
        
        case 'AND'
            a = a+b > 1;
        case 'OR'
            a = a+b > 0;
    end
    
    imagesc(a)
    drawnow
    colormap jet
    
    imwrite(a,['tmp/' num2str( i  ) '.png']); 
                
end

system('rm zStackedThreshCorrected/*.png')
system('cp tmp/*.png zStackedThreshCorrected/')
system('rm tmp/*.png')

%% same thing with combined images

system('rm zStackedYFP/*.png')
system('cp zStackedYFP_unbinned/*.png zStackedYFP/')

%% redo measure on unbinned segmentation

threshFolder = 'zStackedThreshCorrected/';
doMeasures(No,threshFolder,expe);

Me = loadMeasures(No);


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
%imagesc(ind2)

ind  = ind2;


%% Unbing matrices

[i j] = meshgrid(1:size(signal,2),1:size(signal,1));
[ip jp] = meshgrid(linspace(1,size(signal,2),No),1:size(signal,1)) ;

%ind = interp2(i,j,ind,ip,jp,'nearest');
signal = interp2(i,j,signal,ip,jp,'nearest');
trajX = interp2(i,j,trajX,ip,jp,'linear');
trajY = interp2(i,j,trajY,ip,jp,'linear');

N = No;

clf; imagesc(signal); colormap jet

%% transform divisions into unbinned time

for i=1:length(divisions)
   
    divisions(i).motherFrame = (divisions(i).motherFrame-1)*temporalBinning + 1;
    divisions(i).sisterFrame = (divisions(i).sisterFrame-1)*temporalBinning + 1;
    
end

end