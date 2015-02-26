%function [] = segmentImage(k,segPara,inputFolder,filters,openFilter,showImages)
function [] = segmentImage(k,segPara,inputFolder,filters,openFilter,showImages)

b = imread([inputFolder num2str(k) '.png']);
%b = imread(['zStackedYFP_HighPassed/' num2str(k) '.png']);

if( segPara.maxSize > 0 )

    Nf=size(filters,3);
    c = zeros(size(b,1),size(b,2),Nf);

    for i=1:Nf

       tmp = imfilter(b,filters(:,:,i),'replicate');

       tmp(tmp<0)=0;

       %c(:,:,i) = imnorm(medfilt2(tmp,[2 2]));
       c(:,:,i) = imnorm(tmp);


       if(showImages)
        imagesc(c(:,:,i))
        drawnow
       end

    end

    c = mean(c,3);

else

    c = b;
    
end
%c = max(c,[],3);


c = imnorm(double(c));

th = quantile(c(:),segPara.thresh);

b = im2bw(c,th);

if(showImages)
    imagesc(b)
    drawnow
end

b = imfill(b,'holes');
b = imopen(b, openFilter);

if(showImages)
    imagesc(b)
    drawnow
end

mkdirIfNotExist('zStackedThresh')

out = ['zStackedThresh/' num2str(k) '.png'];
imwrite(b,out);

removeCytoplasm(k,inputFolder,segPara.ratioThresh)
