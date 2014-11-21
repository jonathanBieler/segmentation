function coeff = trainSegmentation(inputFolder,segPara,segParaFrames,frames)

doDraw=0;

features = [];
target = [];

for frame = frames

    frame

    para = [];

    para.minSize  = interp1(segParaFrames,[segPara(:).minSize  ],frame);
    para.maxSize  = interp1(segParaFrames,[segPara(:).maxSize  ],frame);
    para.thresh   = interp1(segParaFrames,[segPara(:).thresh   ],frame);
    para.openSize = interp1(segParaFrames,[segPara(:).openSize ],frame);
    para.ratioThresh = interp1(segParaFrames,[segPara(:).ratioThresh ],frame);

    %generateFilters;

    [filters openFilter] = generateFilters(para,doDraw);


    b = ( imread([inputFolder num2str(frame) '.png']) );
    %b = imread(['zStackedYFP_HighPassed/' num2str(k) '.png']);

    Nf=size(filters,3);
    c = zeros(size(b,1),size(b,2),Nf);
    f = zeros( size(b,1)*size(b,2), Nf+4);

    for i=1:Nf

       tmp = imfilter(b,filters(:,:,i),'replicate');

       %tmp(tmp<0)=0;

       %c(:,:,i) = imnorm(medfilt2(tmp,[2 2]));
       c(:,:,i) = imnorm(tmp);


       if(doDraw)
        imagesc(c(:,:,i))
        drawnow
       end

       f( : ,i) = imnorm(tmp(:));

    end

    tmp = imfilter(b,fspecial('gaussian',100,8),'replicate' );
    f(:,end-3) = tmp(:);

    tmp = b-imfilter(b,fspecial('gaussian',100,8),'replicate' );
    f(:,end-2) = tmp(:);

    tmp = imfilter(b,fspecial('log',50,0.2),'replicate' );
    f(:,end-1) = tmp(:);

    tmp = imfilter(b,fspecial('laplacian',0.5),'replicate');
    f(:,end) = tmp(:);

    %c = max(c,[],3);
    c = mean(c,3);

    c = imnorm(double(c));

    th = quantile(c(:),para.thresh);

    seg2 = im2bw(c,th);


    seg2 = imfill(seg2,'holes');
    seg2 = imopen(seg2, openFilter);


    trueSeg = imread(['zStackedThreshCorrected/' n2s(frame) '.png']);
    target = [target; trueSeg(:)];

    features = [features; f];

end

[class,err,post,logp,coeff]  = classify(features,features,target,'linear');