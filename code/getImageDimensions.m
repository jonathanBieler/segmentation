function [N1,N2] = getImageDimensions(expe)


    a = imread(['img/' getImageName(expe.colorNames{1},1)]);
    [N1,N2]=size(a);
