%function [] = write_image(name,size,res=300)
function [] = write_image(name,size,res)

if ( nargin < 3 )    
    
    res = 300;
end

if ( nargin < 3 )    
    size= 0.5;
    res = 300;
end


A = print2array(gcf,res);

if (size ~= 1)
    A = imresize(A,size);
end

imwrite(A,name);