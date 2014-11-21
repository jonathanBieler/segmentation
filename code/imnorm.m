function im = imnorm(im)

im = double(im);

im = im-min(min(min(im)));

if( max(im(:)) > 0 )
    im = im / max(im(:));
else
    disp( 'imnorm: max of image is zero' );
end