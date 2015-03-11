%getOriginalImageName(expe,colorName,i)
function fname = getOriginalImageName(expe,colorName,movie,frame)


%         if movie<10
%             fname = [expe.imgDir '/00' num2str(movie) ' ' colorName '.tif'];
%         elseif movie<100
%             fname = [expe.imgDir '/0' num2str(movie) ' ' colorName '.tif'];
%         else
%             fname = [expe.imgDir '/' num2str(movie) ' ' colorName '.tif'];
%         end

    t = expe.dt*60*60*1000*(frame-1);

    fname = [expe.imgDir '/B - 4(fld 1 wv ' colorName ' - ' colorName '- time ' n2s(frame) ' - ' n2s(t) ' ms).tif'];