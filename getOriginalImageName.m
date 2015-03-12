%getOriginalImageName(expe,colorName,i)
function fname = getOriginalImageName(expe,colorName,movie,frame,naming)

switch naming
    
    case 'biop'

        if movie<10
            fname = [expe.imgDir '/00' num2str(movie) ' ' colorName '.tif'];
        elseif movie<100
            fname = [expe.imgDir '/0' num2str(movie) ' ' colorName '.tif'];
        else
            fname = [expe.imgDir '/' num2str(movie) ' ' colorName '.tif'];
        end
        
    case 'bsf'

        t = expe.dt*60*60*1000*(frame-1);
                
        [movie_number,movie_letter] = ind2sub([6,8],movie);
        movie_letter = char(64 + movie_letter );
        
        movie_name = [movie_letter ' - ' n2s(movie_number)];

        fname = [expe.imgDir '/' movie_name '(fld 1 wv ' colorName ' - ' colorName '- time ' n2s(frame) ' - ' n2s(t) ' ms).tif'];

        
        if strcmp(colorName,'Trans')            
           fname = [expe.imgDir '/' movie_name '(fld 1 wv TL-Brightfield - Cy3- time ' n2s(frame) ' - ' n2s(t) ' ms).tif'];                       
        end
        
        
end
