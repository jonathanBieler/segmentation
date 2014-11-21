
expe = experimentPara();
cd(mainDir)
addpath code/

%%

for i = expe.indexOfFirstMovie:(expe.indexOfFirstMovie + expe.numberOfMovies  -1)
   
    mkdirIfNotExist([mainDir 'movie' num2str(i)]);
    mkdirIfNotExist([mainDir 'movie' num2str(i) '/img/']);
    
    for j = 1:length(expe.colorNames)
        
        colorName = expe.colorNames{j};
        
        if i<10
            fname = [imgDir '/00' num2str(i) ' ' colorName '.tif'];
        elseif i<100
            fname = [imgDir '/0' num2str(i) ' ' colorName '.tif'];
        else
            fname = [imgDir '/' num2str(i) ' ' colorName '.tif'];
        end
    
        fname = ['''' fname ''''];
            

        system(['cp -v ' fname ' ' mainDir 'movie' num2str(i) '/img/']);

    
    end
        
        
end

%% unpack tif files

doDraw = 1 ;

for i = expe.indexOfFirstMovie:(expe.indexOfFirstMovie + expe.numberOfMovies  -1)

    disp(i);
    
    for j = 1:length(expe.colorNames)
        
        colorName = expe.colorNames{j};
    
        if i<10
            fname = [mainDir 'movie' num2str(i) '/img/00' num2str(i) ' ' colorName '.tif'];
        elseif i<100
            fname = [mainDir 'movie' num2str(i) '/img/0' num2str(i) ' ' colorName '.tif'];
        else
            fname = [mainDir 'movie' num2str(i) '/img/' num2str(i) ' ' colorName '.tif'];
        end

        info = imfinfo(fname);
        num_images = numel(info);
        for k = 1:num_images

            zNumber = j;
            fNumber = k;

            A = imread(fname, k);

            if(doDraw && mod(k,5)==0)
                imagesc(A)
                drawnow
            end

            imgName = getImageName(colorName,fNumber);

            name = [mainDir 'movie' num2str(i) '/img/' imgName];
            %name
            imwrite(A,name)

        end
        
        system(['rm "' fname '"']);
    
    end

end
   


