mkdirIfNotExist('testSegmentation')

n=1;
segPara = tPara(n);

generateFilters;

for i=1:length(testImages)

    frame = testImages(i);

    if(length(tPara) == 2) 
        segPara.minSize = interp1([1 N],[tPara(1).minSize tPara(2).minSize],frame);
        segPara.maxSize = interp1([1 N],[tPara(1).maxSize tPara(2).maxSize],frame);
        segPara.thresh = interp1([1 N],[tPara(1).thresh tPara(2).thresh],frame);
        segPara.openSize = interp1([1 N],[tPara(1).openSize tPara(2).openSize],frame);
    end

    disp(100*i/length(testImages))
    segmentImage(testImages(i),segPara,inputFolder,filters,openFilter,doDraw);

    %copy images    
    in = ['zStackedThresh/' num2str(testImages(i)) '.png'];
    out = ['testSegmentation/' num2str(testImages(i)) '_p_' num2str(n) '.png'];
    system(['cp ' in ' ' out]);
    
    a = imread(['zStackedYFP_data/' num2str(testImages(i)) '.png']);        
    b = imread(out);
    %superpose both
    b = (edge(b,'log',0));

    a(b==1) = min(a(:));                
    a =  ind2rgb(uint8(100*imnorm(a).^1),jet); 

    
    imagesc(a)
    drawnow

    imwrite(a,['testSegmentation/' num2str(testImages(i)) '_p_' num2str(n) 'Original.png'])

end

%% generate html page

fileID = fopen('testSegmentation/testSegmentation.html','w');

s = ['<html>\n' ...
     '<head>\n' ...
     '<style type="text/css">\n' ...
     'body {color: #555555; font: 1em/1.3 "Gill Sans", "My Gill Sans", sans-serif;}\n' ...
     'img {float:left; margin:5px; max-width: 1700px; max-height: 1700px;}\n' ...      
     '.mainDiv {border-radius: 0px;background:rgb(240,240,240); padding:5px; margin:10px;\n'...
     '-moz-box-shadow: 3px 3px 5px 0px #ccc; -webkit-box-shadow: 3px 3px 5px 0px #ccc;box-shadow: 3px 3px 5px 0px #ccc;">\n}' ...
     '\n'...
     '.dateTitle {color:rgb(150,150,150); text-shadow: #ffffff 1px -1px 1px; font: 1.2em/1.1 "OFL Sorts Mill Goudy", "Georgia", serif; letter-spacing: 0.05em;}\n'...
     '\n' ...
     'hr {color:rgb(20,20,20);}'...
     '.clear-left {clear:left;}'...
     '</style>\n\n' ...
     '</head>\n\n' ...     
     '<body>\n']; 
 
fprintf(fileID,s);
   
for n=1:1
     
     
    s = ['<div class="mainDiv">' ...
    '<span class="dateTitle">parameter set#' num2str(n) '</span>\n' ...    
    '<hr />\n'];

     fprintf(fileID,s);
     
     
     for i=1:length(testImages)

         src = [num2str(testImages(i)) '_p_' num2str(n) 'Original.png'];
         s = ['<img src="' src '"></img>\n'];
         fprintf(fileID,s);
         
         
         s = ['<div class="clear-left"> </div>'];
         fprintf(fileID,s);
     end
     

    s = ['</div>\n'];  
    fprintf(fileID,s);
     
end
 
 
 
s = ['</body>\n'...
     '</html>\n'];
fprintf(fileID,s);

fclose(fileID);
!open testSegmentation/testSegmentation.html