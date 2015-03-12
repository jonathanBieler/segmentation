if( binsize > 1 )
    
    
    
    
    system(['sudo chmod +w img/*.png']);
    
    if(exist('fullSizeimg','dir'))
        disp('fullSizeimg folder already exist, using that as input'); 
    else
        mkdir('fullSizeimg')
        system(['sudo chmod +w fullSizeimg/']);
        system(['cp -v img/*.png fullSizeimg']);
    end
      
    for j=1:expe.numberOfColors
        for k=1:N

            disp(100*k/N);

            a = imread([outDir 'fullSizeimg/' getImageName(expe.colorNames{j},k) ]);
            a = binImage(a,binsize);
            if(doDraw)
                clf; imagesc(a)
                drawnow
            end
                       
            imwrite(uint16(a/binsize),[outDir 'img/' getImageName(expe.colorNames{j},k) ]);
        end
    end
    if(expe.hasTrans)
        for k=1:N
            a = imread([outDir 'fullSizeimg/' getImageName(expe.transName,k) ]);
            a = binImage(a,binsize);
            if(doDraw)
                clf; imagesc(a)
                drawnow
            end
                       
            imwrite(imnorm(a),[outDir 'img/' getImageName(expe.transName,k) ]);
        end
    end
    
end