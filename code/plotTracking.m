%% video cell position

if(~exist('longTraces','var'))
   longTracesOnly = 0;
end
   

clf;
colormap hot

col = rand(length(traj),3);

col(:,1) = 0.7*col(:,1);
col(:,2) = 0.7*col(:,2);
col(:,3) = col(:,3).^0.5;


for k=1:1:NToTrack
    
    k
    clf;
    %a = imread(['zStackedThresh4/' num2str(k) '.png']);
    a = imread(['zStackedYFP/' num2str(k) '.png']);

    imagesc(double(a).^1);
    hold on
    
    div = divPerframe{k};
    
    for i=1:length(div)
        pos=traj{div(i)};
        
        text(pos(k,1),pos(k,2),'Boom!','color','green','fontSize',18)
        
    end
    
    idx = find(ind(:,k));
    
    
    for i=1:length(idx)
        
               
        j = idx(i);
        
        if( ~longTracesOnly || ismember(j,longTraces) )
        
            plot(trajX(j,k)+0*randn,trajY(j,k),'s','color',col(i,:),'markerSize',10)       

            if(k>1 && trajX(j,k-1)>0 && trajY(j,k-1)>0)

                plot([trajX(j,k-1) trajX(j,k)],[trajY(j,k-1) trajY(j,k)],'color',col(i,:))
            end        

        end
    end
    
    set(gca,'Ydir','normal')
    drawnow    
    pause(pauseTime);
end

hold off

