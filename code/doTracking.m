function [tracks, signal, traj, ind, divisions, divPerframe, trajX, trajY] = doTracking(N, Me,doLinksOnly)


weigth = 0.03*[1/5 1/8 1]; %weights of the different signals (two here) and the area 

%define data given to the tracking algo
nuclei = struct;
for i=1:N    
   m=Me{i};
    
   nuclei(i).carth = cat(1,m.Centroid);
   
   p = zeros( [length(m) size(m(1).MeanIntensity,1)+1 ]); 
   for j=1:length(m)
       
       for k=1:size(m(j).MeanIntensity,1)
          p(j,k) =  m(j).MeanIntensity(k) * weigth(k);
       end
       p(j,end) = m(j).Area * weigth(end);
   end
   
   nuclei(i).properties = p;
end

disp(['Typical coordinates values'])
disp(nuclei(1).carth)
disp(['Typical properties values (signal and area)'])
disp(nuclei(1).properties)
fprintf('If the values in properties are larger \nthan carthesian coordinates, edit the weights doTracking.m\n\n')

%
save nuclei.mat nuclei;

% do the tracking 

addpath ../code/Tracking/
addpath ../code/Tracking/libraries/

doDraw = 0; 

%%

[tracks ] = tracking(nuclei,doLinksOnly);

if( doLinksOnly )
    tracks=0;
    signal=0;
    traj=0;
    ind=0;
    divisions=0;
    divPerframe=0;
    trajX=0;
    trajY=0;
    return
end

disp('done')

%% Build forward traj

%init
c = tracks(2).cluster; 
t = tracks(1).carth; 

ind = zeros(5000,length(tracks));

nTracks = length(c(:,2));

ind(1:nTracks) = c(:,2);

maxGap = 3;

divisions = [];
nDiv = 0;

for f=1:length(tracks)
    
    %f
    %nTracks
    
    
    %pause(1)
    
    
    c = tracks(f).cluster;
    
    tmp = setdiff(c(:,1), ind(1:nTracks,f));
    
    for j=1:length(tmp)
       
        ind(nTracks+1,f)=tmp(j); 
        
        %trace  = backTrace(f,tmp(j),tracks);
        nTracks = nTracks+1;
                       
%        ind(j+nTracks,1:f) = trace(1:f);
        
        %plot(trace(1:f))
        %drawnow
        
    end
       
    for i=1:nTracks
       
        parent = getNextInd(f,ind(i,f),tracks,maxGap);
                
        if(~isempty(parent))
                    
            ind(i,parent(1).frame) = parent(1).index;
            
%             if( parent(1).frame ~= f+1)
%                ind(i,(f+1):parent(1).frame-1) = -1;
%             end
                        
            %trajx = t(ind,1);
            %trajy = t(ind,2);        
        end  
                        
        if(length(parent) > 1)
           %parent(1)
           %parent(2)
           
           for j=2:length(parent)
              %ind(nTracks+1,1:f) = ind(i,1:f);

              ind(nTracks+1,parent(j).frame) = parent(j).index; 
              
              nTracks = nTracks+1;
              
              nDiv = nDiv+1;
              divisions(nDiv).motherFrame = f;
              divisions(nDiv).motherInd = i;
              divisions(nDiv).sisterFrame = parent(j).frame;
              divisions(nDiv).sisterInd = nTracks;
  
           end
        end
        
    end
        
end

ind = ind(1:nTracks,:);

%% get trajectories and signal

[traj signal] = getTrajFromInd(ind,tracks,Me);

%% fill up gaps in trajectories


% for k=1:length(traj)
% 
%     X = traj{k};
% 
%     openGap = 0;
%     wasZero = 0;
% 
%     if(X(1,1) == 0)
%        wasZero=1; 
%     end
% 
%     for i=2:length(X)
% 
%         if(X(i,1) == 0 && ~wasZero && ~openGap)
%             openGap = 1;
%             idx = [];
%             firstPoint(1) = X(i-1,1);
%             firstPoint(2) = X(i-1,2);
%         end
% 
%         if(X(i,1) ~= 0 && openGap)
%             openGap = 0;
%             lastPoint(1) = X(i,1);
%             lastPoint(2) = X(i,2);
%             
%             tmp = interp1([idx(1)-1 i],[firstPoint(1) lastPoint(1)],idx);
%             X(idx,1) = tmp;
%             
%             tmp = interp1([idx(1)-1 i],[firstPoint(2) lastPoint(2)],idx);
%             X(idx,2) = tmp;
%             
%         end
% 
%         if(X(i,1) == 0)
%             wasZero=1; 
%         else
%             wasZero=0; 
%         end
% 
%         if(openGap)
%             idx = [idx i];
%         end
% 
%     end
%     
% %     clf
% %     plot(X(:,2))
% %     pause
%     traj{k}=X;
% 
% end

%plot(X(:,1),X(:,2),'k')

%plot(X(:,1))
%plot(idx,X(idx,1),'r.')


%% get valid divisions

refF = 5;%length of refractory period in frames

valDiv = [];
nValid=1;
for i=1:length(divisions)

    isInRef = 0;
    
    div = divisions(i);

    tmp1 = ind(div.motherInd,:);
    tmp2 = ind(div.sisterInd,:);
    
    dd = find([divisions.motherInd] == div.motherInd);
    for j=1:length(dd)        
        if( div.motherFrame - divisions(dd(j)).motherFrame < refF && div.motherFrame - divisions(dd(j)).motherFrame > 0)
           isInRef = 1; 
        end
    end

    %check trajectory length (don't work anymore)
    if( sum(tmp1==0) < N && sum(tmp2==0) < N)

        %check merge after division
        tmp1 = tmp1(div.motherFrame:end);
        tmp2 = tmp2(div.motherFrame:end);

        m=tmp1./tmp2; m=sum(m==1);
        
        %refactory period
 
        if( 100*m / length(tmp1) < 5 && ~ isInRef) 

            valDiv(nValid) = i;
            nValid = nValid+1;

        end
    end
end



%% get division per frames

divPerframe = cell(1,N);

for i=1:length(valDiv)
    
    div = divisions(valDiv(i));
    tmp = divPerframe{div.motherFrame};
    
    tmp = [tmp div.motherInd];
    
    divPerframe{div.motherFrame} = tmp;
    
end

divDistrib = zeros(1,N);

for i=1:N   
    divDistrib(i) = length(divPerframe{i});
end

% clf
% plot(divDistrib)
% 
% drawnow;
% %pause(3)

%% get cell position in the same format that ind

trajX = zeros(size(ind));
trajY = zeros(size(ind));

for i=1:size(ind,1);
    
    pos = traj{i};
    trajX(i,:) = pos(:,1);
    trajY(i,:) = pos(:,2);
    
end

%%



save tracks.mat tracks
save signal.mat signal
save traj.mat traj
save ind.mat ind
save divisions.mat divisions

