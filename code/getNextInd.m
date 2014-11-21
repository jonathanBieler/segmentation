%function parents = getNextInd(f,ind,tracks,maxGap)
function parents = getNextInd(f,ind,tracks,maxGap)

parents = zeros(0);

np=1;

for k=f+1:min(f+maxGap,length(tracks))

    c = tracks(k).cluster;
    
    idx = find(c(:,2) == ind);
    
%     if(length(idx)>1)
%         'found more than one index'
%         idx = idx(1);
%     end


    for i=1:length(idx)

        if( ~isempty(idx(i)) &&  ~isempty(c(idx(i),3)) && c(idx(i),3) == f  )

            parents(np).index = c(idx(i),1);
            parents(np).frame = k;
            np = np+1;

        end

    end
    
    
end
    