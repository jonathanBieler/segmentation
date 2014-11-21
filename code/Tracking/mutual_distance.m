function dist = mutual_distance(pts1, pts2)
      
%%

    dist = zeros(size(pts1,1),size(pts2,1));
    for i=1:size(pts1,2)
        %size(bsxfun(@minus,pts1(:,i),pts2(:,i).').^2 )
        dist = dist + bsxfun(@minus,pts1(:,i),pts2(:,i).').^2 ;
    end
    
    dist = sqrt(dist);

	return;
    
    %%
end
