function maxp = getDivisionsFromTrace(tmp,q,th)

    tmp = fillTrace(tmp);

    tmp = imnorm(tmp);

    f = [1 1 -5 1 1];

    pe = conv(tmp,f,'same');

    pe(pe<0)=0;

    %[maxp minp] = peakdet(pe,0.11);
    maxp = [];
    detected = find(pe> quantile(pe,q) & pe > th);

    if(~isempty(detected))
        maxp(:,1) = detected;
        maxp(:,2) = pe(maxp(:,1));

        %remove close by peaks
        [v idx] = sort(maxp(:,1),'descend');
        for i=1:idx

            isClose = abs( maxp(idx(i),1) - maxp(:,1)) < 5;
            smaller = maxp(idx(i),2) - maxp(:,2) > 0;

            maxp( isClose & smaller ,1) = 0;     
        end

        nmaxp = [];
        nmaxp(:,1) = maxp(maxp(:,1)>0,1);
        nmaxp(:,2) = maxp(maxp(:,1)>0,2);
        maxp = nmaxp;

    end




