mstd = [];
msum = [];

i = round(traj{idx}(:,1));
j = round(traj{idx}(:,2));


clf;


if( gaussianFilterSize >  0)
    H = fspecial('gaussian',500,gaussianFilterSize);
end

[N1,N2] = getImageDimensions(expe);


w=2*s+1;
ws = 2*s*binsize+1;

if( useFullSizeImages)
    w = 2*s*binsize+1;
end

tot = NtoPlot;

d = round(tot/nR);

if( useFullSizeImages)
    out = zeros(d*ws,nR*ws);
    outseg = zeros(d*ws,nR*ws); 
else

    out = zeros(d*w,nR*w);
    outseg = zeros(d*w,nR*w); 
end

dk = 1;

sel = find(ind(idx,:)>0);

for k=start:1: min(sel(end), start+NtoPlot-1)

    disp(k/sel(end))

    if( useFullSizeImages )
        a = imread(['fullSizeimg/' getImageName(expe.colorNames{colorIndex},k)]);
    else
        a = imread(['img/' getImageName(expe.colorNames{colorIndex},k)]);
    end
    
    
    if(exist(['zStackedThreshCorrectedRefined/' n2s(k) '.png'] ,'file'))
        seg = imread( ['zStackedThreshCorrectedRefined/' n2s(k) '.png'] );
    else
        seg = imread( ['zStackedThreshCorrected/' n2s(k) '.png'] );
    end
    
    [seli selj] = getNeiInd(i(k),j(k),s,N1,N2);
    
    if( useFullSizeImages)
        [selif seljf] = getNeiInd(1+binsize*(i(k)-1),1+binsize*(j(k)-1),binsize*s,binsize*N1,binsize*N2);
        sub_a = a(seljf,selif);
    else
        sub_a = a(selj,seli);
    end

    %if(min(seli)>0 && min(selj>0) && max(seli) < N2 && max(selj) < N1 )

    
    seg = seg(selj,seli);
    
    if( useFullSizeImages )
        seg = imresize(seg, size(sub_a));
    end
    
    
    seg = bwmorph(seg,'remove');
    

    if( useFullSizeImages )
        
               
        ind_out_i = (ws*mod(dk-1,d)+1):(ws*mod(dk-1,d)+ws);
        ind_out_j = (ws*(ceil(dk/d)-1) +1):(ws*(ceil(dk/d)-1) +ws);
        
    else
        
        ind_out_i = (w*mod(dk-1,d)+1):(w*mod(dk-1,d)+w);
        ind_out_j = (w*(ceil(dk/d)-1) +1):(w*(ceil(dk/d)-1) +w);
    
    end
    
    sub_a = double(sub_a);

    if( gaussianFilterSize >  0)        
        sub_a = sub_a - imfilter(sub_a,H,'replicate'); 
    end
    
    mstd  = [mstd; std((sub_a(:)))];
    msum =  [msum; sum((sub_a(:)))];
    
    if(doNormalize)
        sub_a = qnorm(sub_a,0.01);
    end

    out(ind_out_i,ind_out_j,1) =  sub_a;
    outseg(ind_out_i,ind_out_j,1) =  seg;


    dk = dk+1;
    
    if(doDraw)
        clf
        imagesc(out)
        drawnow
    end

end

p = 0.001;

out = out';
outseg = outseg';

out(out==0) = min(out(out>0));

if(showSeg)
    out_ = repmat(out,[1,1,3]);
    tmp = out;
    tmp(outseg==1) = quantile(tmp(:),0.99);
    out_(:,:,1) = tmp;
    q = quantile(out_(:),1-p);
    out_(out_ > q) = q;
else
    out_ = out;
end

out = imnorm(out_);

%% panel A

clf
subplot(3,1,[1 2])

imagesc(out)

%
set(gca,'XTick',[])
set(gca,'YTick',[w/2:w:nR*w])

ts = 0:(d*expe.dt):(expe.dt*expe.numberOfFrames);
ts = ts(1:length([w/2:w:nR*w]));
set(gca,'YTickLabel', round(ts*4)/4 )

title(expe.colorNames{colorIndex})

ylabel('time [h]')
caxis([ quantile(out(:),p) quantile(out(:),1-p) ])

%caxis([ 0 1 ])

%%


hold on

dk = 1;
for k=start:1: min(sel(end), start+NtoPlot-1)
    
        
    if(divMatrix(idx,k)==1)
        
       ind_out_i = (w*mod(dk-1,d)+1):(w*mod(dk-1,d)+w);
       ind_out_j = (w*(ceil(dk/d)-1) +1):(w*(ceil(dk/d)-1) +w);
        
       plot(ind_out_i(1)+40, ind_out_j(1)+40,'*r','lineWidth',1)
    end

    if(peakMatrix(idx,k,colorIndex)==1)
        
       ind_out_i = (w*mod(dk-1,d)+1):(w*mod(dk-1,d)+w);
       ind_out_j = (w*(ceil(dk/d)-1) +1):(w*(ceil(dk/d)-1) +w);
        
       plot(ind_out_i(1)+40, ind_out_j(1)+40,'*r','lineWidth',1)
    end
    
    dk = dk+1;
end

%% panel B

subplot(3,1,3)

hold on
sel = find(ind(idx,:)>0);

t = linspace(0,length(sel)*expe.dt,length(sel));
%s =  imnorm(signalToPlot(idx,sel));

s =  (signalToPlot(idx,sel));

divs = find(divMatrix(idx,sel));
peaks = find(peakMatrix(idx,sel,colorIndex));

for j=1:length(divs)
   plot([t(divs(j)) t(divs(j))],[min(s)-0.04 0.6*max(s(:))],'color','r') 
end

for j=1:length(peaks)
   plot([t(peaks(j)) t(peaks(j))],[min(s)-0.04 0.6*max(s(:))],'color','b') 
end

plot(t,s,'k')
axis([min(t) max(t) min(s) max(s)*1.1])
xlabel('time [h]')

colormap gray