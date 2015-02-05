i = round(traj{idx}(:,1));
j = round(traj{idx}(:,2));

clf;

[N1,N2] = getImageDimensions(expe);
w=2*s+1;

tot = N;

d = round(tot/nR);

out = zeros(d*w,nR*w);
outseg = zeros(d*w,nR*w); 

dk = 1;

sel = find(ind(idx,:)>0);

for k=sel(1):1:sel(end)


    a   = imread( ['img/' getImageName(expe.colorNames{colorIndex},k)] );
    seg = imread( ['zStackedThreshCorrectedRefined/' n2s(k) '.png'] );
    
    [seli selj] = getNeiInd(i(k),j(k),s,N2,N1);

    %if(min(seli)>0 && min(selj>0) && max(seli) < N2 && max(selj) < N1 )

    sub_a = a(selj,seli);
    seg = seg(selj,seli);
    seg = bwmorph(seg,'remove');
    

    ind_out_i = (w*mod(dk-1,d)+1):(w*mod(dk-1,d)+w);
    ind_out_j = (w*(ceil(dk/d)-1) +1):(w*(ceil(dk/d)-1) +w);


    out(ind_out_i,ind_out_j,1) =  sub_a;
    outseg(ind_out_i,ind_out_j,1) =  seg;


    dk = dk+1;

end

p = 0.001;

out = out';
outseg = outseg';

out(out==0) = mean(out(out>0));

if(showSeg)
    out_ = repmat(out,[1,1,3]);
    tmp = out;
    tmp(outseg==1) = quantile(tmp(:),0.99);
    out_(:,:,1) = tmp;
    q = quantile(out_(:),1-p);
    out_(out_ > q) = q;
end

out = imnorm(out_);

%panel A
subplot(3,1,[1 2])

imagesc(out)

%%
set(gca,'XTick',[])
set(gca,'YTick',[w/2:w:nR*w])

ts = 0:(d*expe.dt):(expe.dt*expe.numberOfFrames);
ts = ts(1:length([w/2:w:nR*w]));
set(gca,'YTickLabel', round(ts*4)/4 )

title(expe.colorNames{colorIndex})
%%

caxis([ quantile(out(:),p) quantile(out(:),1-p) ])

hold on

dk = 1;
for k=sel(1):1:sel(end)
    
    ind_out_i = (w*mod(dk-1,d)+1):(w*mod(dk-1,d)+w);
    ind_out_j = (w*(ceil(dk/d)-1) +1):(w*(ceil(dk/d)-1) +w);
        
    if(divMatrix(idx,k)==1)
       plot(ind_out_i(1)+6, ind_out_j(1)+6,'*r','lineWidth',1)
    end

    dk = dk+1;
end

%panel B

subplot(3,1,3)

hold on

t = linspace(0,length(sel)*expe.dt,length(sel));
s =  signalToPlot(idx,sel);

divs = find(divMatrix(idx,sel));
peaks = find(peakMatrix(idx,sel));

for j=1:length(divs)
   plot([t(divs(j)) t(divs(j))],[min(s)-0.04 0.6*max(s(:))],'color','r') 
end

plot(t,s,'k')
axis([min(t) max(t) min(s) max(s)*1.1])
xlabel('time')

colormap gray