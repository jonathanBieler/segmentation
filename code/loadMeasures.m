% Me = loadMeasures(N)
function Me = loadMeasures(N)

Me = cell(1,N);

for k=1:N
    
    name = ['Measures/' num2str(k) '.mat'];
    load(name)
    
    tmp = [];
    for i=1:length(Measurements)
          tmp = [tmp Measurements(i)];
    end
    
    Me{k} = tmp;        
end
