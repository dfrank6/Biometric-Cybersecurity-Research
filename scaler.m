function scaled = scaler(nonscaled)
therange = range(nonscaled,2);
scaled = [];
for i = 1:5
    x = (nonscaled(i,:)-min(nonscaled(i,:)))/therange(i);
    x = (2*x)-1;
    scaled = [scaled;x];
end
end