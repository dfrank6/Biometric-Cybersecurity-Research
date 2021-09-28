function scaled = sscaler(nonscaled)
therange = range(nonscaled);
    scaled = (nonscaled-min(nonscaled))/therange;
    scaled = (2*scaled)-1;
end