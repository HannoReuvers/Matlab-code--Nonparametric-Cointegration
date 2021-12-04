if strcmp(KernelName, 'Epanechnikov')
    MyKernel = @(x) 0.75*(abs(x)<=1).*(1-x.^2);
    muK1 = 1;       % integral over kernel
    muK2 = 0.6;     % integral over kernel squared
end