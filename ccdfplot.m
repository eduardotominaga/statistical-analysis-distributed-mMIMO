function handleCCDF=ccdfplot(x,nBins)
    range=linspace(min(x),max(x),nBins);
    F_bar=zeros(1,length(range));
    for i=1:nBins
        F_bar(i)=sum(x>range(i))/length(x);
    end
    handleCCDF=plot(range,F_bar);
end