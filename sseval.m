function sse = sseval(x,tdata,ydata,baseline)

sse = sum((ydata - fitfunsquare(x,tdata,baseline)).^2);