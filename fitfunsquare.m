function fun = fitfunsquare (x,tdata,b1)
i1 = x(1);
t1 = x(2);

fun = stepfunction(tdata-t1).*(b1 + i1*(tdata-t1).^2) + ...
      stepfunction(t1-tdata).*(b1);
