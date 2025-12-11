function y = sinc(x)
idx=x==0;                                                              

y = sin(x)./x;                                                     
y(idx) = 1;   

end