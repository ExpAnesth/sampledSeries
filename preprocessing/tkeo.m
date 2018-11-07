function e=tkeo(d)
% ** function e=tkeo(d)
% computes the Teager-Kaiser energy of data array d (time series in columns).
% For discrete signals, Teager-Kaiser energy is
% x(n) = x^2(n) - x(n+1)*x(n-1)


% check time series
[n1,n2,n3]=size(d);
if n1<=1 
  if xor(n2>1,n3>1)
    warning('time series of d must be in columns - reshaping')
    d=d(:);
    [n1,n2,n3]=size(d);
  else    
    error('time series of d must be in columns')
  end
end


e=d.^2 + circshift(d,1,1).*circshift(d,-1,1);
