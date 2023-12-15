function [rho,tssptl,tsr] = fcn_spatial_null(ts,D,beta)
tsr = zscore(fcn_ampsurr_ind(ts));
W = exp(-beta*D);
tssptl = tsr*W;
rho = corr(tssptl);

function y = fcn_ampsurr_ind(x)
% this function is a modification of code provided by James Roberts and
% Michael Breakspear

%Creates surrogate multichannel data, by adding random numbers to phase component
%of all channel data, using amplitude adjusted algorithm
[r,c] = size(x);
if r < c
	x = x.';   % make each column a timeseries
end
[n,cc] = size(x);
m = 2^nextpow2(n);
yy=zeros(n,cc);
for i=1:cc    %create a gaussian timeseries with the same rank-order of x
   z=zeros(n,3); gs=sortrows(randn(n,1),1);
   z(:,1)=x(:,i); z(:,2)=[1:n]'; z=sortrows(z,1);
   z(:,3)=gs; z=sortrows(z,2); yy(:,i)=z(:,3);
end
% phsrnd=zeros(m,cc);
% phsrnd(2:m/2,1)=rand(m/2-1,1)*2*pi; phsrnd(m/2+2:m,1)=-phsrnd(m/2:-1:2,1);
% for i=2:cc
% %     phsrnd(:,i)=phsrnd(:,1);
%     phsrnd(2:m/2,i)=rand(m/2-1,1)*2*pi; phsrnd(m/2+2:m,i)=-phsrnd(m/2:-1:2,i);
% end
m = 2^nextpow2(n);
xx = fft(real(yy),m);
phsrnd=zeros(m,cc);
phsrnd(2:m/2,1)=rand(m/2-1,1)*2*pi; phsrnd(m/2+2:m,1)=-phsrnd(m/2:-1:2,1);

for i=2:cc
    phsrnd(2:m/2,i)=rand(m/2-1,1)*2*pi; phsrnd(m/2+2:m,i)=-phsrnd(m/2:-1:2,i);
end
xx = xx.*exp(phsrnd*sqrt(-1));
xx = ifft(xx,m);
xx = real(xx(1:n,:));
y=zeros(n,cc);
for i=1:cc    %reorder original timeseries to have the same rank-order of xx
   z=zeros(n,3); yst=sortrows(x(:,i));
   z(:,1)=xx(:,i); z(:,2)=[1:n]'; z=sortrows(z,1);
   z(:,3)=yst; z=sortrows(z,2); y(:,i)=z(:,3);
end
if r < c
   y = y.';
end
y=real(y);    %small imag. component created by rounding error
