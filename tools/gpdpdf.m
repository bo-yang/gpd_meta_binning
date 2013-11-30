function p=gpdpdf(x,l1,l2)

ln_x=0;
for k=1:x
	ln_x=ln_x+log(k);
end

p=exp(log(l1)+(x-1)*log(l1+x*l2)-(l1+x*l2)-ln_x);

%denorm=1;
%for k=1:x
%	denorm=denorm*k;
%end
%
%p=l1*((l1+x*l2)^(x-1))*exp(-(l1+x*l2))/denorm;
