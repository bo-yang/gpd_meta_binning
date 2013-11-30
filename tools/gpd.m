function [ p ] = gpd(x,lambda1,lambda2)
% gpd: PDF of generalized Poisson distribution.
%	gpd=lambda1*(lambda1+x*lambda2)^(x-1)*exp(-(lambda1+x*lambda2))

logxperm=0;
for i=1:x
	logxperm=logxperm+log(i);
end

logprob=log(lambda1)+(x-1)*log(lambda1+x*lambda2)-(lambda1+x*lambda2)-logxperm;

p=exp(logprob);

end
