function [ P ] = gpdv(X,lambda1,lambda2)
%gpdv: calculate generalized Poisson distribution for given vector.
%   X - a vector containing unique distinct k-mer occurrences.

for i=1:numel(X)
    P(i)=gpd(X(i),lambda1,lambda2);
end

end