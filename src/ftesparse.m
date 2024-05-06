function fte=ftesparse(P,k),

Pk=P^k;
[i,j,nonzeroPk]=find(Pk);
fte=-(1/k)*sum(sparse(i,j,nonzeroPk.*log(nonzeroPk))');
