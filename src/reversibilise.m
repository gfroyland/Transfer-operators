function [R,Phat]=reversibilise(P,p),

%creates the reverse-time transition matrix Phat and
%the symmetrised time-reversible transition matrix R.

p=p/sum(p);
pmat=diag(sparse(p));
pmatinv=diag(sparse(1./p));

Phat=pmatinv*P'*pmat;

R=(P+Phat)/2;

