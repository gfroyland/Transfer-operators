function abstime=absorbtime(P,targetind),
%P is a square stochastic matrix and targetind is a vector of indices that we want
%to hit.  abstime has the same length as the dimension of P.

%make restricted Pind matrix
restrind=setdiff(1:length(P),targetind);
Prestr=P(restrind,restrind);

%compute absorption times
Prestrmat=Prestr-speye(length(restrind));
RHStime=-ones(length(restrind),1);
abstime=zeros(length(P),1);
abstime(restrind)=Prestrmat\RHStime;
abstime(targetind)=zeros(length(targetind),1);


