%%show simulation
load buoytransitmatrix_vSebille_England_F(2012).mat

unif=ones(33673,1);

for i=0:1000
clf
dum=oceanplot((unif).^(1/3),landpoints,lon,lat);
title([num2str(i) ' years'])
%pause
drawnow
unif=(unif'*P)';
end


