function dum=oceanplot(v,landpoints,lon,lat),

dum=0;
maskedv=NaN*zeros(numel(landpoints),1);
maskedv(landpoints==0)=v;
maskedv=reshape(maskedv,numel(lon),numel(lat));
out=real(maskedv)';
imagesc(lon,lat,out);axis xy
%imagesc(lon,lat,sign(real(out)).*sqrt(abs(real(out))));axis xy

