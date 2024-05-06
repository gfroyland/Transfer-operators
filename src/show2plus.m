%
% show 2D-Box-Coverings
%
% showdim2plus(b,p)
%
% p must be a nx1 vector (eg. invariant density)

function h=showdim2plus(b, p)

%pnorm=p/max(abs(p));
col=p';

x=[ b(:,1)-b(:,3) b(:,1)+b(:,3) b(:,1)+b(:,3) b(:,1)-b(:,3) ];
y=[ b(:,2)-b(:,4) b(:,2)-b(:,4) b(:,2)+b(:,4) b(:,2)+b(:,4) ];
h=patch(x',y',col);

shading flat