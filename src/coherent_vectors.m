function [L,leftvecs,rightvecs,S]=coherent_vectors(P,initp,numvecs),

q=(initp'*P)';
ScaleLeft=diag(sparse(initp.^(1/2)));
ScaleRight=diag(sparse(q.^(1/2)));
L=ScaleLeft*P*inv(ScaleRight);

[U,S,V]=svds(L,numvecs);

leftvecs=(U'*inv(ScaleLeft))';
rightvecs=(V'*inv(ScaleRight))';
