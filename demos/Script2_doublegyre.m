%note that this tutorial is not the highest resolution, 
%nor the most cutting-edge results.
%the aim is to present familar dynamical systems
%at a resolution that can be done very quickly.
%I have concentrated on transfer operator topics
%related to fluid flow.

%on the set-oriented side there are gaio algorithms
%for computing invariant manifolds, attracting sets,
%invariant sets, chain recurrent sets...
%on the transfer operator side, there are theory and algorithms
%for computing invariant measures, estimating mixing rates
%and extensions of these ideas to open dynamical sytems,
%random dynamical systems, and infinite-time nonautonomous systems.

%IF YOU DO NOT HAVE GAIO INSTALLED, YOU CAN INSTEAD
%load depth13.mat
%AND SKIP AHEAD TO BEGIN AT THE ``COMPUTE EIGENVECTORS...'' CODE BLOCK. 

addpath(genpath('D:/gaio/GAIO-3.0/GAIO-3.0'))
addpath ../src

% Double Gyre demo 
%

% Double Gyre System
A=0.25; delta=0.25; omega=2*pi; Tstart=0;
v = @(x,T) [-pi*A*sin(pi*(delta*sin(omega*T)*x(:,1).^2+(1-2*delta*sin(omega*T))*x(:,1))).*cos(pi*x(:,2)) pi*A*cos(pi*(delta*sin(omega*T)*x(:,1).^2+(1-2*delta*sin(omega*T))*x(:,1))).*sin(pi*x(:,2)).*(delta*sin(omega*T)*2*x(:,1)+(1-2*delta*sin(omega*T)))];

% map is given by the time integrated vector field
% v is the vector field above
% x is the spatial coordinate
% h is the timestep is h
% n is the number of steps
% Tstart is the initial timenumber of steps to form the
h = 0.01; n = 100; f = @(x) rk4t(v,x,h,n,Tstart); 

% sample points uniformly distributed over [-1,1]^2 (template box)
% in this example, n^2=20^2=400 points are generated per box
n = 20; x = linspace(-1+1/(2*n),1-1/(2*n),n)'; [XX,YY] = meshgrid(x,x);
X = [ XX(:) YY(:) ];

% the tree
% in this example, the initial box is [0,2]\times [0 1]
c = [1 0.5]; r = [1 0.5]; t = Tree(c,r);


%% construct full subdivison

%boxes are subdivided by a factor of two,
%cyclically around the two dimensions, 13 times 
%volume of each box is 2^{-13} of the volume of the original box.
sd = 8; depth = 13;
for i=1:depth,
    t.set_flags('all', sd);
    t.subdivide;
end

b=t.boxes(-1)';

%% compute transition matrix

tic

P = tpgraph(t, f, X, depth)';

toc

%% compute eigenvectors (to get time-asymptotic almost-invariant sets)

[U,V]=eigs(P',3);
p=U(:,1);
p=p/sum(p);

figure;show2plus(b,U(:,2))

%threshold.
[youtmax,ioutmax]=threshold_AI(b,P,U(:,2),2,p)

%third eigenvector gives both eggs vs. chaotic sea
figure;show2plus(b,U(:,3))
%can also threshold with U(:,3) to get both ``eggs''.
[youtmax,ioutmax]=threshold_AI(b,P,U(:,3),2,p)

%simplex plot
plot(U(:,2),U(:,3),'.','markersize',1)

%% compute finite-time almost-invariant sets

%create time-symmetrised transition matrix R
[R,Phat]=reversibilise(P,p);

%compute right eigenvectors of R
[ur,vr]=eigs(R,3);

%and so on...

%% compute singular vectors (t=1)  (can also use t=2.5, t=\pi, etc...)

%comment that here the initial and final domains are the same here, but 
%they need not be.  P can be rectangular.

[L,lv,rv,S]=coherent_vectors(P,ones(8192,1)/8192,3);
figure;
subplot(2,1,1);show2plus(b,lv(:,2))
subplot(2,1,2);show2plus(b,rv(:,2))

%the weighting is important;  you'll see something, but maybe not what you
%should.
[naive_lv,naive_s,naive_rv]=svds(P,3);
figure;
subplot(2,1,1);show2plus(b,naive_lv(:,2))
subplot(2,1,2);show2plus(b,naive_rv(:,2))

%threshold
[youtmax,ioutmax]=threshold_coherent(b,L,lv(:,2),rv(:,2),3,ones(8192,1)/8192)

%compute singular vectors (t=2) (cheap way to get t=2).
[L,lv,rv,S]=coherent_vectors(P^2,ones(8192,1)/8192,3);
figure;
subplot(2,1,1);show2plus(b,lv(:,2))
subplot(2,1,2);show2plus(b,rv(:,2))

%threshold
[youtmax,ioutmax]=threshold_coherent(b,L,lv(:,2),rv(:,2),3,ones(8192,1)/8192)

%the third singular vector (t=2)
figure;
subplot(2,1,1);show2plus(b,lv(:,3))
subplot(2,1,2);show2plus(b,rv(:,3))

%threshold
[youtmax,ioutmax,iout1,iout2]=threshold_coherent(b,L,lv(:,3),rv(:,3),3,ones(8192,1)/8192);


%% compute finite-time entropy

%note that initial and final domains don't have to be the same, ie. P can
%be rectangular.

load depth15

fte=ftesparse(P,5)';
figure; show2plus(b,fte);

%show evolution of high/low-stretching box
hold on;
plot(b(17090,1),b(17090,2),'*')
plot(b(4096,1),b(4096,2),'w*')

ballvec1=zeros(32768,1)';
ballvec1(17090)=1;

ballvec2=zeros(32768,1)';
ballvec2(4096)=1;

figure
for i=0:5,
clf;subplot(2,1,1);show2plus(b,(ballvec1*P^i)');
subplot(2,1,2);show2plus(b,(ballvec2*P^i)');pause
end

%do backward time
Pback=stochasticise(P');
fte=ftesparse(Pback,5)';
figure; show2plus(b,fte);

%SEE FIGURES IN FTE PAPER, INCLUDING SDE-VERSION OF DOUBLE GYRE

%% compute absorption times

%target is RHS of double-gyre system
targetind=find(b(:,1)>1);
abstime=absorbtime(P,targetind),
figure; plot(abstime);
figure;show2plus(b,abstime);

%target is upper half of double-gyre system
targetind=find(b(:,2)>.5)
abstime=absorbtime(P,targetind),
figure; plot(abstime);
figure;show2plus(b,log(abstime));

