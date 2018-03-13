clear

epsilon=0.32;
kappa=1.7;
delta=0.33;

hconv = 0.01:0.01:0.1;
pp = length(hconv);
error = speye(pp,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               Choose A to specify beta=kinetic pressure/magnetic pressure
%               value
%
%               A = 1: "Force-free" equilibrium, i.e. no pressure or
%               constant pressure
%
%               A = 0: Low beta case: plasma neither paramagnetic nor
%               diamagnetic
%
%               A < 0: Towards higher beta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Define up-down asymmetry factor gamma
%
%           The flux condition psi=0 at the bottom point with y coordinate
%           y=gamma*kappa*epsilon is imposed
%
%           If gamma = 1, equilibrium is updown symmetric
%
%           If gamma > 1, equilibrium is more elongated at the bottom
%
%           If gamma < 1, equilibrium is more elongated at the top
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 1.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute geometric coefficients d_1, d_2, d_3, d_4 of the homogeneous
% solution following a reduced version of the procedure proposed in A.J.
% Cerfon and J.P. Freidberg, ``One size fits all” analytic solutions to 
% the Grad–Shafranov equation, Physics of Plasmas 17, 032502 (2010)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix entries for the homogeneous solution

M = [1 (1+epsilon)^2 (1+epsilon)^4 0; 
    1 (1-epsilon)^2 (1-epsilon)^4 0;...
    1 (1-delta*epsilon)^2 ...
    (1-delta*epsilon)^4-4*(1-delta*epsilon)^2*kappa^2*epsilon^2 ...
    kappa*epsilon;
    1 (1-delta*epsilon)^2 ...
    (1-delta*epsilon)^4-4*(1-delta*epsilon)^2*gamma^2*kappa^2*epsilon^2 ...
    -gamma*kappa*epsilon];

%Matrix entries for the particular solutions

B=-[(1+epsilon)^4/8+A*(1/2*(1+epsilon)^2*log(1+epsilon)-(1+epsilon)^4/8);
    (1-epsilon)^4/8+A*(1/2*(1-epsilon)^2*log(1-epsilon)-(1-epsilon)^4/8);
    (1-delta*epsilon)^4/8+A*(1/2*(1-delta*epsilon)^2*log(1-delta*epsilon)...
    -(1-delta*epsilon)^4/8);
    (1-delta*epsilon)^4/8+A*(1/2*(1-delta*epsilon)^2*log(1-delta*epsilon)...
    -(1-delta*epsilon)^4/8)];

% Compute coefficients d_1, d_2, d_3, and d_4

D=M\B;

for i = 1:pp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Mesh domain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% N is the number of grid points in the theta variable
%
params = [D',A,gamma]

fd=@(p) p(:,1).^4/8+D(1)+D(2)*p(:,1).^2 ...
        +D(3)*(p(:,1).^4-4*p(:,1).^2.*p(:,2).^2)+D(4)*p(:,2)...
        +A*(1/2*p(:,1).^2.*log(p(:,1))-p(:,1).^4/8);%Construct implicit distance function - here it is easy, as we know that psi<0 inside the plasma
        
        [p,t]=distmesh2d(fd,@huniform,hconv(i),[1-epsilon,-kappa*epsilon;1+epsilon,kappa*epsilon],[]);%Compute mesh for this boundary
        hold on

b=unique(boundedges(p,t));

n=int2str(length(b));
aie = int2str(i)
s1 = '../../infiles/';
mkdir(strcat(s1,aie))
sb = strcat(s1,aie,'/b.txt');
sp = strcat(s1,aie,'/p.txt');
st = strcat(s1,aie,'/t.txt');
sh = strcat(s1,aie,'/h.txt');
sr = strcat(s1,aie,'/params.txt');

size = hconv(i)
save(sb,'b','-ascii');
save(sp,'p','-ascii');
save(st,'t','-ascii');
save(sh,'size','-ascii');
save(sr,'params','-ascii');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Solve problem numerically with FEM method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

%% [K,F] = assemble(p,t) % K and F for any mesh of triangles: linear phi's
%N=size(p,1);T=size(t,1); % number of nodes, number of triangles
%% p lists x,y coordinates of N nodes, t lists triangles by 3 node numbers
%K=sparse(N,N); % zero matrix in sparse format: zeros(N) would be "dense"
%F=zeros(N,1); % load vector F to hold integrals of phi's times load f(x,y)
% 
%for e=1:T  % integration over one triangular element at a time
%  nodes=t(e,:); % row of t = node numbers of the 3 corners of triangle e
%  Pe=[ones(3,1),p(nodes,:)]; % 3 by 3 matrix with rows=[1 xcorner ycorner]
%  Area=abs(det(Pe))/2; % area of triangle e = half of parallelogram area
%  C=inv(Pe); % columns of C are coeffs in a+bx+cy to give phi=1,0,0 at nodes
%  % now compute 3 by 3 Ke and 3 by 1 Fe for element e
%  grad=C(2:3,:);Ke=Area*grad'*grad/mean(p(nodes,1)); % element matrix from slopes b,c in grad
%  Fe=-Area/3*((1-A)*mean(p(nodes,1))+A/mean(p(nodes,1))); % integral of phi over triangle is volume of pyramid: f(x,y)=CO*R
%  % multiply Fe by f at centroid for load f(x,y): one-point quadrature!
%  % centroid would be mean(p(nodes,:)) = average of 3 node coordinates
%  K(nodes,nodes)=K(nodes,nodes)+Ke; % add Ke to 9 entries of global K
%  F(nodes)=F(nodes)+Fe; % add Fe to 3 components of load vector F
%end   % all T element matrices and vectors now assembled into K and F
% 
%% [Kb,Fb] = dirichlet(K,F,b) % assembled K was singular! K*ones(N,1)=0
%% Implement Dirichlet boundary conditions U(b)=0 at nodes in list b
%K(b,:)=0; K(:,b)=0; F(b)=0; % put zeros in boundary rows/columns of K and F
%K(b,b)=speye(length(b),length(b)); % put I into boundary submatrix of K
%Kb=K; Fb=F; % Stiffness matrix Kb (sparse format) and load vector Fb
% 
%% Solving for the vector U will produce U(b)=0 at boundary nodes
%U=Kb\Fb;  % The FEM approximation is U_1 phi_1 + ... + U_N phi_N
%
%% Exact solution for comparison
%W = p(:,1).^4/8+D(1)+D(2)*p(:,1).^2 ...
%        +D(3)*(p(:,1).^4-4*p(:,1).^2.*p(:,2).^2)+D(4)*p(:,2)...
%        +A*(1/2*p(:,1).^2.*log(p(:,1))-p(:,1).^4/8);
%    
%error(i)=max(abs(U-W));
%end        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%                  Convergence plot
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%figure(2)
%loglog(hconv,error)
%hold on
%loglog(hconv,hconv.^2*error(end)/hconv(end)^2)
