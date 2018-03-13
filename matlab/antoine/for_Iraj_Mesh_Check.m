clear

% epsilon=0.78;
% kappa=2;
% delta=0.35;

epsilon=0.32;
kappa=1.7;
delta=0.33;

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

gamma = 1.2 ;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Mesh domain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% N is the number of grid points in the theta variable
%
hconv = 0.04;

fd=@(p) p(:,1).^4/8+D(1)+D(2)*p(:,1).^2 ...
        +D(3)*(p(:,1).^4-4*p(:,1).^2.*p(:,2).^2)+D(4)*p(:,2)...
        +A*(1/2*p(:,1).^2.*log(p(:,1))-p(:,1).^4/8);%Construct implicit distance function - here it is easy, as we know that psi<0 inside the plasma
        
        [p,t]=distmesh2d(fd,@huniform,hconv,[1-epsilon,-kappa*epsilon;1+epsilon,kappa*epsilon],[]);%Compute mesh for this boundary
        hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  Plotting solutions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct the X-Y grid
[X,Y] = meshgrid(0:.005:1+epsilon+0.5,-gamma*kappa*epsilon-0.1:.005:kappa*epsilon+0.1);

% Construct pressure function (or flux function) Z
Z = X.^4/8+D(1)+D(2)*X.^2 ...
        +D(3)*(X.^4-4*X.^2.*Y.^2)+D(4)*Y...
        +A*(1/2*X.^2.*log(X)-X.^4/8);
% 
figure(1)
contour(X,Y,Z,[0 0],'LineWidth',2,'color','g')% Contour plot boundary of domain
axis equal