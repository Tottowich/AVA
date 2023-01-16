% Modeling a sliding 3D object
% Author: Theodor Jonsson
% Date 12/1/2023
% Testing using vectorized implementation instead of nested loops first.
% The strategy is to construct an adjecency matrix and use that together
% with similar method as in exercise 2 to simulate the system.
% The force function will now be more complicated but it will probabily be
% significantly faster and more general.
%
rng(6)
clear
close all
% ------- GIVEN PROPERTIES -------
Nx = 20; % Number of nodes in x-y-z net.
Ny = 20; % 
Nz = 2; % 
masses = 1; % Masses of the nodes not on the edge of the net.
ks = 500;
kd = 10;
g = 10;
dt = 2e-3;
n_dims = 3;
% Marble properties.
r_marble_max = 4; % Radius of the marble.
density_max = 50;
L = r_marble_max/4; % Evenly distributed particles => sqrt(2) on diagonal.
NM = 4;
density = density_max*rand(NM,1);
r_marble = max(r_marble_max*rand(NM,1),L); % Randomize some
M_marble = pi*r_marble.^3.*density_max/300;
% --------------------------------------
visualize = 1;
record = 0;
name = "Video/Lab6/Lab6MarbleCollisions";

start_x = 0;
start_y = 0;
start_z = 0;
NP = Nx*Ny*Nz; % Total number of particles in the spring grid.
% Time step set-up.
T = 3;
t_steps = T/dt;
ts = 0:dt:T-dt;

% ------- Set up the 2D object --------
% The net object is denoted by the matrix X
x = 0:L:(Nx-1)*L;
y = 0:L:(Ny-1)*L;
z = (Nz-1)*L:-L:0;
[xs,ys,zs] = meshgrid(x,y,z); % Starting at (0,0,0)
                              % goes to ((Nx-1)*L,(Nz-1)*L,(Nz-1)*L))
X_init = cat(4,xs,ys,zs);
X_init = reshape(X_init,[NP n_dims]); % Flatten the matrix.

% Initialize the marble position and velocity
X_marble_init = L+[max(x)-2*L,max(y)-2*L,0].*rand(NM,3);
X_marble_init(:,4) = r_marble;%r_marble + r_marble*rand(NM,1);
X_marble_init(:,3) = (Nz+1)*r_marble;
for i = 1:NM
    fprintf("Marble %d: r = %.3f [m], m = %.3f [kg]\n",i,r_marble(i),M_marble(i))
    fprintf("\t  xyz: [%.2f,%.2f,%.2f]\n",X_marble_init(i,1),X_marble_init(i,2),X_marble_init(i,3))
end
%X_marble_init(1,1:end-1) = X_marble_init(2,1:end-1);
%X_marble_init(1,3) = 17;
%
% X_init(:,3) = X_init(:,3)+start_z;
% X now has Shape (NP x n_dims)
% First Nc entries in X is the top row, Nc+1->2*Nc second row and so on.
V_init = zeros(NP,n_dims);
% Testing with some initial velocity
V_marble_init = zeros(NM,n_dims);
%V_init = V_init+v_init;

% -------------------------------------

% Now we can construct the adjecency matrix for the list of particles.
% Which we can use as a mask to remove the forces from non connected
% particles.
[A,small_diagonals,large_diagonals] = GridAdjacencyMatrix3D(Nx,Ny,Nz);
% Since we want to fix the edges of the net. We determine the indexes of
% these nodes.
fixed = repmat([1;zeros(Ny-2,1);1],Nx,1);
fixed(1:Ny) = 1;
fixed(end-Ny:end) = 1;
fixed = logical(repmat(fixed,Nz,1));
%figure(2)
%p = plot(graph(A),'k.-','XData',X_init(:,1),'YData',X_init(:,2),'ZData',X_init(:,3));
%highlight(p,fixed,"NodeColor",'r','MarkerSize',10)
%daspect([1,1,1]);
%grid on;
%legend("Fixed"+newline+"Loose")
% figure(2)
% scatter3(X_init(:,1),X_init(:,2),X_init(:,3),'b')
% % hold on
% % s = scatter3(X_init(fixed,1),X_init(fixed,2),X_init(fixed,3),'r')
% % title("Fixed and non fixed nodes of the net.")
% % legend("Non-fixed","Fixed")
% % grid on
% % set(s,'Color',)
% % hold off
%
% A = sparse(A);
% 'A' is the adjacency matrix, diagonals are only the diagonal springs.
% 

% ------- Set up the 2D spring --------
% Now construct 3 string matrices.
% These matrices are; spring constant per spring
%                     spring damping coefficients
%                     spring resting length
L_springs = zeros(NP,1,NP);
L_springs(A==1) = L; % Set every connection to L
% But diagonals are longer!
L_springs(small_diagonals==1) = sqrt(2)*L;
% Some diagonals are even larger ! 
L_springs(large_diagonals==1) = sqrt(3)*L;

%
ks_springs = zeros(NP,1,NP);
ks_springs(A==1) = ks;

kd_springs = zeros(NP,1,NP);
kd_springs(A==1) = kd;
%-------------------------------------
springs.L = L_springs;
springs.ks = ks_springs;
springs.kd = kd_springs;
% Masses of the particles.
M = ones(NP,1)*masses; % All particles have the same mass.

tic
[X,X_marble,V,V_marble] = LeapFrogMarbleBounce(X_init,V_init,X_marble_init,V_marble_init,fixed,springs,M,M_marble,g,t_steps,dt);
toc
% t3d = timeit(@() LeapFrogWithSurfaceCheck(X_init,V_init,F,M,circle_surface,t_steps,dt))
% Time to simulate using T=8,dt=2*10^-3,
% Nx=8,Ny=4,Nz=3,Nx_circles=10,Ny_circles=3.
% t3d = 2.6591 s
%
%
%%
figure(1)
if visualize
    VisualizeSpringMarble3D(X,A,X_marble,record,name)
end
%%
% disp("Saved to 'ExampleVideo.avi'")
figure(2)
[E,Ek,Es,Ep] = EnergyCalculationMarble(X,X_marble,V,V_marble,M,M_marble,g,ks_springs,L_springs);
PlotEnergies(E,Ek,Es,Ep,ts,kd)
figure(3);
% Track center of mass.
center_mass_vel = squeeze(sum(V.*M',2))./sum(M);
plot(ts,center_mass_vel(:,1))
grid on
title("Vx center of mass")




