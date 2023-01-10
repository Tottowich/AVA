% Modeling a sliding 3D object
% Author: Theodor Jonsson
% Date 9/1/2023
% Testing using vectorized implementation instead of nested loops first.
% The strategy is to construct an adjecency matrix and use that together
% with similar method as in exercise 2 to simulate the system.
% The force function will now be more complicated but it will probabily be
% significantly faster and more general.
%
clear
close all
% ------- GIVEN PROPERTIES -------
Nx = 9; % Number of nodes in x-y-z net.
Ny = 6; % 
Nz = 1; % 
masses = 0.01; % Masses of the nodes not on the edge of the net.
ks = 1000;
kd = 50;
g = 2;
dt = 2e-3;
L = 1; % Evenly distributed particles => sqrt(2) on diagonal.
n_dims = 3;
% Marble.
r_marble = 1; % Radius of the marble.
M_marble = 0.1; % Mass of the marble.
x_marble = Nx/2;
y_marble = Ny/2;
z_marble = 5;
% --------------------------------------
visualize = 1;
record = 0;
name = "Video/Lab6/Lab6Marble";

start_x = 0;
start_y = 0;
start_z = 0;
NP = Nx*Ny*Nz; % Total number of particles in the spring grid.
% Time step set-up.
T = 2;
t_steps = T/dt;
ts = 0:dt:T-dt;

% ------- Set up the 2D object --------
% The net object is denoted by the matrix X
x = 0:L:(Nx-1)/L;
y = 0:L:(Ny-1)/L;
z = (Nz-1)/L:-L:0;
[xs,ys,zs] = meshgrid(x,y,z); % Starting at (0,0,0)
                              % goes to ((Nx-1)*L,(Nz-1)*L,(Nz-1)*L))
X_init = cat(4,xs,ys,zs);
X_init = reshape(X_init,[NP n_dims]); % Flatten the matrix.
% Initialize the marble position and velocity
X_marble_init = [x_marble y_marble z_marble r_marble];
% X_init(:,3) = X_init(:,3)+start_z;
% X now has Shape (NP x n_dims)

% First Nc entries in X is the top row, Nc+1->2*Nc second row and so on.
V_init = zeros(NP,n_dims);
% Testing with some initial velocity
V_marble_init = zeros(1,n_dims);
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
fixed = repmat(fixed,Nz,1);
figure(1)
scatter3(X_init(:,1),X_init(:,2),X_init(:,3),'b')
hold on
scatter3(X_init(logical(fixed),1),X_init(logical(fixed),2),X_init(logical(fixed),3),'r')
title("Fixed and non fixed nodes of the net.")
legend("Non-fixed","Fixed")
grid on
hold off
%%
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
ms = ones(NP,1)*masses; % All particles have the same mass.
M = diag(ms);

% F = @(X,V) ForceFunction(X,V,ms,g,ks_springs,kd_springs,L_springs);
% [X,V] = LeapFrogWithSurfaceCheck(X_init,V_init,F,M,circle_surface,t_steps,dt);

[X,X_marble,V,V_marble] = LeapFrogMarbleBounce(X_init,V_init,X_marble_init,V_marble_init,springs,M,M_marble,t_steps,dt);
% t3d = timeit(@() LeapFrogWithSurfaceCheck(X_init,V_init,F,M,circle_surface,t_steps,dt))
% Time to simulate using T=8,dt=2*10^-3,
% Nx=8,Ny=4,Nz=3,Nx_circles=10,Ny_circles=3.
% t3d = 2.6591 s
%
%%
figure(1)
if visualize
    VisualizeSpringSystemWithSurface3D(X,A,circle_surface,record,name)
end
%
% disp("Saved to 'ExampleVideo.avi'")
figure(2)
[E,Ek,Es,Ep] = EnergyCalculation(X,V,ms,g,ks_springs,L_springs);
PlotEnergies(E,Ek,Es,Ep,ts,kd)
figure(3);
% Track center of mass.
center_mass_vel = squeeze(sum(V.*ms',2))./sum(ms);
plot(ts,center_mass_vel(:,1))
grid on
title("Vx center of mass")

function [F_net,F_marble] = ForceFunction(X,V,X_marble,V_marble,ms,g,ks,kd,L)
    % This is the force function of the current lab exercise.
    %
    % INPUT
    %   X - (mat) Positions of each particle at current time step.
    %             Shape: (NP x n_dims). Dimensions in order (x,y,z)
    %
    %   V - (mat) Velocities of each particle at current time step.
    %             Shape: (NP x n_dims)
    %   ms - (vec) The masses of each particle.
    %             Shape: (NP x 1)
    %   g - (float) gravitational acceleration constant.
    %
    %   ks - (mat/float) either matrix of shape (NP x 1 x NP) or float.
    %                    If matrix then ks(i,j) indicates coefficient of
    %                    the spring between particle i and particle j. The
    %                    matrix must be symmetric to make sense.
    %   kd - (mat/float) either matrix of shape (NP x 1 x NP) or float.
    %                    If matrix then ks(i,j) indicates damping coefficient 
    %                    of the spring between particle i and particle j. The
    %                    matrix must be symmetric to make sense.
    %
    %   L - (mat/float)  either matrix of shape (NP x NP) or float.
    %                    If matrix then ks(i,j) sym indicates length of
    %                    the spring between particle i and particle j at
    %                    rest.The matrix must be symmetric to make sense.
    %

    % Create a distance and relative velocity tensors of shape (NP x n_dims x NP)
    R = X - permute(X, [3 2 1]); % Relative positions
    V_rels = V-permute(V, [3 2 1]); % Relative velocities
    rs = vecnorm(R,2,2); % Euclidian norm on the second channel to 
                         % get the length of each spring.
                         % This can then be used to construct r_bars.
    r_bars = R./rs; % Shape - (NP x n_dims x NP), will be anti symmetric.
    % Replace NaN with zeros. 
    r_bars(isnan(r_bars))=0;
    % We now want to compute the forces according to the formula (5) of the
    % lab instructions.

    % SPRING
    F_spring = ks.*(rs-L); % Shape (NP x 1 x NP), one spring from each 
                           % particle to another. ks couble be a matrix of
                           % shape (NP x 1 x NP) with different strengths
                           % for each spring.
    % DAMPING
    F_damping = kd.*dot(V_rels,R,2)./rs;
    F_damping(isnan(F_damping))=0;
    % Multiply with the unit vectors of each individual spring.
    F_tensor = -(F_spring+F_damping).*r_bars; % (NP x n_dims x NP)
    % The Entries F(i,:,i) should be zero since this corresponds to the
    % force asserted on particle i on particle i. Which is always zero.
    F_net = sum(F_tensor,3); % Sum along last channel.
                             % Last channel corresponds to each
                             % contribution from each spring
    % F_mat now have the correct shape of (NP x n_dims)
    
    % We must now compute the forces between the bouncing marble and the
    % net.
    centers = X_marble(:,1:end-1);
    radii = X_marble(:,end);
    r = -(centers - permute(X, [3 2 1])); % Shape (N_
    pos_diff = vecnorm(r,2,2);
    r_bars = r./pos_diff;
    inter = squeeze(pos_diff)<=radii; % Logical intersection matrix.
    [inter_circles,inter_particles] = find(inter==1);
    [inter_particles,id] = unique(inter_particles,'first');
    % This has only taken into account the spring system.
    % Now add gravity!
    if size(F_mat,2)==2
        F_g = ms*g*[0 -1]; % (NPx1)x(1x2) => (NPx2)
    else
        F_g = ms*g*[0 0 -1]; % (NPx1)x(1x3) => (NPx3)
    end
    F_mat = F_mat+F_g;
end



