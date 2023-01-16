% Modeling a sliding microscopic 2D object
% Author: Theodor Jonsson
% Date 5/1/2023
% Testing using vectorized implementation instead of nested loops first.
% The strategy is to construct an adjecency matrix and use that together
% with similar method as in exercise 2 to simulate the system.
% The force function will now be more complicated but it will probabily be
% significantly faster and more general.
%
clear
close all
% ------- GIVEN PROPERTIES -------
Nc = 8; % Number of columns of the 2D object.
Nr = 4; % Rows of columns of the 2D object.
masses = 1; % All particles have mass 1.
ks = 100;
kd = 0;
g = 1;
dt = 2e-3;
L = 1; % Evenly distributed particles => sqrt(2) on diagonal.
n_dims = 2;
N_circles = 16; % Number of circles that make up the floor
dist_circle = 0.1; % [0<->1] describing how large portion of radius to seperate.
% --------------------------------------
visualize=0; % 1 for visualize, 0 not.
record = 0; % 1 recording, 0 not.
name = "Video/Lab3/Lab3GridSprings";

r_circle = L; % Randius of the circles which constructs the surface
start_x = 0;
start_y = r_circle;
vx_init = 7;
vy_init = 0;
NP = Nc*Nr; % Total number of particles in the spring grid.
% Time step set-up.
T = 1.2;
t_steps = T/dt;
ts = 0:dt:T-dt;


% ------- Set up the 2D object --------
% The object is denoted by the matrix X
x = meshgrid(0:L:(Nc-1)/L,0:L:(Nr-1)/L)+start_x;
y = flip(meshgrid(0:L:(Nr-1)/L,0:L:(Nc-1)/L)',1)+start_y;
X_init = cat(3,x',y');
X_init = reshape(X_init,[NP n_dims]); % Flatten the matrix.
% X now has Shape (NP x n_dims)

% First Nc entries in X is the top row, Nc+1->2*Nc second row and so on.
V_init = zeros(NP,n_dims);
% Testing with some initial velocity
% V_init(ceil(Nc/2)*3,:) = [5,5]; % Diagonal velocity of the top right particle.
V_init(:,1)=vx_init;
V_init(:,2)=vy_init;
% -------------------------------------

% Now we can construct the adjecency matrix for the list of particles.
% Which we can use as a mask to remove the forces from non connected
% particles.
[A,diagonals] = GridAdjacencyMatrix(Nr,Nc);
% 'A' is the adjacency matrix, diagonals are only the diagonal springs.
% 

% ------- Set up the 2D springs --------
% Now construct 3 string matrices.
% These matrices are; spring constant per spring
%                     spring damping coefficients
%                     spring resting length
L_springs = zeros(NP,1,NP);
L_springs(A==1) = L; % Set every connection to L

% But diagonals are longer!
L_springs(diagonals==1) = sqrt(2)*L;

ks_springs = zeros(NP,1,NP);
ks_springs(A==1) = ks;

kd_springs = zeros(NP,1,NP);
kd_springs(A==1) = kd;
%-------------------------------------

% Masses of the particles.
ms = ones(NP,1)*masses; % All particles have the same mass.
M = diag(ms); % Diagonal matrix.

reps = 1; % 400
%vx_inits = [3,5,7]; % Used for testing multple initial velocities (to find mu).
mus = zeros(reps,length(vx_inits));
for i = 1:length(vx_inits) % All the 
    vx_init = vx_inits(i);
    V_init(:,1)=vx_init;
    for r = 1:reps % We want multiple runs to increase accuracy in mu.
        circle_surface = BuildSurface(N_circles,r_circle,dist_circle,n_dims);
        
        
        % F = @(X,V) ForceFunction(X,V,A,ms,g,ks_springs,kd_springs,L_springs);
        F = @(X,V) ForceFunction(X,V,ms,g,ks_springs,kd_springs,L_springs);
        [X,V] = LeapFrogWithSurfaceCheck(X_init,V_init,F,M,circle_surface,t_steps,dt);
        %timeit(@() LeapFrogWithSurfaceCheck(X_init,V_init,F,M,circle_surface,t_steps,dt));
        if visualize
            figure(1)
            VisualizeSpringSystemWithSurface(X,A,circle_surface,record,name)% close
        end
        [E,Ek,Es,Ep] = EnergyCalculation(X,V,ms,g,ks_springs,L_springs);
        center_mass_vel = squeeze(sum(V.*ms',2))./sum(ms);
        if reps==1 % To avoid slow polling.
            figure(2)
            PlotEnergies(E,Ek,Es,Ep,ts,kd)
            
            % Track center of mass.
            figure(3);
            plot(ts,center_mass_vel(:,1))
            grid on
            title("Vx center of mass")
        end
        vx_diff = center_mass_vel(1,1)-center_mass_vel(end,1);
        % Average acceleration is therefore:
        ax_ave = vx_diff/T;
        % f_mu = f_x =ma_x=mu*mg*cos(theta), theta=0. =>
        mu = ax_ave/g;
        fprintf("\nFriction coefficient using initial x-velocity %.2f: mu = %.3f \n",vx_init,mu)
        mus(r,i) = mu;
    end
end
ave_mu = mean(mus);

% mean(mus) = 0.1665, v = 7.
function F_mat = ForceFunction(X,V,ms,g,ks,kd,L)
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
    F_mat = sum(F_tensor,3); % Sum along last channel.
                             % Last channel corresponds to each
                             % contribution from each spring
    % F_mat now have the correct shape of (NP x n_dims)
    % This has only taken into account the spring system.
    % Now add gravity!
    if size(F_mat,2)==2
        F_g = ms*g*[0 -1]; % (NPx1)x(1x2) => (NPx2)
    else
        F_g = ms*g*[0 0 -1]; % (NPx1)x(1x3) => (NPx3)
    end
    F_mat = F_mat+F_g;
end



