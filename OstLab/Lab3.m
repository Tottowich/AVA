% Modeling a sliding microscopic 2D object
% Author: Theodor Jonsson
% Date 5/1/2023
% Testing using vectorized implementation instead of nested loops first.
% The strategy is to construct an adjecency matrix and use that together
% with similar method as in exercise 2 to simulate the system.
% The force function will now be more complicated but it will probabily be
% significantly faster and more general.
%

% ------- GIVEN PROPERTIES -------
Nc = 8; % Number of columns of the 2D object.
Nr = 6; % Rows of columns of the 2D object.
masses = 1; % All particles have mass 1.
ks = 100;
kd = 2;
g = 1;
dt = 2e-3;
L = 1; % Evenly distributed particles => sqrt(2) on diagonal.
% --------------------------------
NP = Nc*Nr;
n_dims = 2;
% ------- Set up the 2D object --------
% The object is denoted by the matrix X
x = meshgrid(0:L:(Nc-1)/L,0:L:(Nr-1)/L);
y = flip(meshgrid(0:L:(Nr-1)/L,0:L:(Nc-1)/L)',1);
X = cat(3,x',y');
X = reshape(X,[NP n_dims]); % Flatten the matrix.
% X now has Shape (NP x n_dims)
% First Nc entries in X is the top row, Nc+1->2*Nc second row and so on.
V = zeros(NP,n_dims);
% -------------------------------------

% Now we can construct the adjecency matrix for the list of particles.
% Which we can use as a mask to remove the forces from non connected
% particles.

% Figure to check connection.
figure(4)
[A,diagonals] = GridAdjacencyMatrix(Nr,Nc);
gplot(A,X,'b-o')
grid on
axis padded


% ------- Set up the 2D spring --------
% Now construct 3 string matrices.
% These matrices are; spring constant per spring
%                     spring damping coefficients
%                     spring resting length

L_springs = zeros(NP,1,NP);
L_springs(A==1) = L;
L_springs(diagonals==1) = sqrt(2)*L; % NOTE: Assuming equaly distributed.

ks_springs = zeros(NP,1,NP);
ks_springs(A==1) = ks;

kd_springs = zeros(NP,1,NP);
kd_springs(A==1) = kd;
%-------------------------------------

% Masses of the particles.
ms = ones(NP,1)*masses; % All particles have the same mass.
M = diag(ms);

% Testing with some initial velocity
V(Nc,:) = [50,50]; % Diagonal velocity of the top right particle.

% Time step set-up.
T = 0.5;
t_steps = T/dt;
ts = 0:dt:T;

% F = @(X,V) ForceFunction(X,V,A,ms,g,ks_springs,kd_springs,L_springs);
F = @(X,V) ForceFunction(X,V,ms,g,ks_springs,kd_springs,L_springs);
[Xs,Vs] = LeapFrog(X,V,F,M,t_steps,dt);
% VisualizeSpringSystem(Xs,A)
% close
% disp("Saved to 'ExampleVideo.avi'")
[E,Ek,Es,Ep] = EnergyCalculation(Xs,Vs,ms,g,ks_springs,L_springs);
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
    %   g - (float) gravitational constant, typically 9.82
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



