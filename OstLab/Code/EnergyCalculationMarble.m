function [E,Ek,Es,Ep] = EnergyCalculationMarble(X,X_marble,V,V_marble,ms,ms_marble,g,ks,L)
%ENERGYCALCULATIONMARBLE Calculates the energies of the system of springs with the additional marbles.
%   Calculate the various energies of a system.
%
% INPUT
%
%   X - (mat) Positions of each particle at each timestep.
%             Shape: (t_steps x NP x n_dims). Dimensions in order (x,y,z)
%   X_marble - (mat) Positions of each marble at each timestep.
%             Shape: (t_steps x NM x n_dims). Dimensions in order (x,y,z)
%
%   V - (mat) Velocities of each particle at each timestep.
%             Shape: (t_steps x NP x n_dims)
%
%   V_marble - (mat) Velocities of each marble at each timestep.
%             Shape: (t_steps x NM x n_dims)
%
%   ms - (vec) The masses of each particle.
%             Shape: (NP x 1)
%   ms_marble - (vec) The masses of each marble.
%             Shape: (NM x 1)
%
%   g - (float) gravitational constant, typically 9.82
%
%   ks - (mat/float) either matrix of shape (NP x NP) or float.
%                    If matrix then ks(i,j) indicates coefficient of
%                    the spring between particle i and particle j. The
%                    matrix must be symmetric to make sense.
%
%   L - (mat/float)  either matrix of shape (NP x NP) or float.
%                    If matrix then ks(i,j) sym indicates length of
%                    the spring between particle i and particle j at
%                    rest.The matrix must be symmetric to make sense.
%

% OUTPUT
%   
%   E - (mat)  Total energy of the system at each time step.
%              Shape: (t_steps x 1)
%
%   Ek - (mat) Total kinetic energy of the system at each time step.
%              Shape: (t_steps x 1)
%
%   Ep - (mat) Total spring energy of the system at each time step.
%              Shape: (t_steps x 1)
%
%   Es - (mat) Total potential energy of the system at each time step.
%              Shape: (t_steps x 1)
%


[t_steps,NP,n_dims] = size(X);
% SPRING
% Create a distance tensor of shape (t_steps x NP x n_dims x NP)
R = X - permute(X, [1 4 3 2]); % Relative positions, how long are the springs.

rs = vecnorm(R,2,3); % The magnitude of each strings.
% Extent first dimension of ks to match rs and L as (t_steps x NP x 1 x NP)
ks_ext(1,:,1,:) = ks;
% Add timestep dimension to L
L = permute(repmat(L,1,1,1,t_steps),[4,1,2,3]);
spring_energies = squeeze(0.5*ks_ext.*((rs-L).^2)); % The energy of the strings.
% Set all diagonals to zero.
for n =1:n_dims
    spring_energies(:,n,n) = 0;
end


% Divide by two when summing all the springs since it is a symmetric
% matrix.
Es = sum(spring_energies,[2,3])/2;


% KINETIC
vs = vecnorm(V,2,3); % Magnitude of velocity of each particle.
% vs now have the shape (t_steps x NP)
Ek_springs = vs.^2*ms/2; % (t_steps x NP)x(NP x 1) -> (t_steps x 1)
vs_marble = vecnorm(V_marble,2,3);
Ek_marbles = vs_marble.^2*ms_marble/2;
Ek = Ek_springs+Ek_marbles;

% POTENTIAL
h = X(:,:,2*(2==n_dims)+3*(3==n_dims)); % (t_steps x NP)

Ep_springs = g*h*ms; % (t_steps x NP)x(NP x 1) -> (t_steps x 1)
h_marble = X_marble(:,:,2*(2==n_dims)+3*(3==n_dims));
Ep_marble = g*h_marble*ms_marble;
Ep = Ep_springs+Ep_marble;

% TOTAL

E = Ek+Ep+Es;
end

