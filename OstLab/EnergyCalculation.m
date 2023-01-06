function [E,Ek,Es,Ep] = EnergyCalculation(X,V,ms,g,ks,L)
%ENERGYCALCULATION Calculates the energies of the system.
%   Calculate the various energies of a system of springs.
%
% INPUT
%
%   X - (mat) Positions of each particle at each timestep.
%             Shape: (t_steps x NP x n_dims). Dimensions in order (x,y,z)
%
%   V - (mat) Velocities of each particle at each timestep.
%             Shape: (t_steps x NP x n_dims)
%   ms - (vec) The masses of each particle.
%             Shape: (NP x 1)
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
keyboard
spring_energies = squeeze(0.5*ks.*((rs-L).^2)); % The energy of the strings.
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
Ek = vs.^2.*ms/2; % (t_steps x NP)x(NP x 1) -> (t_steps x 1)

% POTENTIAL

ys = X(:,:,2); % (t_steps x NP)
Ep = g*ys.*ms; % (t_steps x NP)x(NP x 1) -> (t_steps x 1) 

% TOTAL

E = Ek+Ep+Es;
end

