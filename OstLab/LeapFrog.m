function [X,V] = LeapFrog(X_init,V_init,F,M,t_steps,dt)
%LEAPFROG Calculate the trajectory of the spring system using 'Leap frog'
% method
%
% Author: Theodor Jonsson
% Date: 05/01/2023
% 
%

% INPUT
%
%   X_init - (mat) Initial position of particles,
%           matrix of shape - (number of particles, number of dim)  or
%                             (NP x n_dims)
%   
%   V_init - (mat) Inital velocity of particles,
%            same size as X_init.
%
%   t_steps - (int) Number of time steps to simulate
%
%   F - (function) Anonymous function using position matrix and velocity
%                  matrix. Return the force matrix by f=F(X,V)
%   
%   M - (mat) Diagonal matrix containing masses of corresponding particle.
%
%   dt - (float) size of the time steps.
%
% OUTPUT
%
%   X - (mat) Positions of the particles for each timestep.
%             Shape: (t_steps x NP x n_dims)
%   
%   V - (mat) Velocities of the particles for each timestep.
%             Shape: (t_steps x NP x n_dims)
%
    NP = size(X_init,1); % Get the number of particles
    n_dims = size(X_init,2); % Get the number of dimensions.
%     if n_dims>2
%         error('Error\n See LeapFrog3 for higher dimensions, LeapFrog is indended for 2 dimensional not: '+n_dims);
%     end
    X = zeros(t_steps,NP,n_dims); % Initialize the position tensor.
    V = zeros(t_steps,NP,n_dims);
    X(1,:,:) = X_init; % Set initial position.
    F_mat = F(X_init,V_init); % Initial force.
    v = V_init-F_mat*dt/2; % Initialize with half Euler step.
    V(1,:,:) = V_init;
    M_inv = inv(M); % Compute inverse of diagonal mass matrix, i.e. 1./M.
    for n = 1:t_steps-1
        xs = squeeze(X(n,:,:));
        F_mat = F(xs,v); % Calculate the force matrix.
        v = v+dt*M_inv*F_mat; % Calculate the next v(n+1/2).
        V(n+1,:,:)=v;
        X(n+1,:,:) = xs+dt*v;
    end
end

