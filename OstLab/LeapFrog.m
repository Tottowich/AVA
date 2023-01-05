function X = LeapFrog(X_init,V_init,F,M,t_steps,dt)
%LEAPFROG Calculate the trajectory of the spring system using 'Leap frog'
% method
%
<<<<<<< HEAD
% Author: Theodor Jonsson
% Date: 05/01/2023
% 
%
=======
>>>>>>> 2330c1ce840d5fa5964b55a933699a3f1f30cacb
% INPUT
%
%   X_init - (mat) Initial position of particles,
%           matrix of shape - (number of particles, number of dim)
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
%
    NP = size(X_init,1); % Get the number of particles
    n_dims = size(X_init,2); % Get the number of dimensions.
    if n_dims>2
        error('Error\n See LeapFrog3 for higher dimensions, LeapFrog is indended for 2 dimensional not: '+n_dims);
    end
    X = zeros(t_steps,NP,n_dims); % Initialize the position tensor.
    X(1,:,:) = X_init; % Set initial position.
    F_mat = F(X_init,V_init); % Initial force.
    v = V_init-F_mat*dt/2; % Initialize with half Euler step.
    M_inv = inv(M); % Compute inverse of diagonal mass matrix, i.e. 1./M.
    for n = 1:t_steps-1
<<<<<<< HEAD
        F_mat = F(squeeze(X(n,:,:)),v); % Calculate the force matrix.
=======
        F_mat = F(X(n,:,:),v); % Calculate the force matrix.
>>>>>>> 2330c1ce840d5fa5964b55a933699a3f1f30cacb
        v = v+dt*M_inv*F_mat; % Calculate the next v(n+1/2).
        X(n+1,:,:) = squeeze(X(n,:,:))+dt*v;
    end
end

