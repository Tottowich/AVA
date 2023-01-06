function [X,V] = LeapFrogWithSurfaceCheck(X_init,V_init,F,M,circle_surface,t_steps,dt)
%LEAPFROGWITHSURFACECHECK Calculate the trajectory of the spring system using 'Leap frog'
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
%   circle_surface - (struct) Structure of circles with 'center' and
%                             'radius', these make up the bottom surface.
%
%   t_steps - (int) Number of time steps to simulate.
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
    if n_dims>2
        error('Error\n See LeapFrog3 for higher dimensions, LeapFrog is indended for 2 dimensional not: '+n_dims);
    end
    X = zeros(t_steps,NP,n_dims); % Initialize the position tensor.
    V = zeros(t_steps,NP,n_dims);
    X(1,:,:) = X_init; % Set initial position.
    F_mat = F(X_init,V_init); % Initial force.
    v = V_init-F_mat*dt/2; % Initialize with half Euler step.
    V(1,:,:) = v;
    M_inv = inv(M); % Compute inverse of diagonal mass matrix, i.e. 1./M.
    for n = 1:t_steps-1
        xs = squeeze(X(n,:,:));
        F_mat = F(xs,v); % Calculate the force matrix.
        v = v+dt*M_inv*F_mat; % Calculate the next v(n+1/2).
        v_new=v;
        x_new=xs+dt*v;
        % Check if collision before proceeding.
        [xs,v] = CheckSurface(x_new,xs,v_new,v_old,circle_surface);
        X(n+1,:,:) = xs+dt*v;
    end
end

function [X,V] = CheckSurface(x_new,x_old,v_new,v_old,circle_surface)
    % This functions is seperate to reduce cluttering.
    %
    % INPUT
    %
    %   x_new - (mat) The proposed new positional state of the grid.
    %
    %   x_old - (mat) The previous position, used to check the direction
    %                 which the particle might have bounced.
    %   v_new - (mat) proposed velocity of the next state.
    %
    %   v_old - (mat) The previous velocity, used to check the where it
    %                 bounced.
    %   circle_surface - (struct) All circles which build up the surface.
    %

    [N_circles,n_dims] = size(circle_surface);
    % We need to check if any of the particles of the grid described with
    % x_new are inside any surface circle. If so, propel that particle above the
    % surface as if it had bounced and change the velocity accordingly.
    
    % NOTE: One particle can only be inside one particle at a time.
    
    % Generate a distance matrix between each particle and each circle.
    


end


