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
    % Initialize the position/velocity tensor.
    X = zeros(t_steps,NP,n_dims);
    V = zeros(t_steps,NP,n_dims);
    % Get the position and radii of the circles that make up the surface.
    centers = circle_surface(:,1:end-1);
    radii = circle_surface(:,end); % Shape: (N_circles x 1)

    X(1,:,:) = X_init; % Set initial position.
    F_mat = F(X_init,V_init); % Initial force.
    v = V_init-F_mat*dt/2; % Initialize with half Euler step.
    V(1,:,:) = v;
    M_inv = inv(M); % Compute inverse of diagonal mass matrix, i.e. 1./M.
    for n = 1:t_steps-1
        xs = squeeze(X(n,:,:)); % Remove singleton dimension.
        F_mat = F(xs,v); % Calculate the force matrix.
        v = v+dt*M_inv*F_mat; % Calculate the next v(n+1/2).
        x_new=xs+dt*v;
        % Check if collision before proceeding.
        r = -(centers - permute(x_new, [3 2 1])); % Shape (N_
        pos_diff = vecnorm(r,2,2);
        r_bars = r./pos_diff;
        % pos_diff has Shape: (NP x N_circles)
        % Check which particles are intersecting
        inter = squeeze(pos_diff)<radii;
        [inter_circles,inter_particles] = find(inter==1);
        % Only intersect with one circle. Using small dt should eliminate this
        % issue but for robustness select first cirle from left to right.
        [inter_particles,id] = unique(inter_particles,'first');
        inter_circles = inter_circles(id);
        % Unfortunatly I could only resort to a loop for this part...
        n_hat = zeros(length(id),n_dims);
        bounce = zeros(length(id),1);
        for i = 1:length(id) % Could not find pairwise indexing for multidimensional array
            n_hat(i,:) = r_bars(inter_circles(i),:,inter_particles(i));
            bounce(i) = 2*(radii(inter_circles(i))-pos_diff(inter_circles(i),:,inter_particles(i)));
        end
        v_par = dot(v(inter_particles,:),n_hat,2).*n_hat;
        v(inter_particles,:) = v(inter_particles,:)-2*v_par;
        x_new(inter_particles,:) = x_new(inter_particles,:)+bounce.*n_hat;
        X(n+1,:,:) = x_new;
        V(n+1,:,:) = v;
    end

end

% function [X,V] = CheckSurface(x_new,x_old,v_new,v_old,circle_surface)
%     % This functions is seperate to reduce cluttering.
%     %
%     % INPUT
%     %
%     %   x_new - (mat) The proposed new positional state of the grid.
%     %
%     %   x_old - (mat) The previous position, used to check the direction
%     %                 which the particle might have bounced.
%     %   v_new - (mat) proposed velocity of the next state.
%     %
%     %   v_old - (mat) The previous velocity, used to check the where it
%     %                 bounced.
%     %   circle_surface - (struct) All circles which build up the surface.
%     %
%     
%     [N_circles,n_dims] = size(circle_surface);
%     [NP,n_dims] = size(x_new);
%     n_dims = n_dims - 1; % Last entry of circle_surface is the radius.
% 
%     % We need to check if any of the particles of the grid described with
%     % x_new are inside any surface circle. If so, propel that particle above the
%     % surface as if it had bounced and change the velocity accordingly.
%     
%     % NOTE: There is a chance that a particle is in the intersection of two
%     %       circles.
%     
%     centers = circle_surface(:,1:end-1); % (N_particles x n_dims)
%     radii = circle_surface(:,end);
%     % Calculate relative position between each particle and each circle.
%     pos_diff = x_new-permute(centers,[3 2 1]);
%     % Generate a distance matrix between each particle and each circle.
%     
% 
% 
% end


