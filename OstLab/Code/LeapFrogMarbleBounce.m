function [X,X_marble,V,V_marble] = LeapFrogMarbleBounce(X_init,V_init,X_marble_init,V_marble_init,fixed,springs,M,M_marble,g,t_steps,dt);
%LeapFrogMarbleBounce Calculate the trajectory of the spring system using 'Leap frog'
% method
%
% Author: Theodor Jonsson
% Date: 10/01/2023
% 
%

% INPUT
%
%   X_init - (mat) Initial position of particles,
%           matrix of shape: (number of particles, number of dim)  or
%                             (NP x n_dims)
%   
%   X_marble_init - (mat) Initial position the marbles with radius of each
%             particle matrix of shape: (NM x n_dims+1)
%                            
%
%   V_init - (mat) Inital velocity of particles,
%            same size as X_init.
%
%   V_marble_init - (mat) Inital velocity of marbles,
%            same size as X_init.
%
%   t_steps - (int) Number of time steps to simulate
%
%   F - (function) Anonymous function using position matrix and velocity
%                  matrix. Return the force matrix by f=F(X,V)
%   
%   M - (mat) Diagonal matrix containing masses of corresponding particle.
%
%   M_marble - (mat) Diagonal matrix containing masses of each marble.
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
    NM = size(X_marble_init,1); % Get the number of marbles.
    n_dims = size(X_init,2); % Get the number of dimensions.
    kd = springs.kd;
    ks = springs.ks;
    L = springs.L;
    X_fixed = X_init(fixed,:);
    % Initialize the position/velocity tensor.
    X = zeros(t_steps,NP,n_dims);
    X_marble = zeros(t_steps,NM,n_dims+1);
    V = zeros(t_steps,NP,n_dims);
    V(1,:,:) = V_init;
    V_marble = zeros(t_steps,NM,n_dims);
    V_marble(1,:,:) = V_marble_init;
    % Get the position and radii of the marbles
    X_marble(1,:,:) = X_marble_init;
    X_marble(:,:,end) = repmat(X_marble_init(:,end)',t_steps,1);
    radii = X_marble_init(:,end); % Shape: (NM x 1)
    M_mat = diag(M);
    M_marble_mat = diag(M_marble);
    M_inv = inv(M_mat); % Compute inverse of diagonal mass matrix, i.e. 1./M.
    M_marble_inv = inv(M_marble_mat);
    
    F = @(x_net,v_net,x_marble,v_marble) ForceFunction(x_net,v_net,x_marble,v_marble,M,M_marble,g,ks,kd,L);
    X(1,:,:) = X_init; % Set initial position of the net.
    [F_net,F_marble] = F(X_init,V_init,X_marble_init,V_marble_init); % Initial force.
    v_net = V_init-M_inv*F_net*dt/2; % Initialize with half Euler step.
    v_marble = V_marble_init - M_marble_inv*F_marble*dt/2;
    v_net(fixed,:) = 0;
    min_dist_marbles = radii+radii';
    min_dist_marbles(eye(NM)==1) = NaN;
    for n = 1:t_steps-1
        xs = squeeze(X(n,:,:)); % Remove singleton dimension.
        xs_marble = squeeze(X_marble(n,:,1:end-1));
        if NM==1
            xs_marble = xs_marble';
        end
        [F_net,F_marble] = F(xs,v_net,xs_marble,v_marble); % Calculate the force matricies.
        % These matricies are standalone from each other.
        %%% MARBLE COLLISION
        
        r = -(xs_marble - permute(xs_marble, [3 2 1]));
        pos_diff = vecnorm(r,2,2);
        r_bars = r./pos_diff; % Shape - (NM x n_dims x NM), will be anti symmetric.
        % Replace NaN with zeros. 
        r_bars(isnan(r_bars))=0;
        inter = squeeze(pos_diff)<=min_dist_marbles;
        [inter_marbles_A,inter_marbles_B] = find(inter==1);
        [inter_marbles_B,id] = unique(inter_marbles_B,'first');
        if ~isempty(id)
            inter_marbles_A = inter_marbles_A(id);
            m_sum = sum(M_marble([inter_marbles_A;inter_marbles_B]))/2;
            v_marble(inter_marbles_A,:) = (M_marble(inter_marbles_A)-M_marble(inter_marbles_B))./m_sum.*v_marble(inter_marbles_A,:)+...
                                         2*M_marble(inter_marbles_B).*v_marble(inter_marbles_B,:)./m_sum;
            v_marble(inter_marbles_B,:) = 2*M_marble(inter_marbles_A)./m_sum.*v_marble(inter_marbles_A,:)+...
                                         (M_marble(inter_marbles_B)-M_marble(inter_marbles_A))./m_sum.*v_marble(inter_marbles_B,:);
            xs_marble(inter_marbles_A,:) = xs_marble(inter_marbles_A,:)+v_marble(inter_marbles_A,:)*dt;
            n_hat = zeros(length(id),n_dims);
            for i = 1:length(id) % Could not find pairwise indexing for multidimensional array
                n_hat(i,:) = r_bars(inter_marbles_A(i),:,inter_marbles_B(i));
            end
            pos_diff = squeeze(pos_diff);
            inds = sub2ind(size(pos_diff),inter_marbles_A,inter_marbles_B);
            xs_marble(inter_marbles_A,:) = xs_marble(inter_marbles_A,:)-(radii(inter_marbles_A)-pos_diff(inds)).*n_hat;
            %keyboard
        end
        %%% END MARBLE COLLISION
        
        
        %%% INTERSECTION FLOOR
        r = -(xs_marble - permute(xs, [3 2 1]));
        pos_diff = vecnorm(r,2,2);
        r_bars = r./pos_diff; % Shape (NM x n_dims x NP)
        inter = squeeze(pos_diff)<=radii;
        if NM==1 % Squeeze removes wrong dimension.
            inter = inter';
        end
        [inter_marbles,inter_nodes] = find(inter==1);
        [inter_nodes,id] = unique(inter_nodes,'first');
        inter_marbles = inter_marbles(id);
        n_hat = zeros(length(id),n_dims);
        v_net_new = v_net+dt*M_inv*F_net; % Calculate the next v(n+1/2).
        v_net_new(inter_nodes,:) = 0;
        v_net_new(fixed,:)=0;
        v_marble_new = v_marble + dt*M_marble_inv*F_marble;
        v_marble_new(inter_marbles,:) = 0;
        for i = 1:length(id) % Could not find pairwise indexing for multidimensional array
            n_hat(i,:) = r_bars(inter_marbles(i),:,inter_nodes(i));
            m_marble = M_marble(inter_marbles(i));
            m_node = M(inter_nodes(i));
            v_marble_new(inter_marbles(i),:) = (m_marble-m_node)/(m_marble+m_node)*v_marble(inter_marbles(i),:)+...
                                               (2*m_node)/(m_marble+m_node)*v_net(inter_nodes(i),:);
            v_net_new(inter_nodes(i),:) = v_net_new(inter_nodes(i),:)+2*m_marble/(m_marble+m_node)*v_marble(inter_marbles(i),:)+...
                                          (m_node-m_marble)/(m_marble+m_node)*v_net(inter_nodes(i),:);
        end
        x_marble_new = xs_marble+dt*v_marble_new;
        if ~isempty(inter_nodes)% && length(unique(inter_marbles))>1
            pos_diff = squeeze(pos_diff);
            if NM==1
                pos_diff = pos_diff';
            end
            inds = sub2ind(size(pos_diff),inter_marbles,inter_nodes);
            if NM~=1
                xs(inter_nodes,:) = xs(inter_nodes,:)+2*(radii(inter_marbles)-pos_diff(inds)).*n_hat;
            else
                xs(inter_nodes,:) = xs(inter_nodes,:)+2*(radii(inter_marbles)-pos_diff(inds))'.*n_hat;
            end
        end
        %%% END INTERSECTION FLOOR
        x_new=xs+dt*v_net_new;
        
        %if ~isempty(inter_marbles) && length(inter_marbles)>1
            %Each node can only collide with one marble at a time for
            %simplicity.
        %    disp("Marble: "+inter_marbles)
        %    disp("Node: "+inter_nodes)
        %end
%         x_marble_new(inter_marbles) = 
        X(n+1,:,:) = x_new;
        X(n+1,fixed,:) = X_fixed; % The edges should be still
        V(n+1,:,:) = v_net_new;
        V(n+1,fixed,:) = 0;
        X_marble(n+1,:,1:end-1) = x_marble_new;
        V_marble(n+1,:,:) = v_marble_new;
        v_marble = v_marble_new;
        v_net = v_net_new;
%         n_hats{n+1} = n_hat;
    end
end

function [F_net,F_marble] = ForceFunction(X,V,X_marble,V_marble,ms,ms_marble,g,ks,kd,L)
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
    F_damping = kd.*dot(V_rels,r_bars,2);
    % F_damping(isnan(F_damping))=0;
    % Multiply with the unit vectors of each individual spring.
    F_tensor = -(F_spring+F_damping).*r_bars; % (NP x n_dims x NP)
    % The Entries F(i,:,i) should be zero since this corresponds to the
    % force asserted on particle i on particle i. Which is always zero.
    F_net = sum(F_tensor,3); % Sum along last channel.
                             % Last channel corresponds to each
                             % contribution from each spring
    % F_mat now have the correct shape of (NP x n_dims)
    
    % This has only taken into account the spring system.
    % Now add gravity!
    F_g = ms*g*[0 0 -1]; % (NPx1)x(1x3) => (NPx3)
    F_net = F_net+F_g;
    F_marble = ms_marble*g*[0 0 -1];

    %if any(any(isnan(F_marble)))
    %    disp("NaN marble found!")
    %    keyboard
    %end
    %if any(any(isinf(F_marble)))
    %    disp("Inf marble found!")
    %    keyboard
    %end
    %if any(any(isnan(F_net)))
    %    disp("NaN net found!")
    %    keyboard
    %end
    %if any(any(isinf(F_net)))
    %    disp("Inf net found!")
    %    keyboard
    %end
end

