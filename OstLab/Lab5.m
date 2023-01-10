% Modeling a sliding 3D object
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
Nx = 8; % Number of particles in x direction
Ny = 4; %
Nz = 3; % 
masses = 1; % All particles have mass 1.
ks = 500;
kd = 25;
g = 10;
dt = 2e-3;
L = 1; % Evenly distributed particles => sqrt(2) on diagonal.
n_dims = 3;
Nx_circles = 10; % Number of circles that make up the floor in the x-direction
Ny_circles = 6; % Number of circles that make up the floor in the y-direction
dist_circle = 0.1; % [0<->1] describing how large portion of radius to seperate.
r_circle = L; % Randius of the circles which constructs the surface
v_init = [2,0,-20];
% --------------------------------------
start_x = 0;
start_y = 0;
start_z = r_circle*2;
NP = Nx*Ny*Nz; % Total number of particles in the spring grid.
% Time step set-up.
T = 1.5;
t_steps = T/dt;
ts = 0:dt:T-dt;

% ------- Set up the 2D object --------
% The object is denoted by the matrix X
x = 0:L:(Nx-1)/L;
y = 0:L:(Ny-1)/L;
z = (Nz-1)/L:-L:0;
[xs,ys,zs] = meshgrid(x,y,z);
xs = xs+start_x;
ys = ys+start_y;
zs = zs+start_z;
% X = meshgrid(x,y,z);
% scatter3(X)
X_init = cat(4,xs,ys,zs);
X_init = reshape(X_init,[NP n_dims]); % Flatten the matrix.

% X now has Shape (NP x n_dims)

% sz = size(xs);
% [ii,jj] = sparse_adj_matrix(sz, L, inf, 1);
% A = sparse(ii, jj, ones(1,numel(ii)), prod(sz), prod(sz));
% A = A-diag(diag(A));
% image(A*256)
% sz = size(x);
% [ii,jj] = sparse_adj_matrix(sz, L, inf, 1);
% A2 = sparse(ii, jj, ones(1,numel(ii)), prod(sz), prod(sz));
% A2 = A2-diag(diag(A2));
% subplot(1,2,1)
% image(A*256);
% subplot(1,2,2)
% image(A2*256)
% gplot(A,[X_init(:,1),X_init(:,2)])
% figure(1)
% plot(graph(full(A)),'k.-','XData',X_init(:,2),'YData',X_init(:,1),'ZData',X_init(:,3));%,'NodeLabel',{});
% axis padded
% xlabel('x')
% ylabel('y')
% zlabel('z')
% figure(2)
% image(full(A)*256)
% plot3(graph(full(A)),X_init(:,1),X_init(:,2),X_init(:,3));%,'NodeLabel',{});
% plot3(X_init(:,1),X_init(:,2),X_init(:,3));
% figure(3)
% scatter3(X_init(:,1),X_init(:,2),X_init(:,3))
% daspect([1,1,1]);
% hold on;
% scatter3(X_init(1,1),X_init(1,2),X_init(1,3),'r')
% scatter3(X_init(2,1),X_init(2,2),X_init(2,3),'g')
% scatter3(X_init(Ny+1,1),X_init(Ny+1,2),X_init(Ny+1,3),'r');
% scatter3(X_init(Ny*Nx,1),X_init(Ny*Nx,2),X_init(Ny*Nx,3),'g')
% for i = 1:size(graph(A).Edges,1)
%     i1 = graph(A).Edges(i,1).EndNodes(1);
%     i2 = graph(A).Edges(i,1).EndNodes(2);
%     line([X_init(i1,1),X_init(i2,1)],[X_init(i1,2),X_init(i2,2)],[X_init(i1,3),X_init(i2,3)]);
% end
% % point = scatter3(X_init(1,1),X_init(1,2),X_init(1,3),'r');
% % for i = 1:Nx*Ny*Nz
% %     %scatter3(X_init(i,1),X_init(i,2),X_init(i,3),'g')
% %     set(point,'XData',X_init(i,1),'YData',X_init(i,2),'ZData',X_init(i,3))
% %     pause(0.5)
% % end
% hold off;
% axis padded
% xlabel('x')
% ylabel('y')
% zlabel('z')
% As = A;
%
% [A,diagonals] = GridAdjacencyMatrix3D(Nx,Ny,Nz);
% %%
% horizontal = repmat([ones(Nx-1, 1); 0], Ny, 1);% Make the first diagonal vector
%                                                % (for horizontal connections)
% 
% horizontal = horizontal(1:end-1);              % Remove the last value
% 
% tl_to_br = [0; horizontal(1:(Nx*(Ny-1)))]; % Make the second diagonal
%                                    % vector (for top left to bottom right)
% vertical = ones(Nx*(Ny-1), 1);         % Make the third diagonal vector
%                                        %   (for vertical connections)
% 
% bl_to_tr = tl_to_br(2:end-1); % Make the fourth diagonal vector
%                                        %   (for bottom left to top right)
% A = diag(horizontal, 1)+...            % Add the diagonals to a zero matrix
%     2*diag(tl_to_br, Nx-1)+...
%     3*diag(vertical, Nx)+...
%     4*diag(bl_to_tr, Nx+1);
% A = A+A.'; % This makes the Adjacency matrix symmetric.
% diagonals = diag(tl_to_br, Nx-1)+diag(bl_to_tr, Nx+1);
% diagonals = diagonals+diagonals';
% %figure(4)
% c = 256/4;
% %image(A*c)
% % 3D Should contain 13 different diagonals of the matrix.
% % This is because it can be a total of 26 connections to a single particle.
% x_direction = repmat([ones(Ny*(Nx-1),1);zeros(Ny,1)],Nz,1);
% x_direction = x_direction(1:end-Ny);
% y_direction = repmat([ones(Ny-1, 1); 0], Nx*Nz, 1);
% % Remove last
% y_direction = y_direction(1:end-1);
% z_direction = repmat(ones(Ny*Nx, 1), Nz-1, 1);
% % z_direction = z_direction(1:end-1);
% % Diagonals withinplane xy
% bl_to_tr_xy = repmat([ones(Ny-1, 1); 0], Nx-1, 1);
% bl_to_tr_xy = repmat([bl_to_tr_xy;zeros(Ny,1)],Nz,1);
% bl_to_tr_xy = bl_to_tr_xy(1:end-Ny-1);
% 
% tl_to_br_xy = repmat([0;ones(Ny-1, 1)], Nx-1, 1);
% tl_to_br_xy = repmat([tl_to_br_xy;zeros(Ny,1)],Nz,1);
% tl_to_br_xy = tl_to_br_xy(1:end-Ny+1);
% % Diagonals xz
% tl_to_br_xz = repmat([ones(Ny*(Nx-1),1);zeros(Ny,1)],Nz-1,1);
% tl_to_br_xz = tl_to_br_xz(1:end-Ny);
% bl_to_tr_xz = repmat([zeros(Ny,1);ones(Ny*(Nx-1),1)],Nz-1,1);
% bl_to_tr_xz = [bl_to_tr_xz;zeros(Ny,1)];
% % Diagonals yz
% tr_to_bl_yz = repmat([ones(Ny-1,1); 0],Nx*(Nz-1),1);
% tr_to_bl_yz = tr_to_bl_yz(1:end-1);
% 
% tl_to_br_yz = repmat([0;ones(Ny-1,1)],Nx*(Nz-1),1);
% tl_to_br_yz(end+1) = 0;
% 
% % multidimensional diagonals
% % Numbers indicate a cube with 1-4 clockwise on the top
% % 
% one_to_seven = repmat([ones(Ny-1,1);0],Nx-1,1);
% one_to_seven = repmat([one_to_seven;zeros(Ny,1)],Nz-1,1);
% 
% one_to_seven = one_to_seven(1:end-Ny-1);
% 
% two_to_eight = repmat([0;ones(Ny-1,1)],Nx-1,1);
% two_to_eight = repmat([two_to_eight;zeros(Ny,1)],Nz-1,1);
% 
% two_to_eight = two_to_eight(1:end-Ny+1);
% 
% three_to_five = repmat([0;ones(Ny-1,1)],Nx-1,1);
% three_to_five = repmat([zeros(Ny,1);three_to_five],Nz-1,1);
% three_to_five = [three_to_five;zeros(Ny+1,1)];
% 
% four_to_six = repmat([ones(Ny-1,1);0],Nx-1,1);
% four_to_six = repmat([zeros(Ny,1);four_to_six],Nz-1,1);
% four_to_six = [four_to_six;zeros(Ny-1,1)];
% 
% 
% size(diag(x_direction,Ny))
% 
% %%
% figure(5)
% size(diag(x_direction,Ny))
% size(diag(z_direction,Nx*Ny))
% A = diag(x_direction,Ny)+...
%     diag(y_direction,1)+...
%     diag(z_direction,Nx*Ny)+...
%     diag(bl_to_tr_xy,Ny+1)+...
%     diag(tl_to_br_xy,Ny-1)+...
%     diag(tl_to_br_xz,Ny*Nx+Ny)+...
%     diag(bl_to_tr_xz,Nx*Ny-Ny)+...
%     diag(tr_to_bl_yz,Ny*Nx+1)+...
%     diag(tl_to_br_yz,Nx*Ny-1)+...
%     diag(one_to_seven,Nx*Ny+Ny+1)+...
%     diag(two_to_eight,Nx*Ny+Ny-1)+...
%     diag(three_to_five,Nx*Ny-Ny-1)+...
%     diag(four_to_six,Nx*Ny-Ny+1);
% small_diagonals = diag(bl_to_tr_xy,Ny+1)+...
%                 diag(tl_to_br_xy,Ny-1)+...
%                 diag(tl_to_br_xz,Ny*Nx+Ny)+...
%                 diag(bl_to_tr_xz,Nx*Ny-Ny)+...
%                 diag(tr_to_bl_yz,Ny*Nx+1)+...
%                 diag(tl_to_br_yz,Nx*Ny-1);
% 
% large_diagonals = diag(one_to_seven,Nx*Ny+Ny+1)+...
%                 diag(two_to_eight,Nx*Ny+Ny-1)+...
%                 diag(three_to_five,Nx*Ny-Ny-1)+...
%                 diag(four_to_six,Nx*Ny-Ny+1);
% % A = A-diag(bl_to_tr_xy,Ny+1);
% % A = A-diag(tl_to_br_xy,Ny);
% 
% image(small_diagonals*128+large_diagonals*256)
% % A = A+A';
% figure(6)
% plot(digraph(full(A)),'k.-','XData',X_init(:,1),'YData',X_init(:,2),'ZData',X_init(:,3));%,'NodeLabel',{});
% axis padded
% xlabel('x')
% ylabel('y')
% zlabel('z')
% %%
% g = digraph(diag(tl_to_br_xy,Ny-1));
% figure(7)
% plot(g,'k.-','XData',X_init(:,1),'YData',X_init(:,2),'ZData',X_init(:,3));%,'NodeLabel',{});
% axis padded
% xlabel('x')
% ylabel('y')
% zlabel('z')
%

% First Nc entries in X is the top row, Nc+1->2*Nc second row and so on.
V_init = zeros(NP,n_dims);
% Testing with some initial velocity
% V_init(ceil(Nc/2)*3,:) = [5,5]; % Diagonal velocity of the top right particle.
V_init = V_init+v_init;

% -------------------------------------

% Now we can construct the adjecency matrix for the list of particles.
% Which we can use as a mask to remove the forces from non connected
% particles.
[A,small_diagonals,large_diagonals] = GridAdjacencyMatrix3D(Nx,Ny,Nz);
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

% Masses of the particles.
ms = ones(NP,1)*masses; % All particles have the same mass.
M = diag(ms);

% Create the floor
circle_surface = BuildSurface2D(Nx_circles,Ny_circles,r_circle,dist_circle,n_dims);

% F = @(X,V) ForceFunction(X,V,A,ms,g,ks_springs,kd_springs,L_springs);
F = @(X,V) ForceFunction(X,V,ms,g,ks_springs,kd_springs,L_springs);
[X,V] = LeapFrogWithSurfaceCheck(X_init,V_init,F,M,circle_surface,t_steps,dt);
% t3d = timeit(@() LeapFrogWithSurfaceCheck(X_init,V_init,F,M,circle_surface,t_steps,dt))
% Time to simulate using T=8,dt=2*10^-3,
% Nx=8,Ny=4,Nz=3,Nx_circles=10,Ny_circles=3.
% t3d = 2.6591 s
%%
figure(1)
record = 1; % 1 to record
VisualizeSpringSystemWithSurface3D(X,A,circle_surface,record)% close
%%
% disp("Saved to 'ExampleVideo.avi'")
figure(2)
[E,Ek,Es,Ep] = EnergyCalculation(X,V,ms,g,ks_springs,L_springs);
PlotEnergies(E,Ek,Es,Ep,ts,kd)
figure(3);
% Track center of mass.
center_mass_pos = squeeze(sum(X.*ms',2));
center_mass_vel = squeeze(sum(V.*ms',2))./sum(ms);
plot(ts,center_mass_vel(:,1))
grid on
title("Vx center of mass")
figure(4)
plot(ts,center_mass_pos(:,1))
grid on
title("X center of mass")

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



