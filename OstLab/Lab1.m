<<<<<<< HEAD
% This is the first exercise. "Preparatory exercise"
% Author: Theodor Jonsson
% Date: 05/01/2023.
% IMPORTANT - See also 'LeapFrog.m' and 'VisualizeSpringSystem.m'
% 
=======
% Visualization
>>>>>>> 2330c1ce840d5fa5964b55a933699a3f1f30cacb
L = 1;
dt = 0.1;
NP = 3; % Number of nodes in the system.
masses = ones(NP,1); % The masses were given as 1 for all nodes.
M = diag(masses); % Construct the diagonal mass matrix.
invM = inv(M); % Compute inverse of M
v_init = 0.01; % Initial magnitude of veloticties.
<<<<<<< HEAD
t_steps = 10; % Number of time steps to simulate
X_init = [-L/2 0;L/2 0; 0 L/2]; % Initial position
=======
t_steps = 1500; % Number of time steps to simulate
%X = zeros(t_steps,NP,2); % Keeps track of position.
X_init = [-L/2 0;L/2 0; 0 L/2]; % Initial position
%V = zeros(t_steps,NP,2); % Keeps track of velocity.
>>>>>>> 2330c1ce840d5fa5964b55a933699a3f1f30cacb
V_init = [v_init 0;0 -v_init;-v_init/2 -v_init/2]; % Initial Velocity of the nodes.

F = @(x_mat,v_mat) 0.5*(0.5-rand(NP,2))-0.1*v_mat; % Random forces (Specified)
X = LeapFrog(X_init,V_init,F,M,t_steps,dt);

<<<<<<< HEAD
VisualizeSpringSystem(X) % Visualize the sequence
=======
%for t = 1:t_steps
    %scatter(X(t,:,1),X(t,:,2))
    %break
%end
%grid on
VisualizeSpringSystem(X)
>>>>>>> 2330c1ce840d5fa5964b55a933699a3f1f30cacb
