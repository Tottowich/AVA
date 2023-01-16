% This is the first exercise. "Preparatory exercise"
% Author: Theodor Jonsson
% Date: 05/01/2023.
% IMPORTANT - See also 'LeapFrog.m' and 'VisualizeSpringSystem.m'
% 
% Visualization
L = 1;
dt = 0.1;
NP = 2; % Number of nodes in the system.
masses = ones(NP,1); % The masses were given as 1 for all nodes.
M = diag(masses); % Construct the diagonal mass matrix.
v_init = 0.01; % Initial magnitude of veloticties.
X_init = [-L/2 0;L/2 0]; % Initial position
t_steps = 1500; % Number of time steps to simulate
V_init = [v_init 0;0 -v_init]; % Initial Velocity of the nodes.

F = @(x_mat,v_mat) 0.5*(0.5-rand(NP,2))-0.1*v_mat; % Random forces (Specified)
[X,V] = LeapFrog(X_init,V_init,F,M,t_steps,dt);

VisualizeSpringSystem(X,[0 1;1 0],0,"Lab1")
