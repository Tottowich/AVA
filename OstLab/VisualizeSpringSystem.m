function  VisualizeSpringSystem(X)
%VISUALIZECHEESE Visualize the simulated spring system'cheese'
%   Function used to animated the system of springs.
%
<<<<<<< HEAD
% Author: Theodor Jonsson
% Date: 05/01/2023
=======
>>>>>>> 2330c1ce840d5fa5964b55a933699a3f1f30cacb
% INPUT
%
%   X - matrix of size (t_steps, #nodes, dimensions)
%
% Initialize animation file before main loop
<<<<<<< HEAD
figure(1); % Create a new figure, if needed.
=======
figure; % Create a new figure, if needed.
>>>>>>> 2330c1ce840d5fa5964b55a933699a3f1f30cacb
set(gcf,'Position',get(0,'Screensize')); % Maximize the window for quality
MOVE = VideoWriter('ExampleVideo.avi'); % Create an video struct, "MOVE.",
% with the output file ExampleVideo.avi
MOVE.Quality = 100; % Set the quality (0-100)
MOVE.FrameRate = 25; % Set frames per seond to PAL standard (animation speed)
open(MOVE);
set(gca,'nextplot','replacechildren');
<<<<<<< HEAD
xlim([min(X(:,:,1),[],'all'),max(X(:,:,1),[],'all')])
ylim([min(X(:,:,2),[],'all'),max(X(:,:,2),[],'all')])
grid on
title("Visualization of spring connected system")
A = ones(size(X,2)); % Adjecency matrix used to connect all nodes.
...
% Begin main time loop
for t = 1:size(X,1)
    figure(1)
    x = squeeze(X(t,:,1));
    y = squeeze(X(t,:,2));
    xy = [x;y]'; % xy should be (NP x n_dims)
    % Update positions etc. and plot
    gplot(A,xy,'b-o') % Plot connected strings using graph plotting.   
    writeVideo(MOVE,getframe(gcf)); % Get a snapshot of the active figure frame
    % End of main time loop
=======
grid on
...
% Begin main time loop
for t = 1:size(X,1)
    scatter(X(t,:,1),X(t,:,2)) % Update positions etc. and plot
    writeVideo(MOVE,getframe(gcf)); % Get a snapshot of the active figure frame
    ...
    % End of main time loop
    ...
>>>>>>> 2330c1ce840d5fa5964b55a933699a3f1f30cacb
end
close(MOVE); % Close and save the avi-file
end

