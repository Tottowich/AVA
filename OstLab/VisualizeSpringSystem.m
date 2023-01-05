function  VisualizeSpringSystem(X)
%VISUALIZECHEESE Visualize the simulated spring system'cheese'
%   Function used to animated the system of springs.
%
% Author: Theodor Jonsson
% Date: 05/01/2023
% INPUT
%
%   X - matrix of size (t_steps, NP, n_dims)
%
% Initialize animation file before main loop
figure(1); % Create a new figure, if needed.
set(gcf,'Position',get(0,'Screensize')); % Maximize the window for quality
MOVE = VideoWriter('ExampleVideo.avi'); % Create an video struct, "MOVE.",
% with the output file ExampleVideo.avi
MOVE.Quality = 100; % Set the quality (0-100)
MOVE.FrameRate = 25; % Set frames per seond to PAL standard (animation speed)
open(MOVE);
set(gca,'nextplot','replacechildren');
% Set the limits to avoid window stuttering.
minminx = min(X(:,:,1),[],'all');
minminy = min(X(:,:,2),[],'all');
maxmaxx = max(X(:,:,1),[],'all');
maxmaxy = max(X(:,:,2),[],'all');
xlim([minminx-1,maxmaxx+1])
ylim([minminy-1,maxmaxy+1])
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
end
grid on
close(MOVE); % Close and save the avi-file
end

