function  VisualizeSpringSystem(X,A)
%VISUALIZECHEESE Visualize the simulated spring system'cheese'
%   Function used to animated the system of springs.
%
% Author: Theodor Jonsson
% Date: 05/01/2023
% INPUT
%
%   X - (mat) matrix of size (t_steps, NP, n_dims)
%
%   A - (mat) Adjacency matrix of the spring system.
%
% Initialize animation file before main loop
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
daspect([1,1,1]);
grid on
hold on
% Initialize the grid of springs.
x = squeeze(X(1,:,1));
y = squeeze(X(1,:,2));
springs_linespec = "k.-";
SPRINGS = plot(graph(A),springs_linespec,'XData',x,'YData',y);%,'NodeLabel',{});
% Begin main time loop
figure(1)
for t = 2:10:size(X,1)
    x = squeeze(X(t,:,1));
    y = squeeze(X(t,:,2));
    % Update the grid.
    set(SPRINGS,'XData',x,'YData',y)

    writeVideo(MOVE,getframe(gcf)); % Get a snapshot of the active figure frame
    % End of main time loop
end
close(MOVE); % Close and save the avi-file
hold off
end

