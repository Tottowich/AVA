function  VisualizeSpringSystemWithSurface(X,A,circle_surface,record,name)
%VISUALIZESPRINSYSTEMWITHSURFACE Visualize the simulated spring system.
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
%   circle_surface - (mat) Matrix of shape (N x n_dims+1). For each circle
%                          there is center (x,y,z) + radius.
%   record - (bool) 0 if not record.
%
%   name - (string) Name of recording
%

% Initialize animation file before main loop
[t_steps,NP,n_dims] = size(X);
set(gcf,'Position',get(0,'Screensize')); % Maximize the window for quality
if record
    MOVE = VideoWriter(name+".avi"); % Create an video struct, "MOVE.",
    % with the output file ExampleVideo.avi
    MOVE.Quality = 60; % Set the quality (0-100)
    MOVE.FrameRate = 25; % Set frames per seond to PAL standard (animation speed)
    open(MOVE);
end
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
% Plot the system as a graph using the adjacency matrix generated for the grid.
% This is done to skip creating a line for each spring connecting the grid.
SPRINGS = plot(graph(A),springs_linespec,'XData',x,'YData',y,'NodeLabel',{});

% Display the circular floor.
N_cirlces = length(circle_surface);
for n = 1:N_cirlces
    % Create a cirlce as a rectangle.
    % Since the rectangle has position based bottom left, and top right.
    % circle_surface contains the center of the circle, find ["bottom
    % left","top right"] for the circle.
    c = circle_surface(n,1:n_dims);
    r = circle_surface(n,end);
    pos = [c(1)-r,c(2)-r,2*r,2*r];
    circles(n) = rectangle(Curvature=[1,1],Position=pos);
    % Create the rectangle which will be a circle.

end
% Begin main time loop
figure(1)
for t = 2:size(X,1)
    x = squeeze(X(t,:,1));
    y = squeeze(X(t,:,2));
    % Update the grid.
    set(SPRINGS,'XData',x,'YData',y)
    if record
        writeVideo(MOVE,getframe(gcf)); % Get a snapshot of the active figure frame
    else
        pause(0.01)
    end
        % End of main time loop
end
hold off
if record
    close(MOVE); % Close and save the avi-file
end

end

