% x = squeeze(X(1,:,1));
% y = squeeze(X(1,:,2));
% xy = [x;y]'; % xy should be (NP x n_dims)
% % Update positions etc. and plot
% figure(1)
% % gplot(A,xy,'b-o')
% % hold on
% % sc = scatter(x,y,'r')
% % hold off
% pl = plot(graph(A),'r-a','XData',x,'YData',y,'NodeLabel',{})

% [ii,jj] = sparse_adj_matrix([Nr,Nc],1,1,1);
% A2 = sparse(ii, jj, ones(1,numel(ii)), NP, NP);
% A2 = full(A2);
% A2 = A2-diag(diag(A2));
% [r,c] = ndgrid(1:Nr,1:Nc);
% figure(1)
% pl = plot(graph(A2),'r-','XData',r(:),'YData',c(:))
% figure(2)
% A2 = flip(A2,1);
% pl = plot(graph(A2),'r-','XData',r(:),'YData',c(:))
% a([1:10]).center = [1:10,0]
% 
% xv = rand(6,1); yv = rand(6,1);
% xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
% x = rand(1000,1); y = rand(1000,1);
% in = inpolygon(x,y,xv,yv);
% plot(xv,yv,x(in),y(in),'.r',x(~in),y(~in),'.b')
% t = 1;
% x = squeeze(X(t,:,:));
% v = squeeze(V(t,:,:));
% centers = circle_surface(:,1:end-1);
% radii = circle_surface(:,end); % Shape: (N_circles x 1)
% r = -(centers - permute(x, [3 2 1])); % Shape (N_
% pos_diff = vecnorm(r,2,2);
% r_bars = r./pos_diff;
% % pos_diff has Shape: (NP x N_circles)
% % Check which particles are intersecting
% inter = squeeze(pos_diff)<radii;
% [inter_circles,inter_particles] = find(inter==1);
% % Only intersect with one circle. Using small dt should eliminate this
% % issue but for robustness select first cirle from left to right.
% [inter_particles,id] = unique(inter_particles,'first');
% inter_circles = inter_circles(id);
% % Unfortunatly I could only resort to a loop for this part...
% n_hat = zeros(length(id),n_dims);
% bounce = zeros(length(id),1);
% for i = 1:length(id) % Could not find pairwise indexing for multidimensional array
%     n_hat(i,:) = r_bars(inter_circles(i),:,inter_particles(i));
%     bounce(i) = 2*(radii(inter_circles(i))-pos_diff(inter_circles(i),:,inter_particles(i)));
% end
% v_par = dot(v(inter_particles,:),n_hat,2);
% v_new = v(inter_particles,:)-2*v_par;
% x_new = x(inter_particles,:)+bounce.*n_hat;
% figure(1)
% PlotFrame(t,X,A,circle_surface)
% quiver(x(inter_particles,1),x(inter_particles,2),n_hat(:,1),n_hat(:,2),0.2)
% hold off
% % Update
% x(inter_particles,:) = x_new;
% % X(t,:,:) = x;
% PlotFrame(t,X,A,circle_surface)
[X,Y,Z] = sphere(100,2);
surf(X,Y,Z,'EdgeColor','k','FaceColor','w','FaceAlpha',0.5)
A = randi([0,1],[10,10]);
A = A+A';
xs = rand([10,3]);
hold on
plot(graph(A),'XData',xs(:,1),'YData',xs(:,2),'ZData',xs(:,3))
axis equal
hold off
% timeit(@() test(x,v,circle_surface))


function test(x,v,circle_surface)
n_dims = 2;
centers = circle_surface(:,1:end-1);
radii = circle_surface(:,end); % Shape: (N_circles x 1)
r = -(centers - permute(x, [3 2 1])); % Shape (N_
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
bounce = zeros(length(id,1));
for i = 1:length(id) % Could not find pairwise indexing for multidimensional array
    n_hat(i,:) = r_bars(inter_circles(i),:,inter_particles(i));
    bounce(i) = 2*(radii(inter_circles(i))-pos_diff(inter_circles(i),:,inter_particles(i)));
end
v_par = dot(v(inter_particles,:),n_hat,2);
v_new = v(inter_particles,:)-2*v_par;
x_new = x(inter_particles,:)+bounce.*n_hat;
end
