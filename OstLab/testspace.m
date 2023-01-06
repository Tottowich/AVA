x = squeeze(Xs(400,:,1));
y = squeeze(Xs(400,:,2));
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