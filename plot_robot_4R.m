function h = plot_robot_4R(theta, L, ax, h)
% h = plot_robot_4R(q,L,ax); passare h al quarto argomento per aggiornare.
validateattributes(theta, {'numeric'}, {'real','finite','vector','numel',4});
validateattributes(L, {'numeric'}, {'real','finite','positive','vector','numel',4});
if nargin < 3 || isempty(ax)
    ax = gca;
end
if ~isscalar(ax) || ~isgraphics(ax,'axes')
    error('plot_robot_4R:InvalidAxes', 'ax deve essere un oggetto axes valido.');
end
L = double(L(:));
angles = [0; cumsum(double(theta(:)))];
P = [0, cumsum(L.*cos(angles(2:5))).'; ...
     0, cumsum(L.*sin(angles(2:5))).'];
frameSize = 0.09*sum(L);
ux = frameSize*cos(angles).';
uy = frameSize*sin(angles).';
offset = 0.025*sum(L);

if nargin < 4 || isempty(h)
    held = ishold(ax);
    hold(ax,'on');
    h.links = plot(ax,P(1,:),P(2,:),'-o','Color',[0.12,0.32,0.72], ...
        'LineWidth',2.5,'MarkerFaceColor',[0.12,0.32,0.72],'MarkerSize',7);
    h.x = quiver(ax,P(1,:),P(2,:),ux,uy,0,'Color',[0.8,0.1,0.1], ...
        'LineWidth',1.2,'MaxHeadSize',0.6);
    h.y = quiver(ax,P(1,:),P(2,:),-uy,ux,0,'Color',[0.1,0.6,0.2], ...
        'LineWidth',1.2,'MaxHeadSize',0.6);
    h.z = plot(ax,P(1,:),P(2,:),'ko','MarkerFaceColor','w','MarkerSize',5);
    h.zDot = plot(ax,P(1,:),P(2,:),'k.','MarkerSize',4,'HandleVisibility','off');
    h.labels = gobjects(5,1);
    for j = 1:5
        label = ['{' num2str(j-1) '}, z_' num2str(j-1) ' \odot'];
        h.labels(j) = text(ax,P(1,j)+offset,P(2,j)+offset,label, ...
            'FontSize',8,'Interpreter','tex','HandleVisibility','off');
    end
    axis(ax,'equal');
    radius = 1.18*sum(L);
    xlim(ax,[-radius,radius]);
    ylim(ax,[-radius,radius]);
    grid(ax,'on');
    xlabel(ax,'x [m]');
    ylabel(ax,'y [m]');
    title(ax,'Terne DH: x rosso, y verde, z uscente dal piano');
    legend(ax,[h.links,h.x,h.y,h.z],{'Link e giunti','x_i','y_i','z_i uscente'}, ...
        'Location','southoutside');
    if ~held
        hold(ax,'off');
    end
else
    set(h.links,'XData',P(1,:),'YData',P(2,:));
    set(h.x,'XData',P(1,:),'YData',P(2,:),'UData',ux,'VData',uy);
    set(h.y,'XData',P(1,:),'YData',P(2,:),'UData',-uy,'VData',ux);
    set(h.z,'XData',P(1,:),'YData',P(2,:));
    set(h.zDot,'XData',P(1,:),'YData',P(2,:));
    for j = 1:5
        set(h.labels(j),'Position',[P(:,j).'+offset,0]);
    end
end
end
