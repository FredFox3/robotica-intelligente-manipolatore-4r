function fig = animate_4R(t, theta, L, desiredXY)
% theta: N x 4; desiredXY: N x 2, campionata agli stessi tempi t.
validateattributes(t, {'numeric'}, {'real','finite','vector','nonempty'});
validateattributes(theta, {'numeric'}, {'real','finite','2d'});
validateattributes(L, {'numeric'}, {'real','finite','positive','vector','numel',4});
t = double(t(:));
theta = double(theta);
L = double(L(:));
N = numel(t);
if size(theta,1) ~= N || size(theta,2) ~= 4 || any(diff(t) <= 0)
    error('animate_4R:InvalidSamples', 'Servono N pose e tempi strettamente crescenti.');
end
if nargin < 4
    desiredXY = [];
end
if ~isempty(desiredXY)
    validateattributes(desiredXY, {'numeric'}, {'real','finite','size',[N,2]});
end
alpha = cumsum(theta,2);
actualXY = [cos(alpha)*L, sin(alpha)*L];
fig = figure('Name','Animazione manipolatore 4R','Color','w');
ax = axes('Parent',fig);
hold(ax,'on');
h = plot_robot_4R(theta(1,:),L,ax);
path = animatedline(ax,'Color',[0.1,0.5,0.8],'LineWidth',1.4);
if ~isempty(desiredXY)
    reference = plot(ax,desiredXY(:,1),desiredXY(:,2),'k--','LineWidth',1.2);
    target = plot(ax,desiredXY(end,1),desiredXY(end,2),'mp', ...
        'MarkerSize',12,'MarkerFaceColor','m');
    movingTarget = plot(ax,desiredXY(1,1),desiredXY(1,2),'o', ...
        'Color',[0.9,0.5,0],'MarkerFaceColor',[0.9,0.5,0],'MarkerSize',6);
    legend(ax,[h.links,reference,path,target,movingTarget], ...
        {'Robot','Traiettoria desiderata','Traiettoria effettiva', ...
        'Target finale','Riferimento corrente'},'Location','southoutside');
else
    legend(ax,[h.links,path],{'Robot','Traiettoria effettiva'},'Location','southoutside');
end

% Campionamento dei fotogrammi; la traccia conserva tutti i campioni.
indices = unique(round(linspace(1,N,min(N,241))));
previous = 0;
for k = indices
    if ~isgraphics(fig)
        return;
    end
    h = plot_robot_4R(theta(k,:),L,ax,h);
    addpoints(path,actualXY(previous+1:k,1),actualXY(previous+1:k,2));
    if ~isempty(desiredXY)
        set(movingTarget,'XData',desiredXY(k,1),'YData',desiredXY(k,2));
    end
    title(ax,sprintf('t = %.2f s | DH: x rosso, y verde, z uscente',t(k)));
    drawnow;
    pause(0.02); % Riproduzione indicativa, non vincolata al tempo reale.
    previous = k;
end
end
