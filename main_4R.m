clearvars;
close all;
clc;

% SI units. Uniform slender rods, centers of mass at half length, parallel z axes.
% Ideal joints: no friction, elasticity, limits or saturation.
params.L = [0.40;0.32;0.24;0.16];
params.m = [2.0;1.5;1.0;0.6];
params.I = params.m.*params.L.^2/12;
params.g = 9.81;
run_model_checks(params);

qStart = [20;-35;25;10]*pi/180;
qSeed = [70;-55;40;-20]*pi/180;
[~,xGoal,yGoal,phiGoal] = fk_4R(qSeed,params.L);
qGoal = ik_4R(xGoal,yGoal,phiGoal,params.L,-1);
% Choose angle representatives close to the initial configuration.
qGoal = qStart + atan2(sin(qGoal-qStart),cos(qGoal-qStart));
moveTime = 5;
endTime = 8; % Hold the reference after the motion to observe settling.
Kp = diag([64,64,49,36]);
Kd = diag([16,16,14,12]);
initialOffset = [4;-3;2;-2]*pi/180;
state0 = [qStart+initialOffset; zeros(4,1)];
sampleTimes = linspace(0,endTime,801);
options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.02);
rhs = @(t,z) closed_loop(t,z,params,Kp,Kd,qStart,qGoal,moveTime);
[t,state] = ode45(rhs,sampleTimes,state0,options);
q = state(:,1:4);
dq = state(:,5:8);
[qr,dqr,ddqr] = quintic_reference(t,qStart,qGoal,moveTime);
errorQ = qr.'-q;
tau = zeros(numel(t),4);
desiredXY = zeros(numel(t),2);
actualXY = zeros(numel(t),2);
for k = 1:numel(t)
    tau(k,:) = computed_torque(t(k),q(k,:).',dq(k,:).',params, ...
        Kp,Kd,qStart,qGoal,moveTime).';
    [~,xd,yd] = fk_4R(qr(:,k),params.L);
    [~,xa,ya] = fk_4R(q(k,:),params.L);
    desiredXY(k,:) = [xd,yd];
    actualXY(k,:) = [xa,ya];
end

finalJointError = norm(errorQ(end,:),inf);
finalVelocityError = norm(dqr(:,end)-dq(end,:).',inf);
finalPositionError = norm(actualXY(end,:)-[xGoal,yGoal]);
[~,~,~,finalPhi] = fk_4R(q(end,:),params.L);
finalAngleError = abs(atan2(sin(finalPhi-phiGoal),cos(finalPhi-phiGoal)));
rmsError = sqrt(mean(errorQ.^2,1));
maxTorque = max(abs(tau),[],1);
assert(all(isfinite(state(:))) && all(isfinite(tau(:))), 'Simulazione non finita.');
assert(finalJointError < 1e-5, 'Convergenza articolare insufficiente.');
assert(finalVelocityError < 1e-5, 'Velocita'' finale non convergente.');
assert(finalPositionError < 1e-5 && finalAngleError < 1e-5, ...
    'La posa finale non coincide con il target.');
fprintf('\nSimulazione completata: errore iniziale massimo = %.3e rad\n',norm(initialOffset,inf));
fprintf('Errore finale: q = %.3e rad, dq = %.3e rad/s\n',finalJointError,finalVelocityError);
fprintf('Errore di posa finale: %.3e m, %.3e rad\n',finalPositionError,finalAngleError);
fprintf('RMS articolare [rad]:  %s\n',mat2str(rmsError,4));
fprintf('Picco |tau_i| [N m]:   %s\n',mat2str(maxTorque,4));

figure('Name','Posizioni articolari','Color','w');
for j = 1:4
    subplot(2,2,j);
    plot(t,q(:,j),'b-',t,qr(j,:).','k--','LineWidth',1.2);
    grid on;
    xlabel('t [s]'); ylabel(sprintf('q_%d [rad]',j));
    legend('Effettivo','Desiderato','Location','best');
end
figure('Name','Errori articolari','Color','w');
plot(t,errorQ,'LineWidth',1.2);
grid on; xlabel('t [s]'); ylabel('e = q_d - q [rad]');
legend('e_1','e_2','e_3','e_4','Location','best');
figure('Name','Coppie di controllo','Color','w');
plot(t,tau,'LineWidth',1.2);
grid on; xlabel('t [s]'); ylabel('\tau [N m]');
legend('\tau_1','\tau_2','\tau_3','\tau_4','Location','best');

figure('Name','Traiettoria cartesiana','Color','w');
plot(desiredXY(:,1),desiredXY(:,2),'k--', ...
    actualXY(:,1),actualXY(:,2),'b-','LineWidth',1.3);
hold on;
plot(xGoal,yGoal,'mp','MarkerSize',12,'MarkerFaceColor','m');
plot(desiredXY(1,1),desiredXY(1,2),'go','MarkerFaceColor','g');
axis equal; grid on; xlabel('x [m]'); ylabel('y [m]');
legend('Traiettoria desiderata','Traiettoria effettiva','Target finale', ...
    'Inizio riferimento','Location','best');
title('Target cartesiano della traiettoria articolare quintica');
figure('Name','Configurazione iniziale e terne DH','Color','w');
plot_robot_4R(state0(1:4),params.L,gca);
drawnow;
animate_4R(t,q,params.L,desiredXY);

function [qr,dqr,ddqr] = quintic_reference(t,q0,qf,T)
% Zero velocity and acceleration at both ends; hold after T.
r = min(1,max(0,t(:).'/T));
s = 10*r.^3-15*r.^4+6*r.^5;
ds = (30*r.^2-60*r.^3+30*r.^4)/T;
dds = (60*r-180*r.^2+120*r.^3)/T^2;
delta = qf-q0;
qr = q0+delta*s;
dqr = delta*ds;
ddqr = delta*dds;
end

function [tau,M,C,G] = computed_torque(t,q,dq,p,Kp,Kd,q0,qf,T)
[qr,dqr,ddqr] = quintic_reference(t,q0,qf,T);
[M,C,G] = dynamics_4R(q,dq,p);
command = ddqr+Kd*(dqr-dq)+Kp*(qr-q);
tau = M*command+C*dq+G;
end

function dz = closed_loop(t,z,p,Kp,Kd,q0,qf,T)
dq = z(5:8);
[tau,M,C,G] = computed_torque(t,z(1:4),dq,p,Kp,Kd,q0,qf,T);
ddq = M\(tau-C*dq-G); % Full dynamics, without substituting the ideal error dynamics.
dz = [dq;ddq];
end

function run_model_checks(p)
qCases = [zeros(4,1), [0.3;-0.7;0.5;-0.2], ...
    [-1.2;0.4;-0.8;1.5], [pi/2;0;0;0], ...
    reshape(sin(1:48),4,12)*pi];
for k = 1:size(qCases,2)
    q = qCases(:,k);
    a = cumsum(q);
    [T,x,y,phi] = fk_4R(q,p.L);
    assert(norm([x;y]-[sum(p.L.*cos(a));sum(p.L.*sin(a))]) < 1e-12);
    R = [cos(phi),-sin(phi);sin(phi),cos(phi)];
    assert(norm(T(1:2,1:2)-R,'fro') < 1e-12);
    assert(norm(T(3:4,:)-[0,0,1,0;0,0,0,1],'fro') < 1e-12);
    for elbow = [-1,1]
        qi = ik_4R(x,y,phi,p.L,elbow);
        check_pose(qi,p.L,[x;y;phi]);
        assert(elbow*sin(qi(2)) >= -1e-12);
    end

    dq = [0.4;-0.3;0.2;-0.1];
    [M,C,G] = dynamics_4R(q,dq,p);
    assert(isequal(size(M),[4,4]) && isequal(size(C),[4,4]) && isequal(size(G),[4,1]));
    assert(norm(M-M.','fro') < 1e-12*max(1,norm(M,'fro')));
    [~,flag] = chol(M);
    assert(flag == 0, 'M non definita positiva.');
    [kinetic,~] = energy_4R(q,dq,p);
    assert(abs(kinetic-dq.'*M*dq/2) < 1e-12*max(1,abs(kinetic)));
    h = 1e-6;
    gradV = zeros(4,1);
    for j = 1:4
        step = zeros(4,1); step(j) = h;
        [~,Vp] = energy_4R(q+step,dq,p);
        [~,Vm] = energy_4R(q-step,dq,p);
        gradV(j) = (Vp-Vm)/(2*h);
    end
    assert(norm(G-gradV,inf) < 1e-6*max(1,norm(G,inf)));
    Mp = dynamics_4R(q+h*dq,dq,p);
    Mm = dynamics_4R(q-h*dq,dq,p);
    skew = (Mp-Mm)/(2*h)-2*C;
    assert(norm(skew+skew.','fro') < 1e-6*max(1,norm(M,'fro')));
    [~,Czero,~] = dynamics_4R(q,zeros(4,1),p);
    assert(norm(Czero,'fro') == 0);
end

% Wrist at the origin, folded 2R, outer and inner 3R boundaries.
specialL = {[1;1;1;0.2], [1;1;2;0.2], [3;0.5;0.5;0.2], [1;1;1e-6;0.2]};
specialPoses = {[0.2;0;0], [0.2;0;0], [2.2;0;0], [-0.200001;0;pi]};
for k = 1:numel(specialL)
    for elbow = [-1,1]
        target = specialPoses{k};
        qi = ik_4R(target(1),target(2),target(3),specialL{k},elbow);
        check_pose(qi,specialL{k},target);
    end
end
% Pose obtained with the first 2R fully folded.
Lf = [1;1;1;0.2];
[~,xf,yf,phif] = fk_4R([0;pi;0;0],Lf);
for elbow = [-1,1]
    qi = ik_4R(xf,yf,phif,Lf,elbow);
    check_pose(qi,Lf,[xf;yf;phif]);
end
must_be_unreachable(2*sum(p.L),0,0,p.L);
must_be_unreachable(0.2,0,0,[3;0.5;0.5;0.2]);
fprintf('Controlli FK/IK, irraggiungibilita'', energie, M, C e G completati.\n');
end

function check_pose(q,L,target)
[~,x,y,phi] = fk_4R(q,L);
assert(norm([x;y]-target(1:2)) < 1e-10*max(1,sum(L)));
assert(abs(atan2(sin(phi-target(3)),cos(phi-target(3)))) < 1e-10);
end

function must_be_unreachable(x,y,phi,L)
recognized = false;
try
    ik_4R(x,y,phi,L);
catch exception
    if strcmp(exception.identifier,'ik_4R:Unreachable')
        recognized = true;
    else
        rethrow(exception);
    end
end
assert(recognized, 'Posa irraggiungibile non riconosciuta.');
end

function [kinetic,potential] = energy_4R(q,dq,p)
% Geometric recursion independent of the assembly of M and G.
a = cumsum(q);
omega = cumsum(dq);
joint = zeros(2,1);
velocity = zeros(2,1);
kinetic = 0;
potential = 0;
for i = 1:4
    u = [cos(a(i));sin(a(i))];
    tangent = [-sin(a(i));cos(a(i))];
    center = joint+p.L(i)*u/2;
    vc = velocity+p.L(i)*omega(i)*tangent/2;
    kinetic = kinetic+(p.m(i)*(vc.'*vc)+p.I(i)*omega(i)^2)/2;
    potential = potential+p.m(i)*p.g*center(2);
    joint = joint+p.L(i)*u;
    velocity = velocity+p.L(i)*omega(i)*tangent;
end
end
