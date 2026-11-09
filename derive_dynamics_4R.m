function model = derive_dynamics_4R()
% Richiede Symbolic Math Toolbox. Non genera file ne' semplifica globalmente.
q = sym('q', [4,1], 'real');
dq = sym('dq', [4,1], 'real');
ddq = sym('ddq', [4,1], 'real');
L = sym('L', [4,1], 'positive');
m = sym('m', [4,1], 'positive');
I = sym('I', [4,1], 'positive');
g = sym('g', 'positive');

alpha = cumsum(q);
joint = sym(zeros(2,1));
kinetic = sym(0);
potential = sym(0);
for i = 1:4
    direction = [cos(alpha(i)); sin(alpha(i))];
    center = joint + L(i)*direction/2;
    Jv = jacobian(center, q);
    velocity = Jv*dq;
    omega = sum(dq(1:i));
    kinetic = kinetic + (m(i)*(velocity.'*velocity) + I(i)*omega^2)/2;
    potential = potential + m(i)*g*center(2);
    joint = joint + L(i)*direction;
end

momentum = jacobian(kinetic, dq).';
M = jacobian(momentum, dq);
G = jacobian(potential, q).';
tauEL = jacobian(momentum,q)*dq + M*ddq ...
    - jacobian(kinetic,q).' + G;
dM = cell(4,1);
for k = 1:4
    dM{k} = diff(M,q(k));
end
C = sym(zeros(4));
for i = 1:4
    for j = 1:4
        for k = 1:4
            C(i,j) = C(i,j) + (dM{k}(i,j) + ...
                dM{j}(i,k) - dM{i}(j,k))*dq(k)/2;
        end
    end
end

% Valutatore in memoria: [M,C,G,tauEL,T,V] = evaluate(q,dq,ddq,L,m,I,g).
evaluate = matlabFunction(M, C, G, tauEL, kinetic, potential, ...
    'Vars', {q,dq,ddq,L,m,I,g}, 'Optimize', false);
qCases = [0, 0.3, -0.9; 0, -0.7, 0.4; 0, 0.5, -0.2; 0, -0.2, 0.8];
for k = 1:3
    p.L = [0.40;0.32;0.24;0.16] + 0.01*(k-1)*[1;2;1;3];
    p.m = [2;1.5;1;0.6] + 0.1*(k-1)*[1;2;3;1];
    p.I = p.m.*p.L.^2/12;
    p.g = 9.81;
    qn = qCases(:,k);
    dqn = [0.3;-0.5;0.2;0.4]*k;
    ddqn = [-0.2;0.1;0.4;-0.3]/k;
    [Ms,Cs,Gs,ts,Ts,Vs] = evaluate(qn,dqn,ddqn,p.L,p.m,p.I,p.g);
    [Mn,Cn,Gn] = dynamics_4R(qn,dqn,p);
    assert(norm(Ms-Mn,'fro') < 1e-9*max(1,norm(Mn,'fro')), 'M non coerente.');
    assert(norm(Cs-Cn,'fro') < 1e-9*max(1,norm(Cn,'fro')), 'C non coerente.');
    assert(norm(Gs-Gn) < 1e-9*max(1,norm(Gn)), 'G non coerente.');
    expectedTau = Mn*ddqn + Cn*dqn + Gn;
    assert(norm(ts-expectedTau) < 1e-9*max(1,norm(expectedTau)), ...
        'Eulero-Lagrange e Christoffel non coerenti.');
    assert(abs(Ts-dqn.'*Mn*dqn/2) < 1e-9*max(1,abs(Ts)), ...
        'Energia cinetica non coerente.');
    a = cumsum(qn);
    heights = cumsum(p.L.*sin(a)) - p.L.*sin(a)/2;
    Vn = p.g*sum(p.m.*heights);
    assert(abs(Vs-Vn) < 1e-9*max(1,abs(Vn)), 'Energia potenziale non coerente.');
end
model = struct('q',q, 'dq',dq, 'ddq',ddq, 'L',L, 'm',m, 'I',I, ...
    'g',g, 'kinetic',kinetic, 'potential',potential, ...
    'M',M, 'C',C, 'G',G, 'tauEL',tauEL, 'evaluate',evaluate);
fprintf('Derivazione simbolica completata; confronto numerico riuscito in 3 casi.\n');
end
