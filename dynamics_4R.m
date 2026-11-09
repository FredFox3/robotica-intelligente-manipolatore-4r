function [M, C, G] = dynamics_4R(theta, dtheta, params)
% M*ddq + C*dq + G = tau. Baricentri a L_i/2; gravita' lungo -y.
validateattributes(theta, {'numeric'}, {'real','finite','vector','numel',4});
validateattributes(dtheta, {'numeric'}, {'real','finite','vector','numel',4});
if ~isstruct(params) || ~isscalar(params) || ...
        ~all(isfield(params, {'L','m','I','g'}))
    error('dynamics_4R:InvalidParams', 'params deve contenere L, m, I e g.');
end
validateattributes(params.L, {'numeric'}, {'real','finite','positive','vector','numel',4});
validateattributes(params.m, {'numeric'}, {'real','finite','positive','vector','numel',4});
validateattributes(params.I, {'numeric'}, {'real','finite','nonnegative','vector','numel',4});
validateattributes(params.g, {'numeric'}, {'real','finite','nonnegative','scalar'});
q = double(theta(:));
dq = double(dtheta(:));
L = double(params.L(:));
m = double(params.m(:));
I = double(params.I(:));
g = double(params.g);

alpha = cumsum(q);
direction = [cos(alpha).'; sin(alpha).'];
tangent = [-sin(alpha).'; cos(alpha).'];
M = zeros(4);
dM = zeros(4,4,4); % dM(:,:,k) = derivata rispetto a q_k.
G = zeros(4,1);
for i = 1:4
    lengths = L;
    lengths(i) = L(i)/2;
    Jv = zeros(2,4);
    for j = 1:i
        Jv(:,j) = tangent(:,j:i)*lengths(j:i);
    end
    % T_i = (m_i*|Jv_i*dq|^2 + I_i*(sum(dq(1:i)))^2)/2.
    M = M + m(i)*(Jv.'*Jv);
    M(1:i,1:i) = M(1:i,1:i) + I(i);
    % V_i = m_i*g*y_ci, quindi G_i = dV/dq_i (segno positivo).
    G = G + m(i)*g*Jv(2,:).';
    for k = 1:i
        dJv = zeros(2,4);
        for j = 1:i
            first = max(j,k);
            dJv(:,j) = -direction(:,first:i)*lengths(first:i);
        end
        dM(:,:,k) = dM(:,:,k) + m(i)*(dJv.'*Jv + Jv.'*dJv);
    end
end

% Christoffel: c_ijk = (d_k M_ij + d_j M_ik - d_i M_jk)/2.
C = zeros(4);
for i = 1:4
    for j = 1:4
        for k = 1:4
            C(i,j) = C(i,j) + 0.5*(dM(i,j,k) + ...
                dM(i,k,j) - dM(j,k,i))*dq(k);
        end
    end
end
end
