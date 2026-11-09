function [T, x, y, phi] = fk_4R(theta, L)
% DH standard: a_i = L_i, alpha_i = d_i = 0; angoli in radianti.
validateattributes(theta, {'numeric'}, {'real','finite','vector','numel',4});
validateattributes(L, {'numeric'}, {'real','finite','positive','vector','numel',4});
theta = double(theta(:));
L = double(L(:));

T = eye(4);
for i = 1:4
    c = cos(theta(i));
    s = sin(theta(i));
    A = [c, -s, 0, L(i)*c; ...
         s,  c, 0, L(i)*s; ...
         0,  0, 1, 0; ...
         0,  0, 0, 1];
    T = T*A;
end
x = T(1,4);
y = T(2,4);
phi = sum(theta); % Orientamento non avvolto; la rotazione e' periodica.
end
