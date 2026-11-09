function theta = ik_4R(x, y, phi, L, elbow)
% Soluzione analitica; elbow = +1 (default) o -1 sceglie il segno di q_2.
validateattributes(x, {'numeric'}, {'real','finite','scalar'});
validateattributes(y, {'numeric'}, {'real','finite','scalar'});
validateattributes(phi, {'numeric'}, {'real','finite','scalar'});
validateattributes(L, {'numeric'}, {'real','finite','positive','vector','numel',4});
if nargin < 5
    elbow = 1;
end
validateattributes(elbow, {'numeric'}, {'real','finite','scalar'});
if abs(elbow) ~= 1
    error('ik_4R:InvalidElbow', 'elbow deve essere +1 oppure -1.');
end
L = double(L(:));
phi = atan2(sin(double(phi)), cos(double(phi)));

% Il polso e' l'estremo del link 3. Normalizzare limita gli errori di scala.
scale = max(L);
l = L/scale;
w = [double(x)/scale; double(y)/scale] - l(4)*[cos(phi); sin(phi)];
r = hypot(w(1), w(2));
tol = 128*eps(max(1, sum(l)));
if ~isfinite(r)
    error('ik_4R:Unreachable', 'Posa fuori dallo spazio raggiungibile.');
end

% Con psi = q_1+q_2+q_3, u = w-l_3*[cos(psi);sin(psi)] deve
% appartenere all'anello 2R: a <= norm(u) <= b.
a = abs(l(1)-l(2));
b = l(1)+l(2);
rhoMin = max(a, abs(r-l(3)));
rhoMax = min(b, r+l(3));
if rhoMin > rhoMax + tol
    error('ik_4R:Unreachable', 'Posa fuori dallo spazio raggiungibile.');
end
if r == 0
    psi = 0; % Qualunque orientamento e' equivalente se il polso e' all'origine.
else
    % rho^2 = r^2+l_3^2-2*r*l_3*cos(delta), delta = psi-atan2(w_y,w_x).
    % L'intervallo ammissibile di cos(delta) evita di fissare un giunto.
    denominator = 2*r*l(3);
    cMin = (r^2 + (l(3)-b)*(l(3)+b))/denominator;
    cMax = (r^2 + (l(3)-a)*(l(3)+a))/denominator;
    cMin = max(-1, min(1, cMin));
    cMax = max(-1, min(1, cMax));
    cDelta = (cMin+cMax)/2;
    delta = atan2(sqrt(max(0, 1-cDelta^2)), cDelta);
    psi = atan2(w(2), w(1)) + delta;
end
u = w-l(3)*[cos(psi); sin(psi)];
rho = hypot(u(1), u(2));
if rho < a-tol || rho > b+tol
    error('ik_4R:Unreachable', 'Il polso residuo non e'' raggiungibile.');
end

if rho <= tol && abs(l(1)-l(2)) <= tol
    % 2R completamente ripiegato: q_1 e' libero.
    q1 = 0;
    q2 = double(elbow)*pi;
else
    rho = max(a,min(b,rho)); % Correzione dei soli sconfinamenti di arrotondamento.
    % Formula di Erone fattorizzata: piu' stabile di sqrt(1-cos(q_2)^2)
    % quando il 2R e' quasi completamente ripiegato o disteso.
    area4 = sqrt(max(0,(b-rho)*(b+rho)*(rho-a)*(rho+a)));
    signedArea4 = double(elbow)*area4;
    q2 = atan2(signedArea4, rho^2-l(1)^2-l(2)^2);
    q1 = atan2(u(2),u(1)) - atan2(signedArea4, rho^2+(l(1)-l(2))*b);
end
q3 = psi-q1-q2;
q4 = phi-psi;
theta = [q1; q2; q3; q4];
theta = atan2(sin(theta), cos(theta));

% Controllo indipendente tramite il prodotto delle trasformazioni DH.
[~, xf, yf, phif] = fk_4R(theta, L);
positionError = hypot(xf/scale-double(x)/scale, yf/scale-double(y)/scale);
angleError = abs(atan2(sin(phif-phi), cos(phif-phi)));
if positionError > 32*tol || angleError > 1e-10
    error('ik_4R:VerificationFailed', 'La soluzione IK non riproduce la posa.');
end
end
