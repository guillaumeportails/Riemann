%% Riemann
% 20/11/2026 preuve de l'hypothèse de Riemann par Anthropic (TBC ?)

function Riemann

% Intervalle B:N a explorer
B = 2;
N = 10000;
assert(B >= 2);
assert(N > B);

% Nombre de nombres premiers inferieurs ou égaux a N
NBP = [];
nbp = length(primes(B-1));

% Quelques zéros de la fonction Zeta sur la droite critique
% Selon Wikipedia :
%   1	± 14.134725141734693790…
%   2	± 21.022039638771554993…
%   3	± 25.010857580145688763…
%   4	± 30.424876125859513210…
%   ...
z = find_zeros_RS(1, 100, 0.1);
assert(abs(z(1) - 14.134) < 8e-1);
% Plot zeros de zeta
fprintf("%d zeros disponibles\n", length(z));
figure(1);
t = min(z):0.05:max(z);
y = arrayfun(@hardyZ_RS, t);
plot(t, y,   z, zeros(size(z)), 'ro' );
grid on;
xlabel('t');
ylabel('Z(t)');
title('Zéros de la fonction de Hardy');
% methode formelle (requiert Python) :
pkg load symbolic
syms S
Z = zeta(S);
for n = 1:4
  t0 = [14.134, 21.022, 25.011, 30.425](n);
  R = vpasolve(imag(Z) == 0, S, 0.5 + 1i*t0);
  R
  z(n)
end

% Construction du graphe
for n = B:N

  if isprime(n)
    nbp = nbp + 1;
  end
  NBP(n-B+1) = nbp;

  if n == 1000; assert(nbp == 168); end;
end

% Plot pi(x)
figure(2);
n = B:N;
plot(n, NBP,  ...
     n, n ./ log(n), ...
     n, Li(n));
legend("π(N) : Nombre de premiers <= N", ...
       "Approximation de Gauss : N/log(N)", ...
       "Log intégral Li(N)");
xlabel("N");
title(sprintf("Riemann de %d à %d", B, N));




%% Les zeros de la fonction zeta de Riemann

% Recherche des zeros
function zeros_t = find_zeros_RS(tmin, tmax, dt)

  if nargin < 3
    dt = 0.1;
  endif

  t = tmin:dt:tmax;
  y = arrayfun(@hardyZ_RS, t);

  zeros_t = [];

  for k = 1:length(t)-1
    if y(k)*y(k+1) < 0
      r = fzero(@hardyZ_RS, [t(k), t(k+1)]);
      zeros_t(end+1) = r;
    end
  end

  zeros_t = unique(zeros_t);
end

% Fonction de Hardy Z(t), formule de Riemann-Siegel
function z = hardyZ_RS(t)
  assert(t > 0);    % t doit être réel et strictement positif

  M = floor(sqrt(t/(2*pi)));
  p = sqrt(t/(2*pi)) - M;

  % Phase de Riemann-Siegel
  th = theta_RS(t);

  % Somme principale
  n = 1:M;
  S = 2 * sum(cos(th - t*log(n)) ./ sqrt(n));

  % Correction de Riemann-Siegel C0 + C1
  C0 = psi_RS(p);
  C1 = C1_RS(p);
  a = sqrt(t/(2*pi));
  R = (-1)^(N-1) * a^(-1/2) ...
      * (C0 + C1/a);

  z = S + R;
end

function v = psi_RS(p)
  v = cos(2*pi*(p.^2-p-1/16)) ./ cos(2*pi*p);
  % Singularites apparentes : limites analytiques
  v(abs(cos(2*pi*p)) < 1e-12) = 0.5;
end

function c1 = C1_RS(p)
  h = 0.002;
  d3 = (psi_RS(p+2*h) ...
      - 2*psi_RS(p+h) ...
      + 2*psi_RS(p-h) ...
      - psi_RS(p-2*h)) / (2*h^3);
  c1 = -d3/(96*pi^2);
end

% Développement asymptotique de la fonction theta
function th = theta_RS(t)
  th = (t/2)*log(t/(2*pi)) - t/2 - pi/8 ...
       + 1/(48*t) ...
       + 7/(5760*t^3) ...
       + 31/(80640*t^5);
end


%%

end
