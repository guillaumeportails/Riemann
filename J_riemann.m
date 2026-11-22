%% Formule explicite de Riemann approchant pi(x) :
% J(x) = Li(x) - 2 * somme Li(x^z) - log(2) + I(x)
% z = les zeros de zeta sur la droite critique coté +i
% I(x) = int(x, inf, 1/(t*(t^2-1)*log(t))

function J = J_riemann(x, zeros_t)
  %
  % zeros_t : vecteur des partie imaginaire des zeros critiques

  assert(x > 1);

  L = log(x);

  % Terme principal Li(x)
  J = Ei(L);                    % Li(x)

  % Contribution des zeros non triviaux
  for n = 1:length(zeros_t)
    z = 0.5 + 1i*zeros_t(n);
    z = Ei(z*L);                % Li(x^z)
    J = J - 2*real(z);
  end

  % Terme correctif des zeros triviaux
  J = J - log(2);

  % Contribution integrale
  I = quadgk(@(t) 1./(t.*expm1(2*t)), L, Inf);

  J = J + I;
end



