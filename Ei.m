%% Exponentielle intégrale  = int(-inf,x, exp(t)/t)
% Ei(w) = gamma + log(w) + somme w^k/(k*k!)
% Li(z) = Ei(log(z))
% Ei(0.372507...) = 0
function y = Ei(w)

  y = -real(expint(-w));
  return;

  euler = 0.5772156649015329;

  t = w;
  s = w;

  for k = 2:10000
    t = t*w*(k-1)/(k*k);
    s = s + t;

    if abs(t) < 1e-14*max(1,abs(s))
      break;
    end
  end

  y = euler + log(w) + s;
end

