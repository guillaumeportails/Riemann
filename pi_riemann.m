
function p = pi_riemann(x, zeros_t)
  % Inversion de Mobius de J(x)

  p = 0;

  K = floor(log(x)/log(2));

  for k = 1:K

    mu = mobius_RS(k);

    if mu != 0
      p = p + mu/k * J_riemann(x^(1/k), zeros_t);
    end

  end

end


function m = mobius_RS(n)

  if n == 1
    m = 1;
    return;
  end

  f = factor(n);

  if length(unique(f)) < length(f)
    m = 0;
  else
    m = (-1)^length(f);
  end

end


