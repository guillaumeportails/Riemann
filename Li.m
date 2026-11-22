%% Logarithme intégral  = int(2,x, t/log(t))
% Li(x) = x/log(x) * (1 + 1!/log(x)^1 + 2!/log(x)^2 + 3!/log(x)^3 ...)
function y = Li(x)
  y = real(-expint(-log(x)));
end

