%% The Bourrin Way

function Riemann2

%pkg install -forge video;
pkg load video;

% 1ers zeros de zeta sur la droite critique  ( abs(imag() )
% cf https://www.lmfdb.org/zeros/zeta/
%  14.1347251417
%  21.0220396388
%  25.0108575801
%  30.4248761259
%  ...
g = dlmread("zerozeta.txt");
g = g(1:50, 2);

% Sanity Checks
assert(abs(real(Ei(0.372507))) < 1e-5);       % L'unique zero de Ei
x = 2 + 3*rand + 1i*3*rand;
assert(abs(-real(expint(-x)) - Ei(x)) < 1e-5);
assert(abs(real(Li(x))-Ei(log(x))) < 1e-5);

% Liste des entiers examinés
X = 2:40;

% PI exact
for i = 1:length(X)
  PI(i) = length(primes(X(i)));
end

% 1ere image de la video
setfig;
plot(X, PI);
axis([min(X) max(X) 0 max(PI)+1.5]);
title(sprintf("PI approché"));
drawnow; refresh;

movie = VideoWriter("PIapproche.mp4");
movie.FrameRate = 5;
open(movie);
writeVideo(movie, getframe(gcf));

% Selon le nombre de termes de correction
for ng = 1:length(g)

  for i = 1:length(X)
    x      = X(i);
    L      = Ei(log(x));
    J      = J_riemann(x,g(1:ng));
    PIa(i) = pi_riemann(x,g(1:ng));

    if mod(x,50) == 0
      printf("Li(%5d)   = %.6f\n", x, L);
      printf("J(%5d)    = %.6f\n", x, J);
      printf("pi approche = %.6f\n", PIa(i));
      printf("pi exact    = %d\n", PI(i));
    end
  end

  % Image pour ce nombre de zeros de correction
  setfig;
  plot(X,PI, X,PIa);
  axis([min(X) max(X) 0 max(PI)+1.5]);
  title(sprintf("PI approché - %2d zeros", ng));
  drawnow; refresh;
  writeVideo(movie, getframe(gcf));
end
close(movie);
end

function setfig
  figure(3);
  xywh = get(gcf, "Position");
  xywh(3) = 800;
  xywh(4) = 600;
  set(gcf, "Position", xywh);
  set(gcf, "Resize", "off");
end