function y = f1(x)
% 5.1 fonctions et les M-fichiers -- Exo 1
%   L'intitulé parle de la fonction 'escalier' mais aucune trace de
%   cette dernière en faisant 'lookfor' ni dans l'intitulé.

y = zeros(size(x));
y(x < 0) = 0;
y(x >= 0 & x <= 2) = 1;
y(x > 2 & x < 18) = 2;
y(x >= 18) = 3;
end

