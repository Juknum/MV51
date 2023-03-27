% exemple de calcul formel avec trace

syms h x;
f = (h - x) * x * (x + h);
fp = diff(f, x);
d = solve(expand(fp));

q = simplify(subs(f, 'x', d));
r = eval(vpa(q));
pretty(eval(vpa(q)));

disp('pour voir le graphe, appuyez sur une touche');
pause;

X = -1:2/500:1;
Y = subs(subs(f, h, 1), x, X);
plot(X, Y);