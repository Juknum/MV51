% exemple de calcul formel avec trace.
syms a b x ;
f  = (a - x )*( x - b );
fp = diff(f, x);
d  = solve(fp);

pretty(simplify(subs(f,x,d)));

disp('pour voir le graphe, appuyez sur une touche');
pause;

X = -1:2/500:1;
Y = subs(subs(f, {a, b}, {-1, 1}), x ,X);
plot(X,Y);