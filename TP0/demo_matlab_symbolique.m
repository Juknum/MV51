disp('exemple 1');
syms x a b;
y = x^2 + a * x + b;

disp(y);
disp(diff(y, x));
disp(diff(y, a));

disp('exemple 2');
syms x h;
limit((sin(x + h) - sin(x)) / h, h, 0)

disp('exemple 3');
disp(int(sin(x) / x, 0, inf));