a = [1 3 5 ; 1 2 4 ; 0 5 1 ];
b = [22; 17 ; 13];

tic;
for k = 1:1000
    x = inv(a)*b; 
end

t = toc;
disp(t);

tic;
for k = 1:1000 
    x = a\b;
end

t = toc;
disp(t);

% Tel que prévu, la fonction inv(a)*b est plus lente que a\b