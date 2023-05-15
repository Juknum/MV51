disp(">> homothetie(3)");disp("ans =");
display(homothetie(3));
disp(">> rotation(pi/2)");disp("ans =");
display(rotation(pi/2));
disp(">> translation([2, 3])");disp("ans =");
display(translation([2, 3]));

% Décommentez les lignes suivantes pour tester la fonction fract2
% (n'oubliez pas de commmetez les lignes de fract3 pour éviter les conflits dans le graphique)

% disp(">> fract2([0, 0; 0, 1], 10, pi/6, .8)");
% fract2([0, 0; 0, 1], 10, pi/6, .8);

% disp(">> fract2([0, .5; 0, .7], 10, pi/6, .8)");
% fract2([0, .5; 0, .7], 10, pi/6, .8);

disp(">> homothetie(3)");disp("ans =");
display(homothetie3(3));
disp(">> rotation3(pi/2, pi/2, pi/2)");disp("ans =");
display(rotation3(pi/2, pi/2, pi/2));
disp(">> translation3([1, 2, 3])");disp("ans =");
display(translation3([1, 2, 3]));

disp(">> fract3([0, 0; 0, 0; 0, 1], 6, pi/6, .8)");
fract3([0, 0; 0, 0; 0, 1], 6, pi/6, .8);

% Renvoie une matrice de transformation homothétique 
% en fonction d'un facteur d'échelle "k".
function H = homothetie(k)
    % Initialise une matrice identité
    H = eye(3);
    % Applique le facteur d'échelle sur les dimensions x et y
    H(1 : dim - 1, 1 : dim - 1) = H(1 : dim - 1, 1 : dim - 1) * k;
end

% Renvoie une matrice de transformation de rotation 
% en fonction d'un angle alpha.
function R = rotation(alpha)
    R = [
        [cos(alpha), - sin(alpha), 0];
        [sin(alpha), cos(alpha), 0];
        [0, 0, 1]
    ];
end

% Applique une translation sur un vecteur de translation "a" pour l'adapter aux
% dimensions de la matrice
function T = translation(a)
    % Initialise une matrice identité
    T = eye(3);
    % Transpose le vecteur de translation pour l'adapter aux 
    % dimensions de la matrice
    a = transpose(a);

    % Affecte les composantes x et y du vecteur de translation aux 
    % éléments correspondants de la matrice
    T(1 : 2, 3) = a;
end

% Génère une fractale à partir d'une matrice M0, avec une 
% certaine profondeur (depth), un angle de rotation (alpha) et un facteur d'échelle (h).
function fract2(M0, depth, alpha, h)
    o = M0(:, 1); % Point initial de la fractale
    a = M0(:, 2); % Point final de la fractale
    t = a - o; % Vecteur de translation entre les deux points

    % Trace une ligne entre les deux points
    line([o(1); a(1)], [o(2); a(2)]);

    if depth > 0
        % Réduit la profondeur de 1 à chaque récursion
        depth = depth - 1;

        % Génère les transformations pour la première récursion
        g1 = translation(t) * rotation(alpha) * homothetie(h) * [a(1) - o(1); a(2) - o(2); 1];
        M_1 = [a(1), o(1) + g1(1); a(2), o(2) + g1(2)];
        fract2(M_1, depth, alpha, h);

        % Génère les transformations pour la deuxième récursion
        g2 = translation(t) * rotation(- alpha) * homothetie(h) * [a(1) - o(1); a(2) - o(2); 1];
        M_2 = [a(1), o(1) + g2(1); a(2), o(2) + g2(2)];
        fract2(M_2, depth, alpha, h);
    end
end

% Renvoie une matrice de transformation homothétique
% en fonction d'un facteur d'échelle "k". (3D)
function H = homothetie3(k)
    H = eye(4);
    % Applique le facteur d'échelle sur les dimensions x, y et z
    H(1 : 3, 1 : 3) = H(1 : 3, 1 : 3) * k;
end

% Renvoie une matrice de transformation de rotation 
% en fonction d'un angle alpha. (3D)
function R = rotation3(alpha_x, alpha_y, alpha_z)
    Rx = [
        1, 0, 0, 0;
        0, cos(alpha_x), - sin(alpha_x), 0;
        0, sin(alpha_x), cos(alpha_x), 0;
        0, 0, 0, 1
    ];

    Ry = [
        cos(alpha_y), 0, sin(alpha_y), 0;
        0, 1, 0, 0;
        - sin(alpha_y), 0, cos(alpha_y), 0;
        0, 0, 0, 1
    ];

    Rz = [
        cos(alpha_z), - sin(alpha_z), 0, 0;
        sin(alpha_z), cos(alpha_z), 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1
    ];

    R = Rx * Ry * Rz;
end

% Applique une translation sur un vecteur de translation "a" pour l'adapter aux
% dimensions de la matrice (3D)
function T = translation3(a)
    T = eye(4);
    a = transpose(a);
    
    % Affecte les composantes x, y et z du vecteur de translation aux 
    % éléments correspondants de la matrice
    T(1 : 3, 4) = a;
end

function fract3(M0, depth, alpha, h)
    o = M0(:, 1); % Point initial de la fractale
    a = M0(:, 2); % Point final de la fractale
    t = a - o; % Vecteur de translation entre les deux points
    
    % Trace une ligne entre les deux points
    line([o(1); a(1)], [o(2); a(2)], [o(3); a(3)]);
    view(3);

    if depth > 0
        depth = depth - 1;

        % Génère les transformations pour la première récursion
        g1 = translation3(t) * rotation3(alpha, alpha, 0) * homothetie3(h) * [a(1) - o(1); a(2) - o(2); a(3) - o(3); 1];
        M_1 = [a(1), o(1) + g1(1); a(2), o(2) + g1(2); a(3), o(3) + g1(3)];
        fract3(M_1, depth, alpha, h);

        % Génère les transformations pour la deuxième récursion
        g2 = translation3(t) * rotation3(- alpha, alpha, 0) * homothetie3(h) * [a(1) - o(1); a(2) - o(2); a(3) - o(3); 1];
        M_2 = [a(1), o(1) + g2(1); a(2), o(2) + g2(2); a(3), o(3) + g2(3)];
        fract3(M_2, depth, alpha, h);

        % Génère les transformations pour la troisième récursion
        g3 = translation3(t) * rotation3(alpha, - alpha, 0) * homothetie3(h) * [a(1) - o(1); a(2) - o(2); a(3) - o(3); 1];
        M_3 = [a(1), o(1) + g3(1); a(2), o(2) + g3(2); a(3), o(3) + g3(3)];
        fract3(M_3, depth, alpha, h);

        % Génère les transformations pour la quatrième récursion
        g4 = translation3(t) * rotation3(- alpha, - alpha, 0) * homothetie3(h) * [a(1) - o(1); a(2) - o(2); a(3) - o(3); 1];
        M_4 = [a(1), o(1) + g4(1); a(2), o(2) + g4(2); a(3), o(3) + g4(3)];
        fract3(M_4, depth, alpha, h);

        % On peut ensuite ajouter davantage de récursions...
    end
end