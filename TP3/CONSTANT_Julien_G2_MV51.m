disp('1) a. affiche_h([1 2 3 4])');
affiche_h([1 2 3 4]);
disp(' ');
disp('1) b. somme_h([1 2 3 4], [5 6 7 8])');disp('ans =');
display(somme_h([1 2 3 4], [5 6 7 8]));
disp(' ');
disp('1) c. produit_h([1 2 3 4], [5 6 7 8])');disp('ans =');
display(produit_h([1 2 3 4], [5 6 7 8]));
disp(' ');
disp('1) d.');
disp('>> conjugue_h([1 2 3 4])');disp('ans =');
display(conjugue_h([1 2 3 4]));
disp('>> norme_h([1 2 3 4])');disp('ans =');
display(norme_h([1 2 3 4]));
disp('>> inverse_h([1 2 3 4])');disp('ans =');
display(inverse_h([1 2 3 4]));

disp(' ');
disp(' ');
disp('2) a. polaire_h([1 2 3 4])');
polaire = polaire_h([1 2 3 4]);
display(polaire);

disp(' ');
disp('b. algebr_h(polaire(1), polaire(2), polaire(3:5))');disp('ans =');
display(algebr_h(polaire(1), polaire(2), polaire(3:5)));

disp(' ');
disp(' ');
disp('3)');
display(proof());

disp(' ');
disp(' ');
disp('4) a. matrice_rotation([1 2 3], pi/6)');disp('ans =');
display(matrice_rotation([1 2 3], pi/6));
disp(' ');
disp('b. h_rot([1, 2, 3, 4])');
hrot = h_rot([1, 2, 3, 4]);
display(hrot);

disp(' ');
disp('c. rot_h(hrot)');disp('ans =');
display(rot_h(hrot));

% 1)
function [] = affiche_h(h)
	fprintf('q = %d + %di + %dj + %dk\n', h(1), h(2), h(3), h(4));
end

function [h] = somme_h(h1, h2)
	h = h1 + h2;
end

function [h] = produit_h(h1, h2)
	v1 = h1(2:4);
	v2 = h2(2:4);
	v = h1(1) * v2 + h2(1) * v1 + cross(v1, v2);
	h(1) = h1(1) * h2(1) - dot(v1, v2);
	h(2:4) = v;
end

function [H] = conjugue_h(h)
	H = [h(1), - h(2:4)];
end

function [H] = norme_h(h)
	H = sqrt(h(1)^2 + h(2)^2 + h(3)^2 + h(4)^2);
end

function [H] = inverse_h(h)
	if ~isequal(h, [0, 0, 0, 0])
		H = conjugue_h(h) / norme_h(h)^2;
	else
		disp('Erreur : division par 0');
	end
end

% 2)
function p = polaire_h(h)
	r = norme_h(h);
	theta = abs(acos((h(1) / r)));
	vect = h(2:4);
	i = vect / (r * sin(theta));
	p = [r, theta, i];
end

function [H] = algebr_h(r, theta, i)
	h = [0, 0, 0, 0];
    h(1) = r * cos(theta);
    h(2:4) = r * i * sin(theta);
    H = h;
end

% 3)
function [res] = proof()
	M = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];

	res = true;

	for u = 1:1:3
		for v = 1:1:3
			for w = 1:1:3
				if (produit_h(produit_h(M(u * 4 + 1 : u * 4 + 4), M(v * 4 : v * 4 + 3)), M(w * 4 + 1 : w * 4 + 4)) ~= produit_h(M(u * 4 + 1 : u * 4 + 4), produit_h(M(v * 4 + 1 : v * 4 + 4), M(w * 4 + 1 : w * 4 + 4))))
					res = false;
				end
			end
		end
	end

end

% 4)
function [I] = matrice_rotation(n, theta)
	I = [
		n(1)^2 * (1 - cos(theta)) + cos(theta), n(1) * n(2) * (1 - cos(theta)) - n(2) * sin(theta), n(1) * n(3) * (1 - cos(theta)) + n(2) * sin(theta);
		n(1) * n(2) * (1 - cos(theta)) + n(1) * sin(theta), n(2)^2 * (1 - cos(theta)) + cos(theta), n(2) * n(3) * (1 - cos(theta)) - n(1) * sin(theta);
		n(1) * n(3) * (1 - cos(theta)) - n(2) * sin(theta), n(2) * n(3) * (1 - cos(theta)) + n(1) * sin(theta), n(3)^2 * (1 - cos(theta)) + cos(theta)
	];
end

function [I] = h_rot(q)
	s = q(1);
	x = q(2);
	y = q(3);
	z = q(4);
    
    r00 = 1 - 2*y^2 - 2*z^2;
    r01 = 2*x * y - 2*s * z;
    r02 = 2*x * z + 2*s * y;
    r10 = 2*x * y + 2*s * z;
    r11 = 1 - 2*x^2 - 2*z^2;
    r12 = 2*y * z - 2*s * x;
    r20 = 2*x * z - 2*s * y;
    r21 = 2*y * z + 2*s * x;
    r22 = 1 - 2*x^2 - 2*y^2;

	I = [r00 r01 r02; r10 r11 r12; r20 r21 r22];
end

function [q] = rot_h(L)
	r00 = L(1, 1);
	r01 = L(1, 2);
	r02 = L(1, 3);
	r10 = L(2, 1);
	r11 = L(2, 2);
	r12 = L(2, 3);
	r20 = L(3, 1);
	r21 = L(3, 2);
	r22 = L(3, 3);

	T = r00 + r11 + r22;
	Tmax = max([r00, r11, r22]);

	if (T > 0)
		s = sqrt(T + 1) / 2;
		x = (r21 - r12) / (4 * s);
		y = (r02 - r20) / (4 * s);
		z = (r10 - r01) / (4 * s);
	
	elseif (Tmax == r00)
		x = sqrt(r00 - r11 - r22 + 1) / 2;
		s = -(r12 - r21) / (4 * x);
		y = (r01 + r10) / (4 * x);
		z = (r02 + r20) / (4 * x);

	elseif (Tmax == r11)
		y = sqrt(r11 - r00 - r22 + 1) / 2;
		s = -(r20 - r02) / (4 * y);
		x = (r01 + r10) / (4 * y);
		z = (r12 + r21) / (4 * y);

	elseif (Tmax == r22)
		z = sqrt(r22 - r00 - r11 + 1) / 2;
		s = -(r01 - r10) / (4 * z);
		x = (r02 + r20) / (4 * z);
		y = (r12 + r21) / (4 * z);
    end
    
    q = [s x y z];
end