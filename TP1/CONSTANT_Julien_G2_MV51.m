% rote.m

function [m] = rote(theta, M)
	theta_deg = theta * pi / 180;
	r = [cos(theta_deg), -sin(theta_deg); sin(theta_deg), cos(theta_deg)];
	m = r * M';
end

rote(90, [1, 0])

% symmetry_oi.m

function res = symmetry_oi(M)
	x = M(1);
	y = M(2);
	res = [x -y];
end

symmetry_oi([4, 3])

% trad_dec_qu.m 

function [alpha, beta] = trad_dec_qu(k)
	alpha = fix(k / 4);
	beta = mod(k, 4);
end

% trad_qu_dec.m 

function k = trad_qu_dec(ab)
	alpha = ab(1);
	beta = ab(2);
	k = 4 * alpha + beta;
end

k = 23;
[a, b] = trad_dec_qu(k);
disp([a, b]); % affiche [5, 3]

ab = [1, 2];
k = trad_qu_dec(ab);
disp(k); % affiche 6

% compose.m

function k = compose(sigma_i, sigma_j)
	alpha = mod((sigma_i(1) + sigma_j(1)), 2);
	beta = mod((-1)^(alpha_j) * beta_i + beta_j, 4);
	k = trad_qu_dec([alpha, beta]);
end

function res = tab(n)
	v = 0 : (2 * n - 1);
	comb = combvec(v, v);
	comb1 = comb(1, :);
	comb2 = comb(2, :);
	arr = arrayfun(@compose, comb1, comb2);
	res = reshape(arr, 2 * n, 2 * n);
end

tab(4)

% inverse.m

function j = inverse(i)
	t = tab(4);
	row_i = t(i + 1, :);
	j = find(row_i == 0) - 1;
end

tableau_d8 = tab(4)

% sous_groupe.m

function sous_groupe(vec)
	while true
		comb = combvec(vec, vec)
		comb1 = comb(1, :);
		comb2 = comb(2, :);

		arr = arrayfun(@compose, comb1, comb2);
		nvec = unique([vec arr]);

		if isequal(nvec, unique(vec))
			vec = nvec;
			break;
		end

		vec = nvec;

	end
	vec(end + 1; 8) = -1;
	res = vec;
end

% tous_les_sous_groupes.m

function res = tous_les_sous_groupes()
	elements = 0 : 2 * 4 - 1;
	max_card = 0;
	sub_groups_arr = [];
	index = 0;

	while max_card < length(elements) - 1:
		index = index + 1;
		vec = nchoosek(elements, index);
		[nvec, ~] = size(vec);

		for row_i = 1 : 1 : nvec
			row_curr = vec(row_i, :);
			ngroup = sous_groupe(row_curr);
			sub_groups_arr = [sub_groups_arr; ngroup];

			if sum(ngroup >= 0) > max_card
				max_card = sum(ngroup >= 0);
			end
		end
	end

	res = unique(sub_groups_arr, 'rows');
end

tous_les_sous_groupes()

% trad_dec_n.m 

function [alpha, beta] = trad_dec_n(k)
  alpha = fix(k / n);
  beta = mod(k, n);
end

% trad_n_dec.m 

function k = trad_n_dec(ab, n)
  alpha = ab(1);
  beta = ab(2);
  k = n * alpha + beta;
end

% compose_n.m

function k = compose_n(sigma_i, sigma_j, n)
  [alpha_i, beta_i] = trad_dec_n(sigma_i);
  [alpha_j, beta_j] = trad_dec_n(sigma_j);

  alpha = mod(alpha_i + alpha_j, n);
  beta = mod((-1)^(alpha_j) * beta_i + beta_j, n);

  k = trad_n_dec([alpha, beta], n);
end

% tab_n.m

function res = tab_n(n)
  v = 0 : (2 * n - 1);
  comb = combvec(v, v);
  comb1 = combinations(1, :);
  comb2 = combinations(2, :);

  arr = arrayfun(@compose_n, comb1, comb2, n * ones(1, length(comb1)));
  res = reshape(arr, 2 * n, 2 * n);
end

% inverse_n.m

function j = inverse_n(i, n)
  t = tab_n(n);
  row_i = t(i + 1, :);
  j = find(row_i == 0) - 1;
end

% sous_groupe_n.m

function res = sous_groupe_n(i, n)
  while true:
    comb = combvec(vec, vec);
    comb1 = comb(1, :);
    comb2 = comb(2, :);

    arr = arrayfun(@compose_n, comb1, comb2, n * ones(1, length(comb1)));
    nvec = unique([vec arr]);

    if isequal(nvec, unique(vec))
      vec = nvec;
      break;
    end

    vec = nvec;
  end
  vec(end + 1 : 2 * n) = -1;
  res = vec;
end

% tous_les_sous_groupes_n.m

function res = tous_les_sous_groupes_n(n)
  elements = 0 : 2 * n - 1;
  max_card = 0;
  sub_groups_arr = [];
  index = 0;

  while max_card < length(elements) - 1:
    index = index + 1;
    vec = nchoosek(elements, index);
    [nvec, ~] = size(vec);

    for row_i = 1 : 1 : nvec
      row_curr = vec(row_i, :);
      ngroup = sous_groupe_n(row_curr, n);
      sub_groups_arr = [sub_groups_arr; ngroup];

      if sum(ngroup >= 0) > max_card
        max_card = sum(ngroup >= 0);
      end
    end
  end

  res = unique(sub_groups_arr, 'rows');
end

M = [0, 0; 2, 2; 2, 1];
plot(M(:, 1), M(:, 2), -.5, -.5, 2.5, 2.5);

% apply_sigma.m

function res = apply_sigma(sigma_k, point, n)
	theta = 360 / n;
	nb_rot = mod(sigma_k, n);
	nb_sym = fix(sigma_k / n);
	npoint = point;

	for i = 1 : nb_sym
		npoint = symmetry_oi(npoint);
	end

	for i = 1 : nb_rot
		npoint = rote(theta, npoint);
	end

	res = npoint;
end

% apply_sigma_motif

function final = apply_sigma_motif(sigma_k, motif, n)
	nmotif = zeros(size(motif));

	for row_i = 1 : length(motif)
		nmotif(row_i, :) = apply_sigma(sigma_k, motif(row_i, :), n);
	end

	final = nmotif;
end

% creat_motif.m

function [] = creat_motif(motif, h, n)
	all_motif = [];

	clf
	hold on

	for i = 1 : length(h)
		nmotif = apply_sigma_motif(h(i), motif, n);
		all_motif = [all_motif; nmotif];
		plot(nmotif(:, 1), nmotif(:, 2));
	end

	axis([min(all_motif(:,1)) -.5 max(all_motif(:,1)) +.5 ... min(all_motif(:,2)) -.5 max(all_motif(:,2)) +.5])

	hold off
end

creat_motif(M, [0, 1, 2, 3, 4, 5, 6, 7], 4);

% triangle

T = [1 1; 1.5 2; 2 1; 1 1];
plot(T(:, 1), T(:, 2), .5, .5, 2.5, 2.5);