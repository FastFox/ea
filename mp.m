function [fopt] = mp(cnf_file, eval_budget)
%function [aopt, fopt] = mp(cnf_file, eval_budget)
% [aopt, fopt] = mc(cnf_file, eval_budget)
%
% Michiel Pepijn MAX-3SAT solver
%
% Author: Johannes W. Kruisselbrink, Edgar Reehuis
% Last modified: September 7, 2011

	% Do you want online plotting? If not, set to false
	doplot = true;
	plotbest = true;

	% Load CNF file and create fitness function handle
	cnf_expr = cnf_read(cnf_file);
	fitnessfct = @(a) evaluate_sat_expr(cnf_expr, a);
	n = size(cnf_expr, 1);

	% Initialize algorithm
	%aopt = rand(n,1) > 0.5;
	%fopt = fitnessfct(aopt);
  	%evals_used = 1;

	evals_used = 0;
	lambda = 50;
	pc = 0.5;
	pm = 1;
	q = 10; % amount of loops in select_tournament
	%P = zeros(lambda, n);
  k = 2; % Tournament size
	P = randn(lambda, n) > 0.5;
	Pnew = zeros(lambda, n);
		
	histf = zeros(1, eval_budget);

	% Initialize population
	for i = 1:lambda
		%P(i, :) = rand(n, 1) > 0.5;
		f(i) = fitnessfct(P(i, :)');
		if(i == 1 | f(i) > fopt | plotbest == false) % 1 == 1, so it will plot every result instead of the best
			aopt = P(i, :);
			fopt = f(i);
		end
		% Statistics
		if (doplot)
			histf(evals_used+1:evals_used+1) = fopt;
			subplot(2,1,1)
			plot(histf(1:evals_used+1))
			subplot(2,1,2)
			bar([1:n],aopt)
			xlim([1 n])
			drawnow();
		end
		evals_used = evals_used + 1;
	end

	for g = 2:ceil(eval_budget / lambda) - 1 % Initialization is the first generation

		for i = 1:lambda
				%p1 = select_tournament(P, f, q)
				%p1 = select_proportional_parent(P, f, lambda);
        p1 = select_tournament_parent(f, k);
        p2 = select_tournament_parent(f, k);
				%p2 = select_proportional_parent(P, f, lambda);
				%p2 = select_tournament(P, f, q);
				
				%c1 = mutation(p1, 1);
				%
				c1 = crossover(P(p1, :), P(p2, :));

				c11 = mutation(c1, 1, n);

				Pnew(i, :) = c11;
		end

		P = Pnew;

		for i = 1:lambda
			f(i) = fitnessfct(P(i, :)');
			if(f(i) > fopt | plotbest == false)
				aopt = P(i, :);
				fopt = f(i);
			end

			% Statistics maintenance and plotting
			if (doplot)
				histf(evals_used+1:evals_used+1) = fopt;
				subplot(2,1,1)
				plot(histf(1:evals_used+1))
				subplot(2,1,2)
				bar([1:n],aopt)
				xlim([1 n])
				drawnow();
			end
			
			evals_used = evals_used + 1;

		end

	end
	
	% Main Monte-Carlo loop
	%while (evals_used < eval_budget)
		%Pnew =  %initialize

		% Generate random bitstring
		%a = rand(n,1) > 0.5;
	%	a = mutation(aopt, 1);
	%	fa = fitnessfct(a);

		% If the new bitstring is better than the one we have, replace the old one
	%	if (fa >= fopt)
	%		aopt = a;
	%		fopt = fa;
	%	end

		% Statistics maintenance and plotting
	%	if (doplot)
	%		histf(evals_used+1:evals_used+1) = fopt;
	%		subplot(2,1,1)
	%		plot(histf(1:evals_used+1))
	%		subplot(2,1,2)
	%		bar([1:n],aopt)
	%		xlim([1 n])
	%		drawnow();
	%	end

   % 		evals_used = evals_used + 1;
%	end
	
end

function s = mutation(s, pm, n)
	for i=1:pm
		a = ceil(rand() * n);
		b = ceil(rand() * n);
		s([a b]) = s([b a]);
	end
end


function r = select_proportional_parent(P, f, lambda)
	
	% Total fitness
	total = sum(f);

	lots = zeros(100, 1);
	j = 1; % Lot index

	% For every solution
	for i = 1:lambda
		for k = j:j + ceil(f(i) / total * 100) % from lot index loop %-times.
			lots(k) = i;
			j = j + 1;
		end
	end

	r = lots(ceil(99 * rand() + 1));
end

function bi = select_tournament_parent(f, k)
  b = -1;
  bi = -1; % Best Index
  n = length(f);
  for i = 1:k
    random = ceil(rand() * n);
    if (f(random) > b)
      bi = random;
      b = f(random);
    end
  end
end

% Uniform crossover
function offspring = crossover(p1, p2)
	offspring = p1;
	for i = 1:length(p1)
		if rand() > 0.5
			offspring(i) = p2(i);
		end
	end
end

function a = select_tournament(P, f, q)
% a = select_tournament(P, f, q)
%
% Select one individual of population P with fitness
% values f using proportional selection.
%
% Author: Johannes W. Kruisselbrink
% Last modified: October 21, 2009

	mu = length(P(:, 1));

	Ptournament = [];
	Ftournament = [];
	for i = 1:q
		p = ceil(rand() * mu);
		Ptournament(i, :) = P(p, :);
		Ftournament(i) = f(p);
	end
	[maxF, maxIndex] = min(Ftournament);
	a = Ptournament(maxIndex, :);

end
