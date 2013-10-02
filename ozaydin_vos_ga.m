function [opt_tour, opt_tour_length] = ga_skeleton(tsp_instance, eval_budget)
% [opt_tour, opt_tour_length] = ga_skeleton(tsp_instance, eval_budget)
%
% Skeleton of a Genetic Algorithm for solving TSP instances
%
% Input:
% - tsp_instance        - string containing the TSP instance name
% - eval_budget         - the number of function evaluations available
%                         for the optimization algorithm (number of generations?)
%
% Output:
% - opt_tour            - a vector containing the best tour found
% - opt_tour_length     - the length of the best tour found
%
% Author: Johannes Kruisselbrink, Edgar Reehuis
% Last modified: February 4, 2011

	% Set true to do online statistics plotting
	doplot = false;

	% Retrieve the city coordinates, distance matrix, and number of cities
	[num_cities, coordinates, distance_matrix] = analyze_tsp(tsp_instance);
	%vpa(distance_matrix, 8)
	% Initialize static parameters
	lambda = 50; % amount of individuals
	pc = 0.8; % amount between 0 and 1, chance to crossover, or just copy the parent
	pm = 1; % Fixed amount of mutations
	q = 70; % amount of loops in select_tournament

	% Statistics data
	evalcount = 0;
	gencount = 0;
	opt_tour_length = Inf;
	opt_tour = NaN(num_cities,1);

	% Initialize population
	for i = 1:lambda

		% Each individual is a random permutation
		%r = ceil(rand(1) * num_cities)
		%P(i,:) = ceil(rand(1, num_cities) * num_cities); % Every city is random, and can be visited none, one or more times.
		%P(i,:) = mod([r:(r+num_cities-1)], num_cities) + 1 % 1 - 38, from random starting point.
		P(i,:) = randperm(num_cities); % random sequence of the cities
		
		% Evaluate the individual
		f(i) = evaluate_tour(distance_matrix, P(i,:));

		% Increase counter after each evaluation and update statistics
		evalcount = evalcount + 1;
		hist_length_best_so_far(evalcount) = opt_tour_length;
		if (f(i) < opt_tour_length)
			opt_tour_length = f(i);
			opt_tour = P(i,:);
		end

	end

	% Evolution loop
	while evalcount < eval_budget
		% Increase generation counter
		gencount = gencount + 1;

		% Generate new population Pnew (recombination, mutation)
		for i = 1:lambda

			% Select parent
			p1 = select_tournament(P, f, q);
			p2 = select_tournament(P, f, q);

			if (rand() < pc)
				[c1, c2] = pmx(p1, p2);
			else
				% No crossover, copy parent
				c1 = p1;
				c2 = p2;

			end

			% Apply mutation
			c1 = mutation(c1, num_cities, pm);

			Pnew(i,:) = c1;
		end

		% Replace old population by new population
		P = Pnew;

		% Evaluate each individual
		for i = 1:lambda

			% Evaluate the tour of each individual
			f(i) = evaluate_tour(distance_matrix, P(i,:));

			% Increase counter after each evaluation and update statistics
			evalcount = evalcount + 1;
			hist_length_best_so_far(evalcount) = opt_tour_length;
			if (f(i) < opt_tour_length)
				opt_tour_length = f(i);
				opt_tour = P(i,:);
			end

		end

		% Generation best statistics
		[hist_generation_opt_tour_length, generation_opt_tour_index] = min(f);
		generation_opt_tour = P(generation_opt_tour_index,:);
		hist_generation_opt_tour(gencount) = hist_generation_opt_tour_length;
			
		if (doplot)
			
			% Plot statistics
			clf

			subplot(2,2,1)
			plot(hist_length_best_so_far(1:evalcount))
			title('Length of best tour found so far')

			subplot(2,2,2)
			plot_tsp_tour(coordinates, opt_tour)
			title('Best tour found so far')

			subplot(2,2,3)
			plot(hist_generation_opt_tour(1:gencount))
			title('Length of best tour the current generation')

			subplot(2,2,4)
			plot_tsp_tour(coordinates, generation_opt_tour)
			title('Best tour of the current generation')

			drawnow()
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

	mu = length(P(:,1));

	Ptournament = [];
	Ftournament = [];
	for i = 1:q
		p = ceil(rand() * mu);
		Ptournament(i,:) = P(p,:);
		Ftournament(i) = f(p);
	end
	[maxF, maxIndex] = min(Ftournament);
	a = Ptournament(maxIndex,:);

end

function [c1, c2] = pmx(p1, p2)
	points = ceil(rand(2,1)*size(p1,2));
	if (points(1) > points(2))
		points = [points(2); points(1)];
	end

	c1(points(1):points(2)) = p1(points(1):points(2));
		
	for j = 1:size(p1,2)
		if j < points(1) || j > points(2)
			k = j;
			while ismember(p2(k),c1) == 1
				k = find(p1 == p2(k));
			end					
			c1(j) = p2(k);
		end
	end

	c2(points(1):points(2)) = p2(points(1):points(2));

	for j = 1:size(p2,2)
		if j < points(1) || j > points(2)
			k = j;
			while ismember(p1(k),c2) == 1
				k = find(p2 == p1(k));
			end					
			c2(j) = p1(k);
		end
	end
end

function s = mutation(s, num_cities, pm)
	for i=1:pm
		a = ceil(rand() * num_cities);
		b = ceil(rand() * num_cities);
		s([a b]) = s([b a]);
	end
end
