%function [aopt, fopt] = mc(cnf_file, eval_budget)
function [fopt] = mc(cnf_file, eval_budget)
% [aopt, fopt] = mc(cnf_file, eval_budget)
%
% Monte Carlo MAX-3SAT solver
%
% Author: Johannes W. Kruisselbrink, Edgar Reehuis
% Last modified: September 7, 2011

	% Do you want online plotting? If not, set to false
	doplot = false;

	% Load CNF file and create fitness function handle
	cnf_expr = cnf_read(cnf_file);
	fitnessfct = @(a) evaluate_sat_expr(cnf_expr, a);
	n = size(cnf_expr, 1);

	% Initialize algorithm
	aopt = rand(n,1) > 0.5;
	fopt = fitnessfct(aopt);
  evals_used = 1;

	% Statistics
  if (doplot)
	  histf = zeros(1, eval_budget);
    histf(1:evals_used) = fopt;
  end

	% Main Monte-Carlo loop
	while (evals_used < eval_budget)

		% Generate random bitstring
		a = rand(n,1) > 0.5;
		fa = fitnessfct(a);

		% If the new bitstring is better than the one we have, replace the old one
		if (fa >= fopt)
			aopt = a;
			fopt = fa;
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
