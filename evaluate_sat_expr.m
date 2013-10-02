function f = evaluate_sat_expr(expr, a)
% f = evaluate_sat_expr(expr, a)
%
% Given a CNF expression and candidate solution, this function computes the number 
% of clauses that are satisfied
%
% Author: Johannes Kruisselbrink, Edgar Reehuis
% Last modified: November 3, 2011

	a = double(a);
	a(a==0) = -1;
	[vars, terms] = size(expr);
	A = repmat(a,1,terms);
	AE = full(A .* expr);
	f = full(sum(max(AE,[],1)));


	if (exist('statistics') == 2) statistics(f); end
end
