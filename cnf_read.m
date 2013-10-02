function [cnf_formula] = cnf_read(cnf_file, varargin)

	fid = fopen(cnf_file);

	if fid < 0
		err_str = ['Error opening file ', cnf_file];
		error(err_str);
	end
	k = 1;

	while ~feof(fid)
		tline = fgetl(fid);
		tline = strtrim(tline);
		if isempty(tline)
			continue;
		end
		if (tline(1) == 'p' || tline(1) == 'P')
			st = regexpi(tline, 'cnf');
			dim = str2num(tline(st+3:end));
			cnf_formula = sparse(dim(1), dim(2));
			continue;
		end
		if (tline(1) ~= 'c') % discard comments
			x = str2num(tline);
			for i=1:length(x)-1
				v = x(i);
				if v > 0
					cnf_formula(v,k) = 1;
				elseif v < 0
					cnf_formula(-v,k) = -1;
				end
			end
			k = k + 1;
		end
	end

	fclose(fid);

end
