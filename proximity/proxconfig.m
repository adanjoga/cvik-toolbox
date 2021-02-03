function fun = proxconfig(opt)
if strcmpi(opt,'euc')
    % Euclidean distance
    fun = @eucdist;
elseif strcmpi(opt,'neuc')
    % Normalized Euclidean distance
    fun = @neucdist; 
elseif strcmpi(opt,'cos')
    % Cosine similarity
    fun = @cosdist;
elseif strcmpi(opt,'pcorr')
    % Pearson's correlation coefficient
    fun = @pcorr;
elseif strcmpi(opt,'scorr')
    % Spearman's correlation coefficient
    fun = @scorr; 
elseif strcmpi(opt,'lap')
    % Laplacian distance
    fun = @lapdist; 
else
    error('Unknown proximity measure');
end