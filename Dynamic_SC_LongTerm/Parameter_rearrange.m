function [X,outputVar] = Parameter_rearrange(Xtilde,X,AmbiguitiesInfo)
    
    if size(Xtilde,1) ~= size(X,1)
        error('Input data must have identical number of rows')
    end
    
    if nargin<3 || ~isfield(AmbiguitiesInfo,'permutation') || ~isfield(AmbiguitiesInfo,'scaling')
        [permutation,scaling,~,~] = Unifying_Ambiguities(Xtilde,X);
        outputVar.permutation = permutation;
        outputVar.scaling = scaling;
    else
        permutation = AmbiguitiesInfo.permutation;
        scaling = AmbiguitiesInfo.scaling;
    end

    X_permute_removal = X*permutation;
    X_ambiguities_removal = X*permutation*scaling;

    all_zero_columns = find(all(X_permute_removal == 0, 1));

    X_permute_removal(:,all_zero_columns) = [];
    X_ambiguities_removal(:,all_zero_columns) = [];


    idx_leftover = ismember(X',X_permute_removal','rows');

    remaining_columns = find(idx_leftover == 0);

    X_remaining = X(:, remaining_columns); 
    X = [X_ambiguities_removal,X_remaining];

    [~,~,outputVar.sim,~] = Unifying_Ambiguities(Xtilde,X);

end