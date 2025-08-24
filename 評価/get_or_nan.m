function v = get_or_nan(S, name)
%GET_OR_NAN 構造体Sに name があれば値、なければ NaN
    if isfield(S, name), v = S.(name); else, v = NaN; end
end
