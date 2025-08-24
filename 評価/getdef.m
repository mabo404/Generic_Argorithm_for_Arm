function v = getdef(S, name, default)
%GETDEF 構造体Sのフィールド name があれば値、なければ default
    if isfield(S, name), v = S.(name); else, v = default; end
end
