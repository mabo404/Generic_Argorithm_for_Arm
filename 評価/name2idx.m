function j = name2idx(nameOrIdx, obj_names)
%NAME2IDX  目的名/番号を列番号に解決
    if isnumeric(nameOrIdx)
        j = nameOrIdx;
    else
        j = find(strcmp(obj_names, string(nameOrIdx)), 1);
    end
    if isempty(j)
        error('目的名 "%s" が見つかりません。', string(nameOrIdx));
    end
end
