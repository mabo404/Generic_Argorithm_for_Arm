function must_have(S, fieldname, who)
%MUST_HAVE セレクタ設定の必須フィールド存在チェック
    if ~isfield(S, fieldname)
        error('%s: %s フィールドが必要です。', who, fieldname);
    end
end
