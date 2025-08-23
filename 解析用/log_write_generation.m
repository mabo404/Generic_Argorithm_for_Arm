function log_write_generation(gen, objs, vars, logdir, meta)
    if ~exist(logdir,'dir'), mkdir(logdir); end

    f1 = fullfile(logdir, sprintf('gen_%04d_objs.csv', gen));
    f2 = fullfile(logdir, sprintf('gen_%04d_vars.csv', gen));
    writematrix(objs, f1);
    writematrix(vars, f2);

    dbg = fullfile(logdir, 'debug_log.csv');
    T = table( ...
        gen, size(objs,1), size(objs,2), size(vars,2), datetime('now'), ...
        'VariableNames', {'gen','numIndividuals','numObjs','numVars','timestamp'} ...
    );

    if nargin>=5 && ~isempty(meta) && ~isempty(fieldnames(meta))
        fn = fieldnames(meta);
        for k = 1:numel(fn)
            try
                T.(fn{k}) = meta.(fn{k});
            catch
                T.(fn{k}) = string(jsonencode(meta.(fn{k})));
            end
        end
    end

    if ~isfile(dbg)
        writetable(T, dbg);
    else
        writetable(T, dbg, 'WriteMode','Append');
    end
end
