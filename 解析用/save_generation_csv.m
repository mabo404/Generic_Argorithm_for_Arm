function save_generation_csv(logdir, g, vars, objs)
    try
        writematrix(vars, fullfile(logdir, sprintf('generation_%03d_vars.csv', g)));
        writematrix(objs, fullfile(logdir, sprintf('generation_%03d_objs.csv', g)));
    catch
    end
end