classdef Log
    methods (Static)
        function append_debug_row(dbgfile, g, tsec, invcnt, viol, best_obj)
            fid = fopen(dbgfile,'a');
            fprintf(fid,'%d,%.6f,%d,%d,%.8g,%.8g,%.8g,%.8g\n', g, tsec, invcnt, viol, best_obj);
            fclose(fid);
        end

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

        function save_generation_csv(logdir, g, vars, objs)
            try
                writematrix(vars, fullfile(logdir, sprintf('generation_%03d_vars.csv', g)));
                writematrix(objs, fullfile(logdir, sprintf('generation_%03d_objs.csv', g)));
            catch
            end
        end
    end
end