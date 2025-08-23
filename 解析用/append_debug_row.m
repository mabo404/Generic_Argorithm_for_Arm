function append_debug_row(dbgfile, g, tsec, invcnt, viol, best_obj)
    fid = fopen(dbgfile,'a');
    fprintf(fid,'%d,%.6f,%d,%d,%.8g,%.8g,%.8g,%.8g\n', g, tsec, invcnt, viol, best_obj);
    fclose(fid);
end