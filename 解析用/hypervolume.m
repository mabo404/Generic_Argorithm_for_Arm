function hv = hypervolume(F, ref_point)
    [N, M] = size(F);

    if M ~= 2 && M ~= 3
        error('hypervolume関数は2次元または3次元のみ対応しています');
    end

    F = sortrows(F);

    if M == 2
        hv = 0;
        prev_f1 = ref_point(1);
        for i = N:-1:1
            width = prev_f1 - F(i,1);
            height = ref_point(2) - F(i,2);
            hv = hv + width * height;
            prev_f1 = F(i,1);
        end
    else
        hv = 0;
        for i = 1:N
            vol = prod(ref_point - F(i,:));
            hv = hv + vol;
        end
    end
end
