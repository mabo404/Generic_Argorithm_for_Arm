classdef SE3
    methods (Static)
        function T = transl(x,y,z) % 平行移動行列Tを返す。引数はx,y,zの平行移動量（想定は実数スカラー）。
            T = [1 0 0 x;
                 0 1 0 y;
                 0 0 1 z;
                 0 0 0 1]; % 4×4の同次変換行列を直接構築。回転は単位、右端列が並進ベクトル、最下行は同次座標の規約[0 0 0 1]。
        end

        function R = rotz(th) % z軸周りの回転（ラジアン）を表す4×4同次行列Rを返す。
            R = [cos(th) -sin(th) 0 0;
                 sin(th)  cos(th) 0 0;
                 0        0       1 0;
                 0        0       0 1]; % 上左3×3がz回りの回転行列、右端列はゼロ、最下行[0 0 0 1]。thがベクトルの場合はサイズ不整合になる設計（=スカラー前提）。
        end

        function R = rotx(th) % x軸周りの回転（ラジアン）を返す。戻り値は4×4同次行列。
            R = [1 0        0       0;
                 0 cos(th) -sin(th) 0;
                 0 sin(th)  cos(th) 0;
                 0 0        0       1]; % 上左3×3がx回りの回転行列。thは（暗黙に）スカラー前提。
        end
    end
end