classdef Util
    methods (Static)
        % function tf = util_isfinite(x, type) % エントリ関数。入力xの内容が有限値かを型に応じて検査して、論理値tfを返す。第2引数 type は任意。
        %     if nargin<2 || isempty(type) || strcmpi(type,'auto') % 自動判定フラグ。第2引数がない/空/'auto' の場合に型推定を行う。strcmpi は大文字小文字を無視。
        %         if isa(x,'quaternion') % quaternionオブジェクトなら優先的にクォータニオン扱いへ。MATLABの姿勢表現オブジェクトに対応。
        %             type = 'quat'; % 型フラグをクォータニオンに。以降のswitchで'quat'分岐へ誘導。
        %         elseif isnumeric(x) && isvector(x) && numel(x)==3 % 数値の3要素ベクトルなら3Dベクトル（位置等）とみなす。
        %             type = 'vec3'; % 型フラグを'vec3'に設定。
        %         elseif isnumeric(x) && isvector(x) && numel(x)==4 % 数値の4要素ベクトルなら（w,x,y,z）形式のクォータニオン配列とみなす。
        %             type = 'quat'; % 型フラグを'quat'に設定。
        %         elseif isnumeric(x) % それ以外の数値は一般の数値行列として扱う。
        %             type = 'mat'; % 型フラグを'mat'に設定。
        %         else % 上記のどれにも当てはまらない場合（文字列・構造体など）。
        %             type = 'other'; % 型フラグを'other'に。後段のotherwiseで汎用チェックを適用。
        %         end
        %     end
        % 
        %     % 2) 型別チェック（空配列は偽）
        %     switch lower(type) % 型フラグを小文字化して分岐。string/charの混在にも概ね対応。
        %         case 'vec3' % 3要素ベクトルの検査に入る。
        %             tf = isnumeric(x) && isvector(x) && numel(x)==3 && ~isempty(x) && all(isfinite(x(:))) && all(isreal(x(:))); % 数値・ベクトル・要素3・空でない・全要素が有限・全要素が実数、をすべて満たすとtrue。設計上、複素数やNaN/Infを明確に排除。
        %         case 'quat' % クォータニオンの検査に入る。数値4要素かquaternionオブジェクトのどちらにも対応。
        %             if isa(x,'quaternion') % quaternionオブジェクトかどうかで分ける。
        %                 tf = ~isempty(x) && all(isfinite(compact(x))); % compact(x)で(w,x,y,z)を数値配列として取り出し、全要素が有限かを検査。空は偽。isrealはquaternionでは暗黙に満たされる想定。
        %             else % こちらは「数値4要素ベクトル」として渡されたクォータニオン。
        %                 tf = isnumeric(x) && isvector(x) && numel(x)==4 && ~isempty(x) && all(isfinite(x(:))) && all(isreal(x(:))); % 数値・ベクトル・要素4・空でない・全要素有限・全要素実数をすべて満たせばtrue。
        %             end
        %         case {'mat','matrix'} % 一般の数値行列。'mat'と'matrix'の同義扱い。
        %             tf = isnumeric(x) && ismatrix(x) && ~isempty(x) && all(isfinite(x(:))) && all(isreal(x(:))); % 数値・行列・空でない・全要素有限・全要素実数ならtrue。ベクトルもismatrixを満たすが、先のswitchでより限定的な'vec3'/'quat'が先に評価される。
        %         otherwise % それ以外の型（構造体・文字列・cell など）は汎用の厳しめチェックへ。
        %             tf = isnumeric(x) && ~isempty(x) && all(isfinite(x(:))); % 数値かつ空でなく、すべて有限ならtrue。isrealを課していないので、複素の有限値は通る（ここは「その他」扱いなので緩め/用途次第で変更可）。
        %     end
        % end
    end
end