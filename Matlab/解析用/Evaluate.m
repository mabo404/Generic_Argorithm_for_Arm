function J = Evaluate(pop)
    %% 1) 入出力と定数定義
    INVALID_FITNESS = 1e6; % 失敗時に返す共通ペナルティ値（J全列に適用）
    SAFETY_THRESHOLD = 1.5; % トルク安全率の下限（これ未満は罰則）
    COLLISION_THRESH = 0.1; % 自己干渉の最小許容距離[m]（下回れば罰則）
    MAX_SAFETY = 1e6; % 安全率の上限クリップ（極端値/ゼロ割対策）
    invalid_row = repmat(INVALID_FITNESS, 1, 4); % 失敗時用の4目的共通ベクトル[J_tau,J_pos,J_ori,J_sc]

    %% 2) 機体設定（質量・最大トルク）
    JOINT_BASIC  = [0.577296, 0.530815, 0.530815, 0.136257, 0.070711, 0.136257]; % 各関節の基本重量[kg]（J1→J6。モーター重量は別途加算）
    MOTOR_WEIGHT = [5.3, 5.3, 5.3, 0.0565, 0.0565, 0.0565]; % 採用モーターの重量[kg]（J1→J6）。m_joint = JOINT_BASIC + MOTOR_WEIGHT で使用
    if ~isempty(JOINT_BASIC) && ~isempty(MOTOR_WEIGHT) % 値が設定されている場合のみ合成を試みる
        if isvector(JOINT_BASIC) && numel(JOINT_BASIC)==6 && all(isfinite(JOINT_BASIC)) && ... % 入力形式チェック：6要素の有限ベクトルか
           isvector(MOTOR_WEIGHT) && numel(MOTOR_WEIGHT)==6 && all(isfinite(MOTOR_WEIGHT)) % 〃（MOTOR_WEIGHT側）
            m_joint = JOINT_BASIC(:) + MOTOR_WEIGHT(:); % 基本重量とモーター重量を合算して上書き
        else  % 形式が不正な場合
            warning('JOINT_BASIC / MOTOR_WEIGHT の形式が不正です（6×1の有限値ベクトルにしてください）。既定の m_joint を使用します。'); % フォーマット不正を警告（処理は継続）
        end % 入力形式チェックの分岐終了
    end % 合成ブロック終了（未設定時は既定の m_joint を使用）

    TAU_MAX_OVERRIDE = [11.768; 11.768; 11.768; 1.862; 1.862; 1.862]; % 採用モーターの最大トルク[N·m]（J1→J6）。指定時は tau_max をこの値で上書き
    if ~isempty(TAU_MAX_OVERRIDE) % 上書きトルクが指定されていれば適用を開始
        if isvector(TAU_MAX_OVERRIDE) && numel(TAU_MAX_OVERRIDE)==6 && all(isfinite(TAU_MAX_OVERRIDE)) % 入力形式チェック：6要素の有限値ベクトルか（行/列は不問）
            tau_max = TAU_MAX_OVERRIDE(:); % 列ベクトル化して最大トルクを上書き適用[N·m]（J1→J6）
        else % 形式が不正だった場合のフォールバック分岐
            warning('TAU_MAX_OVERRIDE の形式が不正です（6×1の有限値ベクトルにしてください）。既定の tau_max を使用します。'); % フォーマット不正を通知（上書きは行わず既定値を使用）
        end % 入力形式チェック（TAU_MAX_OVERRIDE）の分岐終了
    end % tau_max 上書き処理ブロックの終了

    %% 3) 入力チェックと出力配列の確保
    [n,~] = size(pop); % 個体数 n を取得（列数はここでは未使用）
    if size(pop,2) < 5 % 設計変数 L1..L5 が存在するかをチェック（5列未満なら不可）
        error('Evaluate:InvalidPop','pop には少なくとも L1..L5 の5列が必要です。'); % 必須列不足を即時エラー
    end % チェックの終了
    J = zeros(n,4); % 出力行列 J を [n×4] で初期化（[J_tau J_pos J_ori J_sc]）

    %% 4) 時間軸と参照軌道（位置/姿勢）の事前計算
    dt = 0.02; % サンプリング周期 [s]
    T_total = 10; % 解析の総時間 [s]
    T = 0:dt:T_total; % 時間ベクトル（0 〜 T_total を dt 刻み）
    Xe = [-0.3; -0.3; 0.2]; % 終点位置 [m]（x, y, z）
    r = T / T_total; % 正規化時間（0〜1）
    S = 10.*(r.^3) - 15.*(r.^4) + 6.*(r.^5); % 最小ジャークの時刻関数 S(r)
    Qstart = quaternion([-45, 0, 90], 'eulerd', 'ZYX', 'frame');
    QE = quaternion([-45 0 90], 'eulerd', 'ZYX', 'frame');

    % (A) 最短経路化：終点の符号を合わせる（antipodal対策）
    dot_q = sum( compact(Qstart) .* compact(QE) );
    if dot_q < 0
        QE = quaternion(-compact(QE));
        dot_q = -dot_q;
    end

    % (B) 退避策：slerpの特異域（|dot|≈1）に近いときの処理
    if dot_q > 0.9995
        % ほぼ同一姿勢：線形補間→正規化（nlerp相当）
        q1 = compact(Qstart);           % 1×4
        q2 = compact(QE);               % 1×4
        W1 = (1 - S(:));                % M×1
        W2 = S(:);                      % M×1
        Qnum = W1.*q1 + W2.*q2;         % M×4
        Qnum = Qnum ./ vecnorm(Qnum,2,2);   % 行ごとに正規化
        Qtraj_q = quaternion(Qnum);     % 数値→quaternion配列
    else
        % 通常：slerpで補間
        Qtraj_q = slerp(Qstart, QE, S);
    end

    Qtraj_q = Qtraj_q(:);  % 列ベクトル化（M×1）

    for k = 2:numel(Qtraj_q) % k=2以降について四元数の符号連続性を確保（antipodal問題の回避）
        if sum( compact(Qtraj_q(k)) .* compact(Qtraj_q(k-1)) ) < 0 % 直前との内積が負なら符号が反転していると判断
            Qtraj_q(k) = quaternion(-compact(Qtraj_q(k))); % 同じ回転を表す -q に反転して連続な系列に揃える
        end
    end
    Qtraj_num_global = compact(Qtraj_q); % 四元数オブジェクト配列→数値[N×4]に変換（以降の計算用）

    %% 5) 各個体の評価（parfor）
    parfor i=1:n % 各個体を並列に評価（parfor）
        % -- (5.1) 変数の展開と下限チェック（L, 初期姿勢/初期点）
        x = pop(i,:); % i番目個体の設計変数（1行）を取得
        L = [x(1:5),0.1]; % リンク長ベクトル L=[L1..L5, L6(固定0.1m)]

        if any(L(1:5) < 0.1) % L1..L5 の下限0.1mを違反していないか確認
            J(i,:) = invalid_row; % 違反なら当該個体を無効評価として棄却
            continue; % 次の個体へ
        end

        theta0 = deg2rad([-45 -90 90 0 90 0]); % IKの初期関節角[deg]→[rad]（収束性向上の初期値）
        x0 = [-(L(5) + L(6) - L(2)) * cos(pi / 4); % 初期手先位置x：x=y方向(45°)に配置する近似式
              -(L(5) + L(6) - L(2)) * cos(pi / 4); % 初期手先位置y：xと同値（平面対角方向）
              L(1) + L(3) + L(4)]; % 初期手先位置z：縦方向の長さ和で近似

        % -- (5.2) 参照軌道生成（中点ベジェ）
        Pmid = 0.5*(x0 + Xe); % 始点x0と終点Xeの中点を計算（ベジェの制御点の基準）
        Pmid(3) = Pmid(3) + 0.4; % Z方向に+0.4[m]持ち上げてアーチ状に（地面/干渉の回避余裕）
        Ptraj = (1 - S).^2 .* x0 + 2*(1 - S).*S .* Pmid ... % 2次ベジェ補間で位置軌道を生成
            + S.^2 .* Xe; % Sは最小ジャーク時刻関数（0→1）
        M = size(Ptraj, 2); % 軌道点数（時系列のサンプル数）
        if M ~= size(Qtraj_num_global,1) % 位置と姿勢のサンプル数が一致するか検査
            J(i,:) = invalid_row; % 不一致なら当該個体は無効評価
            continue; % 次の個体へ
        end
        Qtraj_num = Qtraj_num_global; % 姿勢参照（数値四元数）をローカルに取得

        % -- (5.3) IK（全時刻）と角度の位相連続化
        [theta, info] = IK.IK_solver(L, Ptraj, Qtraj_q(:), theta0); % 位置軌道Ptrajと姿勢Qtraj_qに対して全時刻の逆運動学を解く（初期角はtheta0）

        theta = unwrap(theta, [], 2); % 時系列方向（2次元目）で角度の位相を連続化（±2πの飛びを除去）
        
        % -- (5.4) 角速度・角加速度の数値微分
        dth = diff(theta, 1, 2) / dt; % 関節角の前進差分で角速度[rad/s]を算出（時系列方向に1階差分）
        dth = [dth, dth(:,end)]; % 末尾の列を複製して列数をthetaと一致させる（長さ合わせ）
        ddth = diff(dth, 1, 2) / dt; % 角速度の前進差分で角加速度[rad/s^2]を算出
        ddth = [ddth, ddth(:,end)]; % 同様に末尾を複製して列数を一致させる

        pos_acc = 0; ori_acc = 0; % 位置誤差・姿勢誤差（2乗）の時間積分用アキュムレータを初期化

        failed = false; % この個体の評価失敗フラグ（途中で異常が出たらtrueにする）

        % -- (5.5) FKで誤差積分（位置・姿勢）
        for k = 1:M % すべての時刻ステップについて順に評価
            [p_k, q_k] = Kin.FK(L, theta(:,k)); % 正運動学：手先位置p_kと姿勢四元数q_kを算出

            q_k = q_k(:).'; % 四元数を1×4の行ベクトルに整形

            p_ref = Ptraj(:,k); % 参照位置（時刻k）の取り出し
            e_pos = p_k - p_ref;
            pos_acc = pos_acc + sum(e_pos.^2); % 位置誤差の二乗を加算（時間積分用）
            q_ref = Qtraj_num(k,:); % 参照姿勢（四元数, 1×4）の取り出し
            if dot(q_ref, q_k) < 0, q_k = -q_k; end % antipodal対策：内積<0なら符号反転で最近傍回転に
            phi = norm(Quat.quat_log_err(q_ref, q_k)); % 姿勢誤差角[rad]（四元数からの回転角）
            ori_acc = ori_acc + phi^2; % 姿勢誤差角の二乗を加算（時間積分用）
        end
        if failed % 途中で異常があった場合
            J(i,:) = invalid_row; % 無効評価をセット
            continue; % 次の個体へ
        end

        J_pos = pos_acc * dt; % 位置誤差二乗の時間積分
        J_ori = ori_acc * dt; % 姿勢誤差角二乗の時間積分

        torque = Dyn.computeTorque(L, theta, ddth, m_joint); % トルク時系列を推定（動力学モデルに基づく）
   
        tau_peak = max(abs(torque), [], 2); % 関節ごとのピークトルク（時刻方向max）
        tau_peak = max(tau_peak, 1e-9); % 0割防止のための下限クリップ
        safety_margin = tau_max ./ tau_peak; % 安全率＝許容トルク/実負荷トルク

        safety_margin = min(safety_margin, MAX_SAFETY); % 異常に大きい値を上限でクリップ（数値安定化）

        % -- (5.6) トルク・安全率ペナルティ
        J_tau = sum(max(0, SAFETY_THRESHOLD - safety_margin).^2); % 安全率が閾値未満の分だけ二乗罰則を加算
        d_min = Collision.minLinkDistance(L,theta); % 全時刻のリンク間最小距離を評価（自己干渉の指標）

        J_sc = max(0, COLLISION_THRESH - d_min)^2; % 最小距離がしきい値未満なら不足分の二乗で罰則

        % -- (5.7) 自己干渉ペナルティと出力パック
        J(i,:) = [J_tau, J_pos, J_ori, J_sc]; % 4目的を所定の列順で格納（[tau pos ori sc]）
    end
end