function [l1, l2, l3, l4] = get_approximate_eigenvalue(query_a, query_b)
% GET_APPROXIMATE_EIGENVALUE 任意の(a,b)における固有値を線形補間で求める
%
%   [l1, l2, l3, l4] = get_approximate_eigenvalue(a, b)
%
%   入力:
%       query_a, query_b : 求めたい点の座標（スカラー、または同じサイズの配列）
%   出力:
%       l1, l2, l3, l4   : 線形補間された固有値
%
%   備考:
%       初回実行時に 'eigenvalues_grid.csv' を読み込み、
%       scatteredInterpolant オブジェクトを作成してメモリに保持します。

    % 永続変数の定義（関数が終了してもメモリに残る変数）
    persistent F1 F2 F3 F4

    % 初回呼び出し時（F1が空のとき）のみデータを読み込む
    if isempty(F1)
        filename = 'eigenvalues_grid.csv';
        
        if ~isfile(filename)
            filename = 'EigenvaluesInterpolation/eigenvalues_grid.csv';
            if ~isfile(filename)
                filename = 'EigenvaluesInterpolation/eigenvalues_grid.csv';
                error('ファイル %s が見つかりません。パスを確認してください。', filename);
            end
        end
        
        % CSV読み込み
        opts = detectImportOptions(filename);
        opts.VariableNamingRule = 'preserve'; % 列名を維持
        T = readtable(filename, opts);
        
        % データ抽出
        a_data = T.a;
        b_data = T.b;
        
        % 補間オブジェクトの作成 (線形補間: 'linear')
        % scatteredInterpolant はデータの並び順を気にせず使えて便利です
        F1 = scatteredInterpolant(a_data, b_data, T.lambda1, 'linear', 'linear');
        F2 = scatteredInterpolant(a_data, b_data, T.lambda2, 'linear', 'linear');
        F3 = scatteredInterpolant(a_data, b_data, T.lambda3, 'linear', 'linear');
        F4 = scatteredInterpolant(a_data, b_data, T.lambda4, 'linear', 'linear');
        
        fprintf('Initialization: CSV loaded and interpolants created.\n');
    end

    % 補間値の計算
    l1 = F1(query_a, query_b);
    l2 = F2(query_a, query_b);
    l3 = F3(query_a, query_b);
    l4 = F4(query_a, query_b);

end