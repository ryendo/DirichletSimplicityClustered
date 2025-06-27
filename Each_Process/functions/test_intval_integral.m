% % INTLABを初期化
% addpath('path_to_intlab');  % INTLABのパスを追加
% startintlab;  % INTLABを起動

% 区間を定義
a = infsup(1, 2);  % 区間 [1, 2]

% 区間関数を定義
f = @(x) x.^2;

% 積分を実行
result = integral(f, a);
disp(['区間積分の結果: ', result])
