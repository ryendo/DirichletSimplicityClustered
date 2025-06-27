% 自動微分変数の作成
x = gradientinit(infsup(2,2.1)); % 初期値 2 を持つ自動微分変数

% 関数の定義
f = x^2 + sin(x);

% 結果の表示
disp('Function value:');
disp(f.x); % 関数値

disp('Derivative value:');
disp(f.dx); % 微分値
