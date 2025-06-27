% MATLABコード：1〜1220の番号を1行に並べたリストをCSVファイルに保存
numbers = 1:1220; % 1から1220までの数値を生成
csvwrite('list_j.csv', numbers'); % 'list_j.csv'というファイル名で保存
