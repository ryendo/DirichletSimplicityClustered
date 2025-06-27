% ディレクトリの設定
dataDir = 'results'; % CSVファイルが保存されているディレクトリ
outputFile = 'bounds_Jun14_1059.csv'; % 出力ファイル名

% CSVファイルのリストを取得
csvFiles = dir(fullfile(dataDir, 'bounds_*.csv'));

% 初期化
combinedData = [];

% 各CSVファイルを読み込み、結合
for k = 1:length(csvFiles)
    % ファイル名の取得
    csvFileName = fullfile(dataDir, csvFiles(k).name);
    
    % データの読み込み
    data = readmatrix(csvFileName);
    
    % データの結合
    combinedData = [combinedData; data]; %#ok<AGROW>
end

% 結合データの保存
writematrix(combinedData, fullfile(dataDir, outputFile));

disp('CSVファイルの統合が完了しました。');
