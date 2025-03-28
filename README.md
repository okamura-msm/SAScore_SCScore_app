# 合成容易性スコア可視化アプリ（SAScore / SCScore）

このStreamlitアプリは、有機分子に対して **SAScore** および **SCScore** を計算・比較するためのツールです。分子構造の複雑さや合成経路の難しさを評価し、直感的に可視化できます。

## 主な機能

- **SAScore**（構造ベースの合成容易性スコア）
- **SCScore**（逆合成ステップ数に基づく合成複雑性スコア）
- SCScoreは1〜10にスケーリングしてSAScoreと比較しやすくした
- SMILESまたはCSVファイルでの入力に対応
- スコア付きの表と分子画像の表示
- 結果をCSV形式でダウンロード可能
- スコア解説セクションもアプリ内に組み込んだ

## フォルダ構成

```
synthetic-accessibility-app/
├── app/
│   ├── SAscore_app.py          # Streamlitアプリ本体
│   ├── sascorer.py             # SAScore計算用モジュール
│   ├── scscore_numpy.py        # SCScore（NumPy実装）
│   ├── fpscores.pkl.gz         # SAScore用の部分構造スコア辞書
│   └── model_files/            # （オプション）SCScoreの学習済モデル
├── sample_data/
│   └── demo_molecules.csv      # デモ用SMILES一覧（準備中）
├── requirements.txt            # 必要なPythonパッケージ
├── LICENSE                     # ライセンス（MIT）
└── README.md                   # このファイル
```

## インストール手順

### 1. リポジトリをクローン
```bash
git clone https://github.com/yourusername/synthetic-accessibility-app.git
cd synthetic-accessibility-app
```

### 2. 必要なライブラリをインストール
```bash
pip install -r requirements.txt
```

`rdkit` が入っていない場合は以下でインストール：
```bash
pip install rdkit-pypi
```

### 3. Streamlitアプリを実行
```bash
streamlit run app/SAscore_app.py
```

---

## SAScoreについて

構造の複雑さや珍しさ、部分構造の出現頻度などからスコアを付ける指標です。  
`sascorer.py` にて計算を行い、必要な辞書ファイル `fpscores.pkl.gz` は `app/` フォルダに同梱されています。

元論文：
- Peter Ertl and Greg Landrum (2010), [J. Med. Chem.]

---

## SCScoreについて

Reaxys反応データに基づき、ニューラルネットが「必要な合成ステップ数」を予測するスコアです。ここでは軽量なNumPy実装を使用しています（TensorFlow不要）。

元論文：
- Coley et al., 2018（MIT） [https://github.com/connorcoley/scscore]

---

## ライセンス

このプロジェクトは MITライセンス のもとで公開されています。

---

## 作成者

作成者：[あなたの名前またはGitHubアカウント]

フィードバックや改善案も歓迎です！

