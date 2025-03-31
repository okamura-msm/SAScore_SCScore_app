import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from PIL import Image
from io import BytesIO
from sascorer import calculateScore
from scscore_numpy import SCScorer

st.set_page_config(page_title="合成容易性スコア比較", layout="wide")
st.title("合成容易性スコア比較アプリ（SAScore / SCScore）")

with st.expander("🔎 SAScore / SCScoreとは？（クリックで開く）"):
    st.markdown("""
### 🧪 SAScore（Synthetic Accessibility Score）

- **構造的な見た目・珍しさ・サイズ・対称性**などに基づき、分子の「作りやすさ」を評価
- RDKitベースのルール＋頻度データにより計算
- スコア範囲：**1（合成容易）〜10（困難）**

**長所**
- 計算が非常に高速
- 多量の分子のスクリーニングに最適
- 見た目の複雑さや官能基の珍しさをよく捉える

**短所**
- 実際の合成可能性と乖離する場合あり（市販品でも高スコアなど）
- 合成経路は考慮されていない

---

### 🧪 SCScore（Synthetic Complexity Score）

- 逆合成的な観点から「何ステップで作れそうか？」をニューラルネットで予測
- Reaxys由来の反応データで学習されたモデルを使用
- スコア範囲：**1（簡単）〜10（複雑）** ※アプリ内では1〜10にスケーリング済

**長所**
- 実際の合成ステップ数と比較的よく相関する
- 構造の規則性や簡便さを反映しやすい

**短所**
- 計算はSAScoreよりやや重い
- 一部の構造（商用品、非標準構造）に過小評価・過大評価が出ることも

---

### 🎯 比較まとめ

| 特性        | SAScore      | SCScore      |
|-------------|--------------|--------------|
| 観点        | 構造の見た目・頻度    | 合成経路の複雑さ     |
| スコア範囲     | 1〜10（低いほど容易） | 1〜10（低いほど容易） |
| 計算速度      | ◎ 非常に高速      | ◯ 中程度        |
| 合成現実性への対応 | △（構造偏重）      | ◯（経路を考慮）     |
| 向いている用途   | 初期スクリーニング    | 詳細評価・最終絞り込み  |
""")

@st.cache_resource
def load_scscore():
    model = SCScorer()
    model.restore()
    return model

sc_model = load_scscore()

def rescale_scscore(score, old_min=1.0, old_max=6.0, new_min=1.0, new_max=10.0):
    score = float(score)
    scaled = (score - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
    return max(min(scaled, new_max), new_min)

def calculate_scores(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        sa_score = calculateScore(mol)
        sc_score_raw = sc_model.get_score_from_smi(smiles)[1]
        sc_score = rescale_scscore(sc_score_raw)
        mw = Descriptors.MolWt(mol)
        return sa_score, sc_score, mw, mol
    except:
        return None, None, None, None

def highlight_score(val):
    if isinstance(val, (int, float)):
        if val < 3.0:
            return 'background-color: lightgreen'
        elif val < 6.0:
            return 'background-color: orange'
        else:
            return 'background-color: red'
    return ''

def safe_mol_to_image(mol, legend=None):
    try:
        return Draw.MolToImage(mol, size=(300, 300), legend=legend)
    except Exception as e:
        print(f"[Warning] MolToImage failed: {e}")
        return Image.new('RGB', (300, 300), color=(255, 255, 255))

def mols_to_image_grid(mols, legends, mols_per_row=6):
    images = [safe_mol_to_image(mol, legend=legend) for mol, legend in zip(mols, legends)]
    rows = [images[i:i + mols_per_row] for i in range(0, len(images), mols_per_row)]
    full_img = Image.new('RGB', (300 * mols_per_row, 300 * len(rows)), color=(255, 255, 255))
    for row_idx, row_imgs in enumerate(rows):
        for col_idx, img in enumerate(row_imgs):
            full_img.paste(img, (300 * col_idx, 300 * row_idx))
    return full_img

input_method = st.radio("入力方法を選択:", ("SMILES を直接入力", "CSVファイルをアップロード"))

if input_method == "SMILES を直接入力":
    smiles_input = st.text_area("1行に1つずつ SMILES を入力してください:")
    if st.button("計算する"):
        smiles_list = [s.strip() for s in smiles_input.splitlines() if s.strip()]
        data = []
        mols = []
        for smi in smiles_list:
            sa, sc, mw, mol = calculate_scores(smi)
            if sa is not None:
                data.append((smi, sa, mw, sc))
                mols.append(mol)

        df = pd.DataFrame(data, columns=["SMILES", "SA Score", "Molecular Weight", "SC Score"])
        st.dataframe(df.style
            .applymap(highlight_score, subset=["SA Score"])
            .applymap(highlight_score, subset=["SC Score"]))

        if mols:
            st.subheader("分子構造（スコア付き）")
            legends = [f"SA: {sa:.2f}, SC: {sc:.2f}" for _, sa, _, sc in data]
            img = mols_to_image_grid(mols, legends)
            st.image(img)

elif input_method == "CSVファイルをアップロード":
    uploaded_file = st.file_uploader("CSVファイルを選択（SMILES列を含む必要あり）")
    if uploaded_file is not None:
        df_in = pd.read_csv(uploaded_file)
        if "SMILES" not in df_in.columns:
            st.error("CSVに 'SMILES' カラムが必要です。")
        else:
            df_out = []
            mols = []
            for i, row in df_in.iterrows():
                smi = row["SMILES"]
                sa, sc, mw, mol = calculate_scores(smi)
                if sa is not None:
                    df_out.append((smi, sa, mw, sc))
                    mols.append(mol)

            df = pd.DataFrame(df_out, columns=["SMILES", "SA Score", "Molecular Weight", "SC Score"])
            st.dataframe(df.style
                .applymap(highlight_score, subset=["SA Score"])
                .applymap(highlight_score, subset=["SC Score"]))

            if mols:
                st.subheader("分子構造（スコア付き）")
                legends = [f"SA: {sa:.2f}, SC: {sc:.2f}" for _, sa, _, sc in df_out]
                img = mols_to_image_grid(mols, legends)
                st.image(img)

            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("結果をCSVでダウンロード", csv, "scores.csv", "text/csv")
