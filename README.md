# CAPS-Maker

2つのDNA配列（野生型・変異型など）からCAPSマーカーを探索する、日本語GUI付きのオープンソースツールです。[primerblast-oss](https://github.com/light-suzuki/primerblast-oss)の制限酵素切断シミュレーションとPrimer3/BLAST+処理を利用します。

**注意：** 本ソフトはin silico予測です。実際のPCR・制限酵素消化・ゲル検証を置き換えるものではありません。

## 機能

- 2つのFASTAファイルをアップロード、または配列を直接貼り付け
- グローバルアラインメントでSNP・InDelを検出
- 全制限酵素の探索、または指定酵素のみ（複数選択）
- A/Bのどちらが切断されるか、AA・AB・BBの断片長と切断部位数を一覧化
- 予測ゲルパターンから共優性・優性型・判別困難を判定（遺伝的優性とは異なる）
- 最小断片長、ゲル評価、推奨ゲル濃度、仮想ゲル（SVG）表示/保存
- Primer-BLASTモードではPCR産物長の下限・上限を指定、Primer3とBLAST+による特異性検証
- CSV（UTF-8 BOM付き、Excel対応）/TSV出力、プライマー発注表

## 必要条件

- Python 3.10以降、Git
- Primer-BLASTモードのみ、primer3_core・blastn、および作成済みのローカルBLAST nucleotide DB
- Streamlit、Biopython、pandas、primerblast-oss（requirements.txtに記載）

### Ubuntu / WSL2で起動

~~~bash
sudo apt update
sudo apt install -y python3-venv git primer3 ncbi-blast+

git clone https://github.com/light-suzuki/CAPS-Maker.git
cd CAPS-Maker
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt

streamlit run app.py
~~~

GUIが起動したらブラウザで http://localhost:8501 を開きます。Windows側のブラウザからWSL上のStreamlitに接続することもできます。

### Primer-BLASTモードの準備

ゲノムFASTAから、あらかじめBLAST DBを作成します。

~~~bash
makeblastdb -in reference_genome.fa -dbtype nucl -parse_seqids -out reference_db
~~~

GUIの「BLAST DBのローカルパス」に reference_db を指定してください。野生型/変異型の2配列とは別に、ゲノム全体での特異性確認用のDBが必要です。複数のDBを使うときは1行1パスで入力します。

## 使い方

1. 左右に配列A/Bを入力します（各1本、現行GUIでは8,000 bp以内）。
2. モードを選びます。
   - **制限酵素スキャン**：入力DNA配列の全長を仮想消化します。PCR産物の切断予測ではなく、候補酵素の探索に用います。
   - **Primer-BLAST付きマーカー設計**：検出したSNP・InDelを選び、PCR産物長の下限/上限、BLAST DB、Primer3候補数を指定します。両アレルで塩基配列が完全一致する共通プライマーだけを残して、予測PCR産物を切断解析します。
3. 対象酵素を「利用可能な全酵素」または「指定酵素」から選びます。
4. 「CAPSマーカーを探索」をクリックします。
5. 候補表をフィルタし、CSV/TSVおよびゲル画像（SVG）をダウンロードできます。

### 出力例

| 酵素 | 切断アレル | A断片(bp) | B断片(bp) | AB断片(bp) | マーカー型 |
|---|---|---|---|---|---|
| EcoRI | Aのみ | 300 + 200 | 500 | 500 + 300 + 200 | 共優性（3遺伝子型を識別） |

数値は説明用です。実際の出力には切断部位数、バンド分離、ゲル評価、Primer-BLASTの特異性情報も含まれます。

## マーカー型の意味

- 共優性：AA・AB・BBの予測バンドパターンを3つとも区別できる。
- 優性型：ホモ接合体同士は区別できても、ヘテロ接合体が一方のホモ接合体と同じに見える。
- 判別困難：ホモ接合体同士のゲルパターンも明確に区別できないなど、診断に不適。

これらは検出パターンの区別であり、変異の生物学的な優性/劣性とは異なります。分離の程度はゲルモデルで推定しており、実験上の分離を保証しません。

## 科学的な注意事項

- 完全消化と両アレルの増幅を仮定します。小断片は見えないことがあり、ゲル濃度、DNA量、酵素条件で変わります。
- 初版では両アレルに完全一致するプライマーのみを採用します。ミスマッチを許容するPCRは扱いません。
- BLASTがspecificと評価しても、検索が不完全なら特異性は保証されません。Bアレルの局所配列にプライマーが結合することと、Bゲノム全体の検索は別です。
- 仮想ゲルは標的PCR産物の消化モデルであり、オフターゲット産物由来のバンドは重ねていません。
- dCAPS（人工ミスマッチによる制限酵素認識部位作成）は、現行GUIでは未対応です。

## テスト

~~~bash
python -m pip install pytest
python -m pytest -q
~~~

## 今後の拡張候補

- 実験室で保有する酵素リストの保存と酵素価格・反応条件による順位付け
- 複数SNPの一括解析・複数Primer3設計
- dCAPSフォールバック
- 実際のPCR/消化結果の記録・予測との照合
- プライマー3′末端付近の多型チェックの強化
- 数万塩基以上の座位に対応する高速アラインメント

## ライセンス

MIT License。primerblast-ossもMIT Licenseで公開されています。依存ライブラリ・データのライセンスは各配布元の規定に従ってください。
