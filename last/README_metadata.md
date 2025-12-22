# Metadata-Driven Layout

## 出力ディレクトリの基本構造

`dwl_config.yaml` の `paths.out_dir` 直下に、短縮名の組み合わせと日付でネストしたランディレクトリが作成されます。`trisbst_3spc_fromDwl.sh` が `metadata/` サブディレクトリを事前に生成し、NCBI Datasets から取得したサマリ JSON を配置します。

```
<out_dir>/
  Homsap1_Chlre2_Scere3/
    20251109/
      metadata/
        Homsap1_GCA_000001405.29.json
        Chlre2_GCA_000002595.3.json
        Scere3_GCA_000146045.2.json
        metadata_manifest.json
      Homsap12Chlre2Scere3_many2one_20251109.maf
      ...
```

- `metadata/*.json`: それぞれのアクセッションごとの `datasets summary genome accession` 出力。ファイル名は `<short_name>_<accession>.json`。
- `metadata_manifest.json`: ラン全体のメタデータ。`collect_run_summary.py` はここだけを情報源として org2/org3 のアクセッションや short_name を解決します。

## metadata_manifest.json の主なフィールド

| フィールド | 説明 |
| --- | --- |
| `date` | ラン日付（コマンド引数の DATE）。 |
| `combo_dir` | `<org1_short>_<org2_short>_<org3_short>` 形式のラン識別子。 |
| `run_dir` / `metadata_dir` | 絶対パス。後続解析はこれらを辿れば必要ファイルにアクセスできます。 |
| `organisms[]` | org1/org2/org3 の配列。各要素には `slot`（org1/2/3）, `role`（outgroup/ingroup）, `accession`, `directory_name`, `short_name`, `fasta_path`, `gff_path`, `raw_organism_name`, `ncbi_full_name`, `metadata_json` が含まれます。 |

### collect_run_summary.py との連携

- org2/org3 の short_name を `metadata_manifest.json` から取得し、TSV や PDF のファイル探索に利用します。
- アクセッション ID・生物名・サマリ JSON パスも manifest から直接取得するため、ファイル名の解釈や追加の `datasets` 呼び出しは不要になりました。

## 運用メモ

- `trisbst_3spc.sh` 側でも `dwl_config.yaml` の `paths.out_dir` を参照するため、ランナーと解析で同じ設定値を共有できます。
- 既存のランを再解析する場合は、`metadata/metadata_manifest.json` が存在することを確認してください（欠けている場合は `trisbst_3spc_fromDwl.sh` を再実行し、メタデータのみ再生成してからパイプラインを回せます）。

