# panCB：泛基因组核心 block（core_block）鉴定流程（v6.0）

panCB 以参考基因组为锚，按染色体对多基因组进行成对比对与初筛、对齐优先级排序、核心块迭代识别，最终合并到个体全基因组坐标并提取序列，同时生成统计报告与交互式可视化结果。

## 快速开始
- 参考与输入
  - `-r/--ref`: 参考基因组名称（不含后缀）
  - `-g/--genome_dir`: 包含所有样本 `{name}.fasta` 的目录
  - `-accs`: 样本列表文件，每行一个样本名称（需与 `genome_dir` 中的 fasta 名称一致）
  - `-o/--output`: 结果输出根目录

- 运行示例
```bash
# 标准模式（每轮BLAST验证，准确性最高）
python panCB.py \
  -r REF_NAME \
  -g /path/to/genomes \
  -accs /path/to/accs.txt \
  -o /path/to/output \
  --processes-nucmer 10 --processes-core 5 --resume --step all   # 或 --step core

# 快速迭代模式（推荐，固定每2轮验证+增量验证+动态flanking）
python panCB.py \
  -r REF_NAME \
  -g /path/to/genomes \
  -accs /path/to/accs.txt \
  -o /path/to/output \
  --processes-nucmer 10 --processes-core 5 --fast --resume
```

- 日志：运行日志写入 `output/logs/run_panCB_YYYYmmdd_HHMMSS.log`，并在控制台同步输出。默认记录关键信息（INFO 级别及以上），包括步骤、错误、统计等。如需详细的调试日志，可使用 `--debug` 参数。
- 临时文件清理：流程完成后会自动删除所有染色体的临时目录（`output/genomeFilter_ChrXX/tmp`），以节省磁盘空间。清理过程会记录到日志中，清理失败不会影响流程结果。

---

## 依赖环境

- 外部工具（需在 PATH 中可用）：
  - MUMmer（`nucmer`, `delta-filter`, `show-coords`）
  - samtools（用于 `faidx`）
  - bedtools（`getfasta`, `merge`, `sortBed` 等）
  - TRF（Tandem Repeats Finder，命令 `trf`）
  - BLAST+（`makeblastdb`, `blastn`）
- Python 依赖：
  - Biopython（`Bio`）
  - matplotlib（用于生成趋势图）

---

## 命令行参数

- 必需参数：
  - `-r, --ref`：参考基因组名称（不含后缀）
  - `-g, --genome_dir`：基因组目录，包含 `{name}.fasta`
  - `-o, --output`：输出目录
  - `-accs`：样本名称列表文件

- 可选参数（默认值）：
  - `--processes-nucmer`（5）：nucmer 比对阶段的并行进程数，控制同时运行的比对任务数量
  - `--processes-core`（5）：核心块迭代阶段的并行进程数，控制同时处理的染色体数量
  - `--resume`：断点续跑，若目标产物已存在且非空则跳过
  - `--step`（all）：执行步骤；`all`=同时运行 core 与 alternative；`core`=仅运行 core
  - `--blastn-threads`（1）：核心迭代阶段 `blastn` / `blastn-short` 使用的线程数
    - **性能优化建议**：对于8核CPU，推荐 `--processes-core 4 --blastn-threads 2`（总8线程）；对于16核CPU，推荐 `--processes-core 4 --blastn-threads 4`（总16线程）
    - **重要**：核心迭代阶段总CPU使用 = `processes_core × blastn_threads`，建议不超过CPU物理核心数
  - `--fast`：启用快速迭代模式（v6.0优化版）
    - **默认行为**：不添加此参数时，使用标准迭代模式（每轮执行完整 BLAST 验证）
    - **v6.0核心优化**：
      - **固定验证间隔**：每2轮执行一次BLAST验证（R2, R4, R6...），避免误差累积
      - **增量验证**：每次验证只验证新增样本，而非全量验证所有样本
      - **动态flanking**：根据block长度自适应计算flanking值（block长度的0.5%）
    - **工作原理**：
      ```
      R1: 样本1 ∩ 样本2 → 数学计算更新坐标
      R2: R1 ∩ 样本3   → 增量BLAST验证（只验证样本3）+ 动态flanking
      R3: R2 ∩ 样本4   → 数学计算更新坐标
      R4: R3 ∩ 样本5   → 增量BLAST验证（只验证样本4-5）+ 动态flanking
      ...
      ```
      
## 硬编码参数说明

为简化使用，以下参数已在代码中硬编码，无需通过命令行指定：

### NUCmer 比对参数（genome_alignment.py）
- **min_match = 50**：nucmer 最小匹配长度 `-l`
- **min_cluster = 100**：nucmer 聚簇阈值 `-c`
- **max_gap = 30**：nucmer 允许缺口 `-g`
- **nucmer_threads = 20**：nucmer 线程数 `-t`
- **identity = 90**：`delta-filter -i` 身份度阈值（%），第一步粗筛
- **min_len = 50**：核心块最小长度阈值（bp），用于初筛与迭代

### 核心块识别参数（core_blocks.py）
- **min_len = 50**：核心块最小长度阈值（bp）
- **blastn_evalue = "1e-20"**：长序列 `blastn` 的 e-value 阈值
- **blastn_short_evalue = "1"**：`blastn-short` 的 e-value 阈值
- **blastn_short_word_size = 7**：`blastn-short` 的 `word_size`
- **bedtools_merge_distance = 100**：`bedtools merge -d` 的最大间距
- **trf_params = "2 7 7 80 10 50 500"**：TRF（Tandem Repeats Finder）参数

这些参数值是经过优化的推荐配置，适用于大多数泛基因组分析场景。如需修改这些参数，请直接编辑相应的源代码文件。

### 两步身份度过滤说明

panCB 采用两步身份度过滤策略：

1. **第一步（delta-filter）**：身份度 ≥ 90%（硬编码，不可调整）
   - 位置：`genome_alignment.py` 中的 `delta-filter -i 90`
   - 目的：粗筛，快速过滤掉低质量比对
   - 生成文件：`*.filtered.delta` → `*.filtered.coords`

2. **第二步（core_block_filter）**：身份度 ≥ 98%（可通过 `--core-identity` 调整）
   - 位置：`genome_alignment.py` 中的 `core_block_filter` 函数
   - 目的：精筛，确保核心block的高质量
   - 生成文件：`*.filtered.1227`

**调整建议**：
- 如果核心区域太少，可降低 `--core-identity` 到 95-97
- 如果核心区域质量不够高，可提高 `--core-identity` 到 99

---

## 流程说明（panCB pipeline）

下面以 `panCB.py` 驱动的 5 个步骤为主线，说明关键产物与目录。

### 步骤 0：索引准备（自动）
- 脚本会为 `accs` 中所有 `{name}.fasta` 预生成 `samtools faidx` 索引（`.fai`）。
- 日志命名空间：`panCB.index`。

### 步骤 1：按染色体拆分（`Chr_split/`）
- 功能：从每个 `{name}.fasta` 中按参考基因组染色体 ID 拆分，生成 per-chr fasta 切片。
- 入口：`genome_alignment.split_all_genomes_by_chr`
- 产物：
  - `output/Chr_split/ChrXX/{sample}_ChrXX.fasta`
  - 若样本缺失该染色体，则写入 1 bp 的占位序列 `N`（可断点续跑检测非空）。
- 日志命名空间：`panCB.split`。

### 步骤 2：成对比对与初筛（染色体级任务）
- 入口：`genome_alignment.run_all_alignments_for_chr`
- 对参考与其它样本在同一染色体上进行 NUCmer 比对→`delta-filter`→`show-coords`，随后进行“核心块初筛”。
- 关键参数（已硬编码）：min_match=50, min_cluster=100, max_gap=30, nucmer_threads=20, identity=90, min_len=50
- 产物（以染色体 `ChrXX` 为例）：
  - `output/genomeFilter_ChrXX/`（每染色体工作根目录）
    - `REF_RP_QRY.*`：NUCmer 产出（`.delta`、`.filtered.delta`、`*.filtered.coords` 等）
    - `coreB/`：核心结果目录
      - `REF_RP_QRY.filtered.coords.{min-len}.filtered.1227`：核心候选（初筛）
- 初筛规则见 `genome_alignment.core_block_filter`（覆盖率、身份度、最小长度、坐标去重/合并等）。
- 并发（两阶段并行策略）：
  - 阶段1（比对）：按 `--processes-nucmer` 控制 nucmer 比对任务的并行数量，所有样本×染色体的比对任务作为独立任务并行处理。
  - 阶段2（迭代）：按 `--processes-core` 控制染色体级别的并行数量，同时处理多条染色体的核心块迭代识别。
- 日志命名空间：`panCB.align`。

### 步骤 3：对齐优先级排序（不生成 alignSum 目录）
- 入口：`core_blocks.generate_alignment_priority_for_chr`
- 逻辑：统计每个样本在该染色体的核心长度总和，据此在内存排序；不生成 `alignSum` 目录与 `.acc/.acc-nh` 文件。
- 顺序记录：写入 `output/genomeFilter_ChrXX/coreB/REF.ChrXX.iter`（仅样本名，一行一个，参考基因组不在该文件中）。
- 并行执行：阶段2中，单条染色体的排序会在其比对结束后立即进行。
- 日志命名空间：`panCB.stats`。

### 步骤 4：核心块迭代识别（`coreB/`）
- 入口：`core_blocks.run_iterative_core_blocks_for_chr`
- 初始化：取排序前两名样本 H1、H2 作为 R1 的 `(ref, qry)`；之后依次将剩余样本加入，迭代生成 `Rk`。
- 单次迭代核心识别：`core_blocks.extract_core_blocks_once`
  - 依赖：`samtools`, `trf`, `blastn`, `sortBed`, `bedtools`
  - 中间日志/检查：
    - `output/genomeFilter_ChrXX/tmp/ChrXX.warning`（累计追加）
    - `output/genomeFilter_ChrXX/tmp/ChrXX.tmp.*`（临时文件前缀）
    - `coreB/ChrXX.{min-len}.{outf}-trans` 与 `coreB/ChrXX.{min-len}.{outf}-inters`（转座/重叠处理记录）
  - 主要产物：
    - `coreB/ChrXX.{min-len}.Rk_RP_qryRecord.coords`（各轮结果，`k=1..`）
  - 迭代趋势可视化：
    - `coreB/core_block_trend/ChrXX_trend.tsv`（记录每轮迭代的block数量）
    - `coreB/core_block_trend/ChrXX_trend.png`（block数量变化趋势折线图）
- 结束：记录最后一轮 `R*` 的路径于日志。
- 并行执行：受主控进程池调度，同一时间可对多条染色体并行迭代（数量上限为 `--processes-core`）。
- 日志命名空间：`panCB.core`。
- **重要说明**：迭代过程中block数量通常呈上升趋势，这是正常现象。详见 `迭代过程详解.md`。

### 步骤 5：生成染色体级核心坐标（`ChrXX.Os.*`）
- 入口：`sequence_process.extract_chr_level_coords`
- 触发：每条染色体的迭代完成后立即执行，可单独断点续跑。
- 作用：把最终 `R*` 结果拆分为参考 + 所有样本的 BED（三列：`chr  start  end`）。
- 产物：`output/genomeFilter_ChrXX/coreB/ChrXX.Os.RP_qryRecord.{sample}`
- 日志命名空间：`panCB.seq`。

### 步骤 5.5：生成染色体级别核心block变化趋势（`coreB/core_block_trend/`）
- 触发：步骤4迭代过程中自动记录
- 产物：
  - `output/genomeFilter_ChrXX/coreB/core_block_trend/ChrXX_trend.tsv`（记录每轮迭代的block数量）
  - `output/genomeFilter_ChrXX/coreB/core_block_trend/ChrXX_trend.png`（block数量变化趋势折线图）
- 日志命名空间：`panCB.core`。

### 步骤 6：合并坐标与序列提取（`pan_core_block/`）
- 入口：`sequence_process.merge_chr_coords_and_extract_sequences`
- 输入：所有 `ChrXX.Os.*` 染色体级坐标
- 跨染色体合并为全基因组 BED：
  - `output/genomeFilter_ChrXX/coreB/ChrXX.Os.RP_qryRecord.{sample}`
  - `output/pan_core_block/RP_qryRecord.{sample}`
- 提取序列（按样本）：
  - `output/pan_core_block/fasta/RP_qryRecord.{sample}.fasta`（来自 `bedtools getfasta`）
- 日志命名空间：`panCB.seq`。
- 断点续跑：检查 `pan_core_block/fasta/` 目录下是否存在非空的 fasta 文件。

### 步骤 6.5：可选 - alternative_block（按染色体输出）
- 触发方式：在命令行指定 `--step all`（默认）。若设为 `--step core`，则跳过本步骤。
- 入口：`alternative_block.run_alternative_for_all_chr`
- 产物：`output/genomeFilter_ChrXX/alterB/`
- 断点续跑：检查 `alterB/` 目录是否存在。

### 步骤 7：核心长度统计（`pan_core_block/report/`）
- 入口：`report.summarize_core_lengths`
- 产物：
  - `output/pan_core_block/report/core_len_per_chr.tsv`（列：`chr  sample  core_len`）
  - `output/pan_core_block/report/core_len_total.tsv`（列：`sample  total_core_len`）
- 日志命名空间：`panCB.report`。
- 断点续跑：检查 TSV 文件是否存在且非空。

### 步骤 8：详细文本报告和可视化报告
- 入口：`report.write_detailed_summary` + `report.generate_all_visualizations`
- 产物：
  - `output/pan_core_block/report/core_summary.txt`（逐样本×逐染色体的块数量、总长度、平均长度）
  - `output/pan_core_block/report/core_summary.html`（交互式HTML统计报告）
- 功能：
  - TXT报告：详细的文本统计信息
  - HTML报告：柱状图 + 可点击染色体矩形 + 详情表格
- 日志命名空间：`panCB.report`。
- 断点续跑：检查 TXT 和 HTML 文件是否存在且非空。

### 步骤 9：基因组级别趋势报告
- 入口：`core_blocks.generate_genome_trend_report`
- 触发：流程末尾自动执行，读取所有染色体的 TSV 趋势文件
- 产物：
  - `output/pan_core_block/report/core_block_trend_report.html`（交互式HTML报告）
- 功能：
  - 上部柱状图：展示各染色体最终block数量
  - 中部指示器：可点击选择染色体
  - 下部内容：点击染色体后显示迭代趋势折线图（上）和详细数据表格（下）
- 日志命名空间：`panCB.core`。
- 断点续跑：检查 HTML 文件是否存在且非空。


    - `ChrXX.RP_qryRecord.{sample}.alternative.bed`
    - `ChrXX.RP_qryRecord.{sample}.sorted.len`
  - `output/genomeFilter_ChrXX/alterB/alter-FASTA/`
    - `ChrXX.RP_qryRecord.{sample}.sorted.fasta`
    - `ChrXX.RP_qryRecord.{sample}.alternative.fasta`

注意：panCB 不会生成原版中的 `genomeFilter/coreB/alter-PanRice/alternative-block/` 目录，也不会进行跨染色体的 alternative 合并，所有结果均为染色体级。

---

## 输出目录结构总览

以下展示典型输出（省略若干中间文件名）：

```text
output/
├─ logs/
│  └─ run_panCB_YYYYmmdd_HHMMSS.log
├─ checkpoints/
│  └─ resume_state.json                    # 断点状态文件
├─ Chr_split/
│  └─ Chr01/
│     ├─ REF_Chr01.fasta
│     ├─ SAMPLE1_Chr01.fasta
│     └─ ...
├─ genomeFilter_Chr01/
│  ├─ REF_RP_SAMPLE1.delta
│  ├─ REF_RP_SAMPLE1.filtered.delta
│  ├─ REF_RP_SAMPLE1.filtered.coords
│  ├─ coreB/
│  │  ├─ REF_RP_SAMPLE1.filtered.coords.{min-len}.filtered.1227
│  │  ├─ Chr01.{min-len}.R1_RP_qryRecord.coords      # 迭代结果文件（Level 1）
│  │  ├─ Chr01.{min-len}.R2_RP_qryRecord.coords
│  │  ├─ Chr01.{min-len}.R*_RP_qryRecord.coords      # 最终迭代结果
│  │  ├─ Chr01.{min-len}.R*_RP_qryRecord.coords-trans
│  │  ├─ Chr01.{min-len}.R*_RP_qryRecord.coords-inters
│  │  ├─ Chr01.Os.RP_qryRecord.REF                   # 染色体坐标文件（Level 3）
│  │  ├─ Chr01.Os.RP_qryRecord.SAMPLE1
│  │  ├─ REF.Chr01.iter                              # 迭代顺序文件
│  │  └─ core_block_trend/                           # 趋势文件目录（Level 2）
│  │     ├─ Chr01_trend.tsv                          # 趋势TSV数据
│  │     └─ Chr01_trend.png                          # 趋势折线图
│  └─ tmp/  （注：流程完成后会自动删除）
│     ├─ Chr01.warning
│     └─ Chr01.tmp.*
│  └─ alterB/
│     ├─ alter-BED/
│     │  ├─ Chr01.RP_qryRecord.SAMPLE1.sorted.bed
│     │  ├─ Chr01.RP_qryRecord.SAMPLE1.sorted.bed.named
│     │  ├─ Chr01.RP_qryRecord.SAMPLE1.alternative.bed
│     │  └─ Chr01.RP_qryRecord.SAMPLE1.sorted.len
│     └─ alter-FASTA/
│        ├─ Chr01.RP_qryRecord.SAMPLE1.sorted.fasta
│        └─ Chr01.RP_qryRecord.SAMPLE1.alternative.fasta
├─ genomeFilter_Chr02/
│  └─ ...（同上结构）
└─ pan_core_block/
   ├─ RP_qryRecord.REF                               # 全基因组坐标文件
   ├─ RP_qryRecord.SAMPLE1
   ├─ fasta/
   │  ├─ RP_qryRecord.REF.fasta                      # 核心序列文件
   │  └─ RP_qryRecord.SAMPLE1.fasta
   └─ report/
      ├─ core_len_per_chr.tsv                        # 染色体级别长度统计
      ├─ core_len_total.tsv                          # 全基因组长度统计
      ├─ core_summary.txt                            # 详细文本报告
      ├─ core_summary.html                           # 交互式HTML统计报告
      └─ core_block_trend_report.html                # 迭代趋势HTML报告（Level 4）
```

### 文件含义说明
- `REF_RP_QRY.filtered.coords`：`show-coords` 输出的简化坐标表（过滤后），用于初筛的输入。
- `*.filtered.coords.{min-len}.filtered.1227`：初筛后的核心候选区块（每条记录包含参考与查询坐标）。
- `ChrXX.{min-len}.Rk_RP_qryRecord.coords`：第 k 轮迭代生成的核心区块联合坐标。
- `ChrXX.{min-len}.Rk_RP_qryRecord.coords-trans`：迭代中“转座/错序”剔除记录。
- `ChrXX.{min-len}.Rk_RP_qryRecord.coords-inters`：迭代中“重叠切分/删除”记录。
- `ChrXX.Os.RP_qryRecord.{sample}`：按样本在该染色体上的 BED 坐标（3 列：`chr  start  end`）。
- `pan_core_block/RP_qryRecord.{sample}`：跨染色体合并后的全基因组 BED 坐标。
- `pan_core_block/fasta/RP_qryRecord.{sample}.fasta`：由 BED 抽取的核心序列。
- `report/core_len_per_chr.tsv`：每样本×每染色体的核心长度统计。
- `report/core_len_total.tsv`：每样本全基因组核心长度总和。
- `report/core_summary.txt`：文本化的详细统计报告。
- `report/core_summary.html`：交互式HTML统计报告，支持按染色体查看详细数据。
- `report/core_block_trend_report.html`：基因组级别的迭代趋势交互式HTML报告，点击染色体可查看折线图和数据表格。
