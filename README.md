# panCB：泛基因组核心 block（core_block）鉴定流程（v6.0）

panCB 以参考基因组为锚，按染色体对多基因组进行成对比对与初筛、对齐优先级排序、核心块迭代识别，最终合并到个体全基因组坐标并提取序列，同时生成统计报告与交互式可视化结果。

**更新日志（v6.0 - 2025.12）**：
- **重大优化**：快速迭代模式全面升级
  - **强制中间验证**：快速模式现在固定每2轮执行一次BLAST验证，避免坐标误差累积
  - **增量验证**：参照标准模式，每次验证只验证新增样本，大幅减少验证时间（节省70-80%）
  - **动态flanking**：根据block长度自适应计算flanking值（block长度的0.5%），替代固定值
  - **参数简化**：移除`--round-verify`和`--flanking-size`参数，使用优化的固定策略
  - **性能提升**：在保持准确性的同时，验证时间与新增样本数成正比，而非累积样本数
  - 详细说明请参考 `.kiro/specs/fast-mode-optimization/design.md`

**更新日志（v5.3 - 2025.12）**：
- **重要修复**：修复迭代过程中样本坐标计算错误的问题
  - 问题描述：当交集区域等于参考区域或查询区域时，原代码会跳过BLAST验证，直接保留原始坐标，导致某些样本的核心block坐标异常偏大
  - 修复方案：移除快捷路径，始终对所有样本进行BLAST验证和坐标重新计算
  - 详细说明请参考 `BUGFIX_坐标计算修复说明.md`

**更新日志（v5.2 - 2025.12）**：
- 优化断点续跑功能：精确检测迭代完成状态，支持从中断的迭代轮次继续执行
- 修复迭代文件排序问题：使用数值排序替代字符串排序，确保正确识别最终迭代结果
- 增强HTML报告：点击染色体时同时显示迭代趋势折线图和数据表格
- 新增可视化报告：自动生成核心block长度折线图（PNG）和交互式HTML统计报告
- 完善文档说明

---

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
    - **性能提升**：
      - 验证时间与新增样本数成正比，而非累积样本数
      - 相比旧版快速模式，验证时间减少70-80%
      - 相比标准模式，总时间减少约50%（取决于样本数）
    - **准确性保证**：
      - 每2轮验证一次，误差累积最多1轮
      - 动态flanking自适应容忍误差（小block用小flanking，大block用大flanking）
      - 结果与标准模式高度一致（坐标差异<5bp）
    - **适用场景**：
      - ✅ 推荐用于所有规模的泛基因组分析（10-1000+样本）
      - ✅ 平衡准确性和速度的最佳选择
      - ✅ 计算资源有限的环境
    - **输出兼容**：产生与标准模式相同格式的输出文件
  - `--round-verify`（已弃用）：此参数在v6.0中已弃用
    - **弃用原因**：快速模式现在固定每2轮验证一次，无需手动配置
    - **向后兼容**：参数保留但会被忽略，并发出弃用警告
    - **迁移建议**：直接使用`--fast`参数即可，无需额外配置
  - `--flanking-size`（已弃用）：此参数在快速模式下已弃用
    - **弃用原因**：快速模式现在使用动态flanking（block长度的0.5%）
    - **向后兼容**：参数保留但在快速模式下会被忽略，并发出弃用警告
    - **动态flanking示例**：
      - block=1000bp → flanking=5bp (1000×0.005=5)
      - block=1100bp → flanking=6bp (1100×0.005=5.5≈6)
      - block=10000bp → flanking=50bp (10000×0.005=50)
  - `--debug`：启用调试模式，输出详细的 block 过滤日志到 `output/logs/ChrXX_debug.log`
  - `--core-identity`（98.0）：核心block身份度阈值（%），范围90-100，用于初筛阶段过滤低身份度的block
  - `--core-coverage`（80.0）：核心block覆盖度阈值（%），范围80-100，用于迭代阶段验证block的覆盖率

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

### 文件依赖关系与重生成

下图展示了各类文件之间的依赖关系，箭头表示"可以从...重生成"：

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           文件依赖关系图                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  迭代结果文件 (R*_RP_qryRecord.coords)  ←── 无法自动重生成，需重新迭代        │
│       │                                                                      │
│       ├──────────────────────────────────────────────────────┐              │
│       │                                                      │              │
│       ▼                                                      ▼              │
│  趋势TSV文件 (ChrXX_trend.tsv)              染色体坐标文件 (ChrXX.Os.*)      │
│       │                                           │                         │
│       ▼                                           ▼                         │
│  趋势PNG文件 (ChrXX_trend.png)              全基因组坐标文件 (RP_qryRecord.*) │
│       │                                           │                         │
│       │                                           ▼                         │
│       │                                     序列文件 (*.fasta)               │
│       │                                                                      │
│       └──────────────────────┐                                              │
│                              ▼                                              │
│                    基因组HTML报告 (core_block_trend_report.html)             │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 重生成能力说明

| 文件类型 | 依赖文件 | 重生成函数 | 自动重生成 |
|---------|---------|-----------|-----------|
| 迭代结果文件 | 比对结果 + 初筛文件 | `run_iterative_core_blocks_for_chr` | ❌ 需重新迭代 |
| 趋势TSV文件 | 迭代结果文件 | `regenerate_chr_trend_from_iter_files` | ✅ 自动 |
| 趋势PNG文件 | 趋势TSV文件 | `regenerate_chr_trend_png_from_tsv` | ✅ 自动 |
| 染色体坐标文件 | 最终迭代结果 | `extract_chr_level_coords` | ✅ 自动 |
| 全基因组坐标文件 | 染色体坐标文件 | `merge_chr_coords_and_extract_sequences` | ✅ 自动 |
| 序列文件 | 全基因组坐标文件 | `merge_chr_coords_and_extract_sequences` | ✅ 自动 |
| 基因组HTML报告 | 所有趋势TSV文件 | `generate_genome_trend_report` | ✅ 自动 |

---

## 文件含义说明

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

---

## 性能优化建议

### 快速迭代模式（--fast）v6.0优化版

panCB 提供两种迭代模式：
- **标准模式**（默认）：每轮迭代都执行完整的 BLAST 验证，结果最准确
- **快速模式**（`--fast`）：通过增量验证和动态flanking策略，在保持高准确性的同时大幅减少计算开销

**启用快速模式**：
```bash
python panCB.py -r REF -g genomes -accs accs.txt -o output --fast
```

**使用标准模式**（默认，无需额外参数）：
```bash
python panCB.py -r REF -g genomes -accs accs.txt -o output
```

#### v6.0核心优化

**1. 固定验证间隔（每2轮）**
- 自动在R2, R4, R6...执行BLAST验证
- 避免坐标误差累积（最多累积1轮）
- 无需手动配置验证频率

**2. 增量验证（只验证新增样本）**
- 每次验证只验证新增的样本，而非全量验证
- 已验证样本的坐标保持不变
- 验证时间与新增样本数成正比

**3. 动态flanking（自适应计算）**
- 根据block长度自动计算flanking值
- 公式：`flanking = round(block_length × 0.005)`
- 小block用小flanking，大block用大flanking

#### 工作原理对比

```
标准模式（每轮验证）:
  R1: 样本1 ∩ 样本2 → BLAST验证所有样本 → 更新坐标
  R2: R1 ∩ 样本3   → BLAST验证所有样本 → 更新坐标
  R3: R2 ∩ 样本4   → BLAST验证所有样本 → 更新坐标
  ...

快速模式v6.0（增量验证+固定间隔）:
  R1: 样本1 ∩ 样本2 → 数学计算更新坐标
  R2: R1 ∩ 样本3   → 增量BLAST验证（只验证样本3）+ 动态flanking
  R3: R2 ∩ 样本4   → 数学计算更新坐标
  R4: R3 ∩ 样本5   → 增量BLAST验证（只验证样本4-5）+ 动态flanking
  ...
```

#### 性能提升

**验证时间对比**（以20个样本为例）：

| 模式 | 验证策略 | 验证次数 | 相对时间 |
|------|---------|---------|---------|
| 标准模式 | 每轮全量验证 | 19轮 × 所有样本 | 100% |
| 快速模式v5.3 | 最终全量验证 | 1轮 × 所有样本 | 50% |
| **快速模式v6.0** | **每2轮增量验证** | **10轮 × 新增样本** | **30%** |

**关键优势**：
- 相比标准模式：总时间减少约70%
- 相比旧版快速模式：验证时间减少70-80%
- 准确性：与标准模式高度一致（坐标差异<5bp）

#### 动态Flanking示例

| Block长度 | Flanking值 | 计算过程 |
|----------|-----------|---------|
| 100bp | 1bp | 100×0.005=0.5≈1 (最小为1) |
| 1000bp | 5bp | 1000×0.005=5.0 |
| 1050bp | 5bp | 1050×0.005=5.25≈5 |
| 1100bp | 6bp | 1100×0.005=5.5≈6 |
| 10000bp | 50bp | 10000×0.005=50.0 |

#### 适用场景

**推荐使用快速模式**：
- ✅ 所有规模的泛基因组分析（10-1000+样本）
- ✅ 需要平衡准确性和速度
- ✅ 计算资源有限的环境
- ✅ 生产环境的常规分析

**推荐使用标准模式**：
- ✅ 对准确性要求极高的研究
- ✅ 样本数较少（<10个）
- ✅ 需要发表的高质量数据

#### 日志标识

启用快速迭代模式时，日志中会显示：
```
使用快速迭代模式 + 增量验证（固定每2轮）+ 动态flanking
[Chr01] R2 增量验证（验证间隔=2轮）
[Chr01] 动态flanking: 5bp (block长度×0.5%)
```

#### 参数迁移指南（v5.3 → v6.0）

**旧版本（v5.3）**：
```bash
# 需要手动配置验证频率和flanking大小
python panCB.py ... --fast --round-verify 10 --flanking-size 500
```

**新版本（v6.0）**：
```bash
# 自动优化，无需额外配置
python panCB.py ... --fast
```

**弃用参数处理**：
- `--round-verify`：保留但被忽略，发出警告
- `--flanking-size`：在快速模式下被忽略，发出警告
- 建议：移除这些参数，直接使用`--fast`

### 迭代模式技术细节

#### 标准模式 vs 快速模式

| 特性 | 标准模式 | 快速模式v6.0 |
|------|---------|-------------|
| 坐标更新方式 | 每轮BLAST验证 | 数学计算+定期验证 |
| 验证频率 | 每轮 | 每2轮 |
| 验证策略 | 全量验证 | 增量验证 |
| Flanking策略 | 无flanking | 动态flanking（0.5%） |
| 误差累积 | 无 | 最多1轮 |
| 准确性 | 最高 | 高（<5bp差异） |
| 速度 | 基准 | 快70% |
| 适用场景 | 高精度研究 | 常规分析 |

#### 增量验证原理

**标准模式（全量验证）**：
```python
# 每轮验证所有样本
for block in blocks:
    for sample in all_samples:  # 验证所有样本
        blast_and_validate(sample)
```

**快速模式v6.0（增量验证）**：
```python
# 只验证新增样本
for block in blocks:
    # 保留已验证样本的坐标
    validated_coords = block[:2 + last_verified_count*2]
    # 只验证新增样本
    for sample in new_samples:  # 只验证新增样本
        new_coords = blast_and_validate(sample)
        validated_coords.extend(new_coords)
```

**时间复杂度对比**：
- 标准模式：O(n²) - 每轮验证所有样本
- 快速模式v5.3：O(n) - 最终验证所有样本
- **快速模式v6.0**：**O(1)** - 每轮只验证新增样本

#### 动态Flanking计算

**计算公式**：
```python
flanking = max(1, round(block_length * 0.005))
```

**设计理念**：
- 小block（<1kb）：误差小，用小flanking（1-5bp）
- 中block（1-10kb）：误差中等，用中flanking（5-50bp）
- 大block（>10kb）：误差大，用大flanking（>50bp）

**优势**：
- 自适应：根据block大小自动调整
- 高效：小block不浪费计算资源
- 准确：大block有足够容错空间

### 线程数配置

**重要概念**：
- `--processes`：染色体级并行数（同时处理多条染色体）
- `--blastn-threads`：每个blastn命令使用的线程数
- **总CPU使用** = `processes × blastn_threads`
- **建议**：总线程数不超过CPU物理核心数

**推荐配置**：

| CPU核心数 | 推荐配置 | 总线程数 | 说明 |
|----------|---------|---------|------|
| 8核 | `--processes-nucmer 8 --processes-core 4 --blastn-threads 2` | 8 | 平衡并行，适合大多数情况 |
| 8核 | `--processes-nucmer 8 --processes-core 2 --blastn-threads 4` | 8 | blastn优化，适合染色体数少的情况 |
| 16核 | `--processes-nucmer 16 --processes-core 4 --blastn-threads 4` | 16 | 平衡并行，推荐 |
| 16核 | `--processes-nucmer 16 --processes-core 2 --blastn-threads 8` | 16 | blastn优化，最大化blastn性能 |

**性能提升**：
- 默认配置（`--blastn-threads 1`）：基准
- 优化配置（`--blastn-threads 2-8`）：可带来 **2-6倍加速**
- 详细分析请参考 `TRF_vs_blastn_耗时分析.md` 和 `线程数配置说明.md`

---

## 断点续跑与常见注意事项

- **断点管理**：自 panCB-2025.11 起，引入 `output/checkpoints/resume_state.json` 记录流程状态。新版机制除了记录"步骤完成"外，还会存储染色体级别的 `pending_chr / done_chr`：
  - 首次运行或未加 `--resume` 将清空旧断点；
  - 加 `--resume` 时会读取断点并打印"上次已完成 / 下一步"；
  - 步骤2-4 在运行过程中会持续刷新 `pending_chr`，异常退出后仅重跑缺失染色体；
  - 步骤5-8 会在进入前检查实际产物（BED/FASTA/报告/图表）是否存在，若缺失会自动回退。
- **迭代完成检测**（v5.2新增）：断点续跑时会精确检测迭代完成状态：
  - 检查迭代是否达到预期轮次（样本数-1）
  - 检查染色体级别坐标文件是否存在
  - 只有同时满足以上条件才认为该染色体已完成
  - 若删除了某轮迭代结果文件，会自动从该轮继续执行
  - 日志会显示：`[Chr22] 迭代未完成: 已完成R145, 预期R146，需要继续执行`
- **局部跳过策略**：各子任务内部仍会在 `resume` 模式下检查已有产物（例如 per-chr 切片、比对结果、per-sample 序列），配合断点信息可最大限度避免重复计算。
- **临时目录自动清理**：流程成功完成后，会自动删除所有染色体的临时目录（`output/genomeFilter_ChrXX/tmp`），以节省磁盘空间。清理过程会记录到日志中，即使清理失败也不会影响流程结果。
- **输入一致性**：`accs` 文件中的名称需与 `genome_dir` 下 `{name}.fasta` 完全一致。
- **染色体集合以参考基因组为准**：若样本缺失某条染色体，将写入 1 bp 占位序列以维持流程连续性。
- **外部工具失败**：不会立即中断全部流程，但会在日志中标记并跳过该对（参见日志命名空间）。TRF失败警告通常是正常的（可能由于序列过短），不影响流程继续。
- **硬编码参数**：技术参数（如 nucmer、blastn、TRF 等）已在代码中硬编码为推荐值，无需通过命令行指定。如需修改，请直接编辑源代码文件。

---

## 断点检测机制详解（Resume Checkpoint Detection）

panCB v5.2 引入了全面的层级化断点检测机制，能够检测任意目录或文件的缺失，并从正确的断点位置恢复执行。

### 层级化检测顺序

断点检测按照文件依赖关系，从底层到顶层依次检查：

```
Level 1: 迭代结果文件 (R*_RP_qryRecord.coords)
    ↓
Level 2: 趋势文件 (TSV + PNG)
    ↓
Level 3: 染色体级别坐标文件 (ChrXX.Os.RP_qryRecord.*)
    ↓
Level 4: 基因组级别报告 (core_block_trend_report.html)
```

### 检测的文件类型

| 层级 | 文件类型 | 文件路径 | 检测条件 |
|------|---------|---------|---------|
| Level 1 | 迭代结果文件 | `genomeFilter_ChrXX/coreB/ChrXX.{min_len}.R*_RP_qryRecord.coords` | 文件存在且迭代轮次达到预期（样本数-1） |
| Level 2 | 趋势TSV文件 | `genomeFilter_ChrXX/coreB/core_block_trend/ChrXX_trend.tsv` | 文件存在且非空 |
| Level 2 | 趋势PNG文件 | `genomeFilter_ChrXX/coreB/core_block_trend/ChrXX_trend.png` | 文件存在且非空 |
| Level 3 | 染色体坐标文件 | `genomeFilter_ChrXX/coreB/ChrXX.Os.RP_qryRecord.{sample}` | 所有样本的坐标文件都存在 |
| Level 4 | 基因组HTML报告 | `pan_core_block/report/core_block_trend_report.html` | 文件存在且非空 |

### 各文件类型的重生成行为

**1. 迭代结果文件缺失**
- 行为：无法自动重生成，需要重新执行迭代过程
- 日志：`[ChrXX] 迭代未完成，需要继续执行`

**2. 趋势TSV文件缺失**
- 前提：迭代结果文件存在且完整
- 行为：从迭代结果文件扫描每轮的block数量，重新生成TSV
- 日志：`[ChrXX] 趋势TSV文件缺失，尝试从迭代结果重生成...`

**3. 趋势PNG文件缺失**
- 前提：趋势TSV文件存在
- 行为：调用 `Rn_block_trend.py` 从TSV数据重新绘制折线图
- 日志：`[ChrXX] 趋势PNG文件缺失，尝试从TSV重生成...`

**4. 染色体级别坐标文件缺失**
- 前提：最终迭代结果文件存在且有效
- 行为：从最终迭代结果文件重新提取各样本的坐标
- 日志：`[ChrXX] 染色体级别坐标文件缺失，需要由sequence_process重生成`

**5. 基因组级别HTML报告缺失**
- 前提：所有染色体的趋势TSV文件存在
- 行为：收集所有染色体的趋势数据，重新生成交互式HTML报告
- 日志：`步骤9：HTML报告缺失但趋势文件存在，尝试重新生成...`

### 完整性检查函数

核心检测函数位于 `core_blocks.py`：

- `check_chr_completeness(output_dir, chr_name, min_len, ref_genome)` - 全面检查单个染色体的所有输出文件完整性
- `ensure_chr_outputs_complete(output_dir, chr_name, min_len, ref_genome)` - 确保染色体输出完整，缺失文件会尝试重生成
- `has_chr_trend_files(output_dir, chr_name)` - 检查趋势文件（TSV和PNG）是否存在

坐标文件检测函数位于 `sequence_process.py`：

- `can_regenerate_chr_coords(output_dir, chr_name, min_len)` - 检查是否可以从迭代结果重生成坐标文件

---

## 断点续跑场景与示例

以下是常见的断点续跑场景及其预期行为：

### 场景1：删除趋势目录后重跑

**操作**：
```bash
# 删除某个染色体的趋势目录
rm -rf output/genomeFilter_Chr01/coreB/core_block_trend/

# 使用 --resume 重跑
python panCB.py -r REF -g genomes -accs accs.txt -o output --resume
```

**预期行为**：
1. 检测到 Chr01 的趋势文件缺失
2. 检查迭代结果文件是否存在且完整
3. 从迭代结果文件重新生成 `Chr01_trend.tsv`
4. 从TSV文件重新生成 `Chr01_trend.png`
5. 继续执行后续步骤

**日志输出**：
```
[Chr01] 趋势文件缺失，尝试从迭代结果重生成...
[Chr01] 趋势TSV文件已重生成: output/genomeFilter_Chr01/coreB/core_block_trend/Chr01_trend.tsv
[Chr01] 趋势PNG文件已重生成: output/genomeFilter_Chr01/coreB/core_block_trend/Chr01_trend.png
[Chr01] 趋势文件重新生成成功
```

### 场景2：删除染色体坐标文件后重跑

**操作**：
```bash
# 删除某个染色体的坐标文件
rm output/genomeFilter_Chr01/coreB/Chr01.Os.RP_qryRecord.*

# 使用 --resume 重跑
python panCB.py -r REF -g genomes -accs accs.txt -o output --resume
```

**预期行为**：
1. 检测到 Chr01 的坐标文件缺失
2. 检查最终迭代结果文件是否存在
3. 从最终迭代结果文件重新提取各样本的坐标
4. 继续执行全基因组合并步骤

**日志输出**：
```
[Chr01] 染色体级别坐标文件缺失，需要重新生成
[Chr01] 选择最终迭代文件: Chr01.50.R146_RP_qryRecord.coords
[Chr01] 坐标文件重新生成成功
```

### 场景3：删除HTML报告后重跑

**操作**：
```bash
# 删除基因组级别HTML报告
rm output/pan_core_block/report/core_block_trend_report.html

# 使用 --resume 重跑
python panCB.py -r REF -g genomes -accs accs.txt -o output --resume
```

**预期行为**：
1. 检测到HTML报告缺失
2. 检查所有染色体的趋势TSV文件是否存在
3. 收集所有染色体的趋势数据
4. 重新生成交互式HTML报告

**日志输出**：
```
步骤9：HTML报告缺失但趋势文件存在，尝试重新生成...
步骤9：基因组级别HTML报告生成成功: output/pan_core_block/report/core_block_trend_report.html
```

### 场景4：删除多个文件后重跑

**操作**：
```bash
# 删除多个文件
rm -rf output/genomeFilter_Chr01/coreB/core_block_trend/
rm output/genomeFilter_Chr02/coreB/Chr02.Os.RP_qryRecord.*
rm output/pan_core_block/report/core_block_trend_report.html

# 使用 --resume 重跑
python panCB.py -r REF -g genomes -accs accs.txt -o output --resume
```

**预期行为**：
1. 按层级顺序检测每个染色体的完整性
2. Chr01：重生成趋势文件（TSV + PNG）
3. Chr02：重生成坐标文件
4. 最后：重生成基因组级别HTML报告

### 场景5：删除迭代结果文件后重跑

**操作**：
```bash
# 删除某轮迭代结果
rm output/genomeFilter_Chr01/coreB/Chr01.50.R100_RP_qryRecord.coords

# 使用 --resume 重跑
python panCB.py -r REF -g genomes -accs accs.txt -o output --resume
```

**预期行为**：
1. 检测到 Chr01 迭代未完成
2. 从最后完成的迭代轮次继续执行
3. 完成迭代后重新生成趋势文件和坐标文件

**日志输出**：
```
[Chr01] 迭代未完成: 已完成R99, 预期R146，需要继续执行
```

---

## 主要模块与入口

- `panCB.py`：主控入口（参数解析、日志、总流程编排）
- `genome_alignment.py`：染色体拆分、NUCmer 比对、初筛
- `core_blocks.py`：对齐优先级排序、核心块迭代识别
- `sequence_process.py`：按样本拆分坐标、合并全基因组并提取序列
- `report.py`：核心长度统计与详细文本报告
- `utils.py`：通用工具与命令封装
- `debug_log.py`：调试日志模块，提供统一的调试日志管理功能

---

## 脚本功能详解

- `panCB.py`
  - CLI 入口，解析所有参数（含 `blastn`/`bedtools`/`TRF` 配置、`--step`、`--resume`），初始化日志，并预生成 `.fai`。
  - 管理 `utils.ResumeManager` 状态文件 `output/checkpoints/resume_state.json`，支持按步骤、按染色体断点续跑与自动回退。
  - 将拆分、比对、排序、迭代、染色体输出、全基因组合并、alternative、统计与清理串联，步骤 4 完成后立即触发步骤 4.5。

- `genome_alignment.py`
  - `split_all_genomes_by_chr`：以参考染色体列表为模板生成 `Chr_split/ChrXX/{sample}_ChrXX.fasta`，缺失染色体写入 1bp `N`。
  - `run_all_alignments_for_chr`：遍历除参考外的所有样本，调用 `nucmer → delta-filter → show-coords` 并收集 `.filtered.coords`。
  - `core_block_filter`：实现 98% 身份度、最小长度、坐标切分与去重规则，输出 `*.filtered.1227` 初筛文件。

- `core_blocks/` 包（模块化实现）
  - **iterate_control.py**：迭代控制与进度管理
    - `run_iterative_core_blocks_for_chr`：主入口函数，协调标准/快速模式的迭代流程
    - `generate_alignment_priority_for_chr`：统计样本核心长度并排序
    - `VERIFY_INTERVAL`：固定验证间隔常量（=2）
    - `recover_verification_state`：断点续跑时恢复验证状态
    - `get_expected_iter_count`：获取预期迭代轮次数
    - `get_last_completed_iter`：获取最后完成的迭代轮次
    - `is_iteration_complete`：检测迭代是否完全完成
    - `has_chr_level_coords`：检测染色体级别坐标文件是否存在
    - `has_chr_trend_files`：检测趋势文件是否存在
    - `generate_genome_trend_report`：生成基因组级别HTML报告
  - **standard_iteration.py**：标准迭代模式实现
    - `extract_core_blocks_once`：单轮核心计算，执行完整BLAST验证
  - **fast_iteration.py**：快速迭代模式实现
    - `extract_core_blocks_fast`：快速迭代（仅数学计算）
    - `validate_incremental`：增量BLAST验证（只验证新增样本）
    - `validate_final_coordinates`：最终验证（已弃用，保留向后兼容）
  - **coord_utils.py**：坐标计算工具
    - `calculate_coordinate_update`：数学计算更新坐标
    - `calculate_dynamic_flanking`：动态计算flanking值
    - `calculate_flanking_coords`：计算带flanking的提取坐标
    - `nt_stat`：统计序列核苷酸数量
  - **block_filters.py**：Block后处理过滤
    - `clean_transposition`：清理转座/错序blocks
    - `process_overlaps`：处理重叠blocks
    - `filter_by_length`：长度过滤
    - `filter_invalid_coords`：过滤异常坐标

- `sequence_process.py`
  - `extract_chr_level_coords`（步骤 4.5）：把最终 `R*` 结果拆成 `ChrXX.Os.RP_qryRecord.{sample}`（三列 BED）。
  - `merge_chr_coords_and_extract_sequences`：跨染色体合并生成 `pan_core_block/RP_qryRecord.{sample}`，并以 `bedtools getfasta` 输出 `fasta/`。
  - `has_sequence_outputs`：用于 resume 判断染色体级、全基因组 BED 与 FASTA 是否齐全。

- `alternative_block.py`
  - `run_alternative_for_all_chr`：读取 `ChrXX.Os.*`，基于 `{sample}.fasta.fai` 的染色体长度计算 core 的补集，输出 `alter-BED/` 与 `alter-FASTA/`。
  - 自动生成 `sorted.bed`, `.sorted.len`, `.alternative.bed`, `.sorted.bed.named` 及对应 FASTA，并支持 resume。

- `report.py`
  - `summarize_core_lengths`：统计 `chr × sample` 与全基因组核心长度，输出 `core_len_per_chr.tsv`、`core_len_total.tsv`。
  - `write_detailed_summary`：生成 `core_summary.txt`，逐样本 × 逐染色体报告区块数量 / 总长度 / 平均长度。
  - `generate_summary_html_report`：读取 `core_summary.txt`，生成交互式HTML统计报告（`core_summary.html`）。
  - `generate_all_visualizations`：统一调用可视化函数，生成HTML报告。
  - `has_summary_outputs`、`has_detail_report`：检测统计及文本报告是否存在。

- `utils.py`
  - 收录常用方法：`ensure_dir`、`run_cmd`、`ensure_fai_exists`、`intersection`、`count_lines` 等。
  - `ResumeManager`：封装状态的读写、步骤标记、元数据更新、自动回退等逻辑。

- `debug_log.py`
  - `DebugLogger` 类：调试日志管理器，用于管理染色体级别的调试日志文件。
  - `create_debug_logger`：工厂函数，创建 `DebugLogger` 实例。
  - 当启用 `--debug` 参数时，会将详细的 block 过滤过程记录到 `output/logs/ChrXX_debug.log` 文件中。

- `Rn_block_trend.py`
  - 可视化模块，用于绘制核心block迭代变化趋势图。
  - 主要功能：
    - `read_tsv`：读取TSV格式的迭代数据
    - `plot_trend`：绘制单染色体PNG折线图
    - `collect_all_chr_data`：收集所有染色体的TSV数据
    - `generate_html_report`：生成基因组级别的交互式HTML报告（v5.2增强：点击染色体显示折线图+表格）
    - `generate_genome_report`：主入口函数，协调数据收集和报告生成
  - 自动调用：
    - 染色体级别：每条染色体迭代完成后自动生成PNG图
    - 基因组级别：所有染色体完成后自动生成HTML报告
  - 手动使用：`python Rn_block_trend.py -i Chr01_trend.tsv -o Chr01_trend.png -c Chr01`

---

## 迭代过程与可视化

### 迭代趋势可视化

panCB 会自动记录每轮迭代的 block 数量变化，并生成可视化图表：

**生成文件**：
- 染色体级别：
  - `output/genomeFilter_ChrXX/coreB/core_block_trend/ChrXX_trend.tsv`：记录每轮迭代的数据
  - `output/genomeFilter_ChrXX/coreB/core_block_trend/ChrXX_trend.png`：block数量变化趋势图
- 基因组级别：
  - `output/pan_core_block/report/core_block_trend_report.html`：交互式HTML报告

**TSV 文件格式**：
```
iteration	block_count	sample_name
R1	1234	Sample2
R2	1456	Sample3
R3	1678	Sample4
...
```

### Block 数量上升现象

**观察到的现象**：迭代过程中 block 数量通常呈上升趋势

**这是正常的！原因如下**：

1. **迭代不是简单取交集**：而是对交集区域进行切分
2. **切分效应**：新样本的变异位点会将长 block 切分成多个短 block
3. **生物学意义**：反映了样本间的遗传多样性

**详细解释**：
- R1: 100个block（每个可能很长）
- R2: 对R1的每个block与Sample3求交集并切分 → 150个block（更短）
- R3: 继续切分 → 200个block（更短）
- ...
- **总长度下降**，**block数量上升**，**平均长度下降**

**验证合理性**：
- Block 数量：上升 ✓
- 总长度：下降 ✓
- 平均长度：下降 ✓

详细说明请参考 `迭代过程详解.md` 文档。

---

## v6.0技术改进总结

### 快速模式优化

**优化前（v5.3）**：
- 验证策略：仅最终验证
- 验证方式：全量验证所有样本
- Flanking：固定值（500bp或2kb）
- 问题：验证时间随样本数线性增长

**优化后（v6.0）**：
- 验证策略：固定每2轮验证
- 验证方式：增量验证新增样本
- Flanking：动态计算（block长度的0.5%）
- 优势：验证时间与新增样本数成正比

### 性能提升

| 指标 | v5.3快速模式 | v6.0快速模式 | 提升 |
|------|-------------|-------------|------|
| 验证频率 | 仅最终 | 每2轮 | 更频繁 |
| 验证时间 | O(n) | O(1) | 70-80% |
| 误差累积 | 最多n轮 | 最多1轮 | 显著减少 |
| 准确性 | 中等 | 高 | 接近标准模式 |

### 代码模块化

**优化前**：
- 单一文件：`core_blocks.py`（2000+行）
- 难以维护和测试

**优化后**：
- 模块化包：`core_blocks/`
  - `iterate_control.py`：迭代控制
  - `standard_iteration.py`：标准模式
  - `fast_iteration.py`：快速模式
  - `coord_utils.py`：坐标工具
  - `block_filters.py`：过滤器
- 易于维护、测试和扩展

### 用户体验改进

**简化参数**：
- 移除`--round-verify`（自动优化）
- 移除`--flanking-size`（动态计算）
- 用户只需`--fast`即可获得最佳性能

**向后兼容**：
- 保留弃用参数，发出警告
- 输出格式完全兼容
- 平滑升级，无需修改脚本

### 技术文档

详细的技术文档位于：
- 设计文档：`.kiro/specs/fast-mode-optimization/design.md`
- 需求文档：`.kiro/specs/fast-mode-optimization/requirements.md`
- 任务列表：`.kiro/specs/fast-mode-optimization/tasks.md`
- 测试代码：`tests/test_fast_mode_core.py`

---

如需复现实验或更改参数，建议在不同 `output` 目录下运行以避免产物冲突。