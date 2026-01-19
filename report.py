import os
import logging
import json
from typing import List, Dict, Any
from collections import defaultdict
import utils


def summarize_core_lengths(output_dir: str, chr_list: List[str], ref_genome: str, resume: bool = False) -> str:
    """
    统计最终核心block长度：
    - 对每条染色体，读取 genomeFilter_ChrXX/coreB/ChrXX.Os.RP_qryRecord.{name}
    - 统计各样本在该染色体上的总核心长度
    - 输出：
      output/pan_core_block/report/core_len_per_chr.tsv (chr, sample, length)
      output/pan_core_block/report/core_len_total.tsv (sample, total_length)
    返回 report 目录路径。
    """
    logger = logging.getLogger("panCB.report")
    report_dir = os.path.join(output_dir, "pan_core_block", "report")
    os.makedirs(report_dir, exist_ok=True)

    per_chr_path = os.path.join(report_dir, "core_len_per_chr.tsv")
    total_path = os.path.join(report_dir, "core_len_total.tsv")

    if resume and os.path.isdir(report_dir) and os.path.isfile(per_chr_path) and os.path.isfile(total_path):
        logger.info("[resume] 检测到已有 report 目录，跳过统计: %s", report_dir)
        return report_dir

    # 收集样本集合
    names_set = set()
    for chr_name in chr_list:
        chr_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}", "coreB")
        if not os.path.isdir(chr_dir):
            continue
        for fname in os.listdir(chr_dir):
            if fname.startswith(f"{chr_name}.Os.RP_qryRecord."):
                names_set.add(fname.split(".Os.RP_qryRecord.")[-1])

    # 统计 per-chr
    total_len: Dict[str, int] = {}
    with open(per_chr_path, "w") as fo:
        fo.write("chr\tsample\tcore_len\n")
        for chr_name in chr_list:
            chr_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}", "coreB")
            for name in sorted(list(names_set)):
                bed = os.path.join(chr_dir, f"{chr_name}.Os.RP_qryRecord.{name}")
                lens = 0
                if os.path.isfile(bed):
                    with open(bed) as f:
                        for line in f:
                            cols = line.strip().split()
                            if len(cols) >= 3:
                                try:
                                    lens += int(cols[2]) - int(cols[1])
                                except ValueError:
                                    pass
                fo.write(f"{chr_name}\t{name}\t{lens}\n")
                total_len[name] = total_len.get(name, 0) + lens

    with open(total_path, "w") as fo:
        fo.write("sample\ttotal_core_len\n")
        for name in sorted(total_len.keys()):
            fo.write(f"{name}\t{total_len[name]}\n")

    logger.info("统计完成: %s(%dB), %s(%dB)", per_chr_path, utils.file_size(per_chr_path), total_path, utils.file_size(total_path))
    return report_dir


def has_summary_outputs(output_dir: str) -> bool:
    """
    检查统计报告文件是否已生成。
    
    Args:
        output_dir: 输出根目录
    
    Returns:
        如果core_len_per_chr.tsv和core_len_total.tsv都存在且非空，返回True
    """
    report_dir = os.path.join(output_dir, "pan_core_block", "report")
    per_chr_path = os.path.join(report_dir, "core_len_per_chr.tsv")
    total_path = os.path.join(report_dir, "core_len_total.tsv")
    return utils.file_size(per_chr_path) > 0 and utils.file_size(total_path) > 0


def write_detailed_summary(output_dir: str, chr_list: List[str], ref_genome: str, resume: bool = False) -> str:
    """
    输出详细汇总报告（文本，UTF-8）：
    格式示例：
    泛基因组核心block分析汇总报告

    参考基因组: ref

    分析样本数: 3

    分析样本: query2, query1, ref

    核心区块统计:

      query2 - Chr25:

        核心区块统计信息

        参考基因组: ref

        比较基因组: query2

        染色体: Chr25

        区块数量: 187

        总长度: 18879483 bp

        平均长度: 100959.80 bp
    """
    logger = logging.getLogger("panCB.report")
    report_dir = os.path.join(output_dir, "pan_core_block", "report")
    os.makedirs(report_dir, exist_ok=True)
    summary_path = os.path.join(report_dir, "core_summary.txt")

    if resume and os.path.isdir(report_dir) and os.path.isfile(summary_path):
        logger.info("[resume] 检测到已有详细报告，跳过: %s", summary_path)
        return summary_path

    # 收集样本集合
    names_set = set()
    for chr_name in chr_list:
        chr_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}", "coreB")
        if not os.path.isdir(chr_dir):
            logger.warning("[%s] 缺少核心目录: %s", chr_name, chr_dir)
            continue
        for fname in os.listdir(chr_dir):
            if fname.startswith(f"{chr_name}.Os.RP_qryRecord."):
                names_set.add(fname.split(".Os.RP_qryRecord.")[-1])
    
    if not names_set:
        logger.warning("未找到任何 ChrXX.Os.RP_qryRecord.* 文件，无法生成详细报告。请确保步骤4.5已正确执行。")
        return summary_path

    names = sorted(list(names_set))

    with open(summary_path, "w", encoding="utf-8") as fo:
        fo.write("泛基因组核心block分析汇总报告\n\n")
        fo.write(f"参考基因组: {ref_genome}\n\n")
        fo.write(f"分析样本数: {len(names)}\n\n")
        fo.write("分析样本: " + ", ".join(names) + "\n\n")
        fo.write("核心区块统计:\n\n")

        for sample in names:
            for chr_name in chr_list:
                bed = os.path.join(output_dir, f"genomeFilter_{chr_name}", "coreB", f"{chr_name}.Os.RP_qryRecord.{sample}")
                block_count = 0
                total_len = 0
                if os.path.isfile(bed):
                    with open(bed) as f:
                        for line in f:
                            cols = line.strip().split()
                            if len(cols) >= 3:
                                try:
                                    block_count += 1
                                    total_len += int(cols[2]) - int(cols[1])
                                except ValueError:
                                    pass
                avg_len = (float(total_len) / block_count) if block_count > 0 else 0.0

                fo.write(f"  {sample} - {chr_name}:\n\n")
                fo.write("    核心区块统计信息\n\n")
                fo.write(f"    参考基因组: {ref_genome}\n\n")
                fo.write(f"    比较基因组: {sample}\n\n")
                fo.write(f"    染色体: {chr_name}\n\n")
                fo.write(f"    区块数量: {block_count}\n\n")
                fo.write(f"    总长度: {total_len} bp\n\n")
                fo.write(f"    平均长度: {avg_len:.2f} bp\n\n")

    logger.info("详细汇总报告: %s(%dB)", summary_path, utils.file_size(summary_path))
    return summary_path


def has_detail_report(output_dir: str) -> bool:
    report_dir = os.path.join(output_dir, "pan_core_block", "report")
    summary_path = os.path.join(report_dir, "core_summary.txt")
    return utils.file_size(summary_path) > 0


def parse_core_summary_txt(txt_path: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    解析core_summary.txt文件，按染色体分组整理数据。
    
    Args:
        txt_path: TXT文件路径
    
    Returns:
        字典，键为染色体名称，值为样本数据列表
    """
    chr_data = defaultdict(list)
    current_sample = None

    try:
        with open(txt_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # 识别样本-染色体标识行（如"HMXRS1 - Chr01:"）
                if ' - Chr' in line and line.endswith(':'):
                    sample_chr = line[:-1].split(' - ')
                    if len(sample_chr) == 2:
                        compare_genome = sample_chr[0].strip()
                        chromosome = sample_chr[1].strip()
                        current_sample = {
                            "比较基因组": compare_genome,
                            "染色体": chromosome,
                            "区块数量": None,
                            "总长度（bp）": None,
                            "平均长度（bp）": None
                        }
                    continue

                # 提取样本详情字段
                if current_sample is not None and ': ' in line:
                    key, value = [item.strip() for item in line.split(': ', 1)]
                    if key == "区块数量":
                        current_sample["区块数量"] = int(value)
                    elif key == "总长度":
                        current_sample["总长度（bp）"] = int(value.replace(' bp', ''))
                    elif key == "平均长度":
                        current_sample["平均长度（bp）"] = round(float(value.replace(' bp', '')), 2)
                        # 保存完整的样本数据
                        if all(v is not None for v in current_sample.values()):
                            chr_name = current_sample["染色体"]
                            chr_data[chr_name].append(current_sample.copy())
                        current_sample = None

    except Exception as e:
        logging.getLogger("panCB.report").error("解析TXT文件失败: %s", e)

    return dict(chr_data)


def generate_summary_html_report(output_dir: str, chr_list: List[str], ref_genome: str, resume: bool = False) -> str:
    """
    读取core_summary.txt文件，生成交互式HTML报告。
    
    Args:
        output_dir: 输出根目录
        chr_list: 染色体列表
        ref_genome: 参考基因组名称
        resume: 是否启用断点续跑
    
    Returns:
        生成的HTML文件路径，失败返回空字符串
    """
    logger = logging.getLogger("panCB.report")
    
    report_dir = os.path.join(output_dir, "pan_core_block", "report")
    txt_path = os.path.join(report_dir, "core_summary.txt")
    html_path = os.path.join(report_dir, "core_summary.html")
    
    if resume and os.path.isfile(html_path) and utils.file_size(html_path) > 0:
        logger.info("[resume] 跳过已有HTML报告: %s", html_path)
        return html_path
    
    if not os.path.isfile(txt_path):
        logger.warning("TXT文件不存在: %s", txt_path)
        return ""
    
    # 解析TXT文件
    chr_data = parse_core_summary_txt(txt_path)
    
    if not chr_data:
        logger.warning("未从TXT文件中解析到有效数据")
        return ""
    
    # 准备图表数据
    chromosomes = sorted(chr_data.keys())
    # 计算每个染色体的总区块数（所有样本的区块数之和的平均值，或取第一个样本）
    block_counts = []
    for chr_name in chromosomes:
        if chr_data[chr_name]:
            # 取该染色体下所有样本的区块数量的平均值
            counts = [s["区块数量"] for s in chr_data[chr_name] if s["区块数量"] is not None]
            block_counts.append(int(sum(counts) / len(counts)) if counts else 0)
        else:
            block_counts.append(0)
    
    chart_datasets = [
        {
            "label": "核心区块数量",
            "data": block_counts,
            "backgroundColor": "rgba(52, 152, 219, 0.6)",
            "borderColor": "rgba(52, 152, 219, 1)",
            "borderWidth": 1,
            "borderRadius": 4
        }
    ]
    
    # 生成HTML内容
    chr_blocks_html = ''.join([f'<div class="chr-block" data-chr="{chr_name}">{chr_name}</div>' for chr_name in chromosomes])
    
    html_content = f'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>泛基因组核心block统计分析结果</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{
            font-family: "Arial", "Microsoft YaHei", sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            line-height: 1.6;
            background-color: #f5f7fa;
        }}
        .section {{
            margin-bottom: 30px;
            padding: 20px;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            background: white;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
            border-bottom: 1px solid #e0e0e0;
            padding-bottom: 10px;
        }}
        .header-info {{
            background: #ecf0f1;
            padding: 15px;
            border-radius: 6px;
            margin-bottom: 20px;
        }}
        .header-info p {{
            margin: 5px 0;
        }}
        #chart-container {{
            width: 100%;
            height: 400px;
            margin-top: 20px;
        }}
        .chr-blocks {{
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
            margin-top: 20px;
        }}
        .chr-block {{
            padding: 15px 25px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 6px;
            cursor: pointer;
            transition: all 0.3s ease;
            font-weight: bold;
            box-shadow: 0 2px 5px rgba(0,0,0,0.2);
        }}
        .chr-block:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 10px rgba(0,0,0,0.3);
        }}
        .chr-block.active {{
            background: linear-gradient(135deg, #e74c3c 0%, #c0392b 100%);
        }}
        #detail-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            text-align: center;
        }}
        #detail-table th, #detail-table td {{
            padding: 12px 15px;
            border: 1px solid #e0e0e0;
        }}
        #detail-table th {{
            background-color: #34495e;
            color: white;
            font-weight: bold;
        }}
        #detail-table tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        #detail-table tr:hover {{
            background-color: #e8f4f8;
        }}
        #table-placeholder {{
            margin-top: 20px;
            padding: 40px;
            text-align: center;
            color: #7f8c8d;
            background-color: #f8f9fa;
            border-radius: 6px;
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            color: #95a5a6;
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <h1>🧬 泛基因组核心block统计分析报告</h1>
    
    <div class="section">
        <div class="header-info">
            <p><strong>参考基因组:</strong> {ref_genome}</p>
            <p><strong>分析染色体数:</strong> {len(chromosomes)}</p>
        </div>
    </div>

    <div class="section">
        <h2>📊 各染色体核心区块数量统计</h2>
        <div id="chart-container">
            <canvas id="blockCountChart"></canvas>
        </div>
    </div>

    <div class="section">
        <h2>🔬 染色体详情查询</h2>
        <h3>请点击下方染色体矩形块查看详情</h3>
        <div class="chr-blocks">
            {chr_blocks_html}
        </div>
        <div id="table-container">
            <div id="table-placeholder">👆 未选择染色体，请点击上方染色体矩形块</div>
        </div>
    </div>

    <div class="footer">
        Generated by panCB Pipeline - Core Block Statistics Reporter
    </div>

    <script>
        const ctx = document.getElementById('blockCountChart').getContext('2d');
        new Chart(ctx, {{
            type: 'bar',
            data: {{
                labels: {json.dumps(chromosomes)},
                datasets: {json.dumps(chart_datasets, ensure_ascii=False)}
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        beginAtZero: true,
                        title: {{
                            display: true,
                            text: '核心区块数量（个）',
                            font: {{ size: 14 }}
                        }},
                        ticks: {{ precision: 0 }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: '染色体',
                            font: {{ size: 14 }}
                        }}
                    }}
                }},
                plugins: {{
                    legend: {{ display: false }},
                    tooltip: {{
                        callbacks: {{
                            label: function(context) {{
                                return '区块数量：' + context.raw + ' 个';
                            }}
                        }}
                    }}
                }}
            }}
        }});

        const chrData = {json.dumps(chr_data, ensure_ascii=False)};
        const chrBlocks = document.querySelectorAll('.chr-block');
        const tableContainer = document.getElementById('table-container');

        chrBlocks.forEach(block => {{
            block.addEventListener('click', function() {{
                chrBlocks.forEach(b => b.classList.remove('active'));
                this.classList.add('active');
                
                const currentChr = this.getAttribute('data-chr');
                const sampleList = chrData[currentChr];
                
                if (!sampleList || sampleList.length === 0) {{
                    tableContainer.innerHTML = '<div id="table-placeholder">该染色体暂无数据</div>';
                    return;
                }}
                
                let tableHtml = '<h3>' + currentChr + ' 染色体核心block统计详情</h3>';
                tableHtml += '<table id="detail-table">';
                tableHtml += '<thead><tr><th>比较基因组</th><th>染色体</th><th>区块数量（个）</th><th>总长度（bp）</th><th>平均长度（bp）</th></tr></thead>';
                tableHtml += '<tbody>';
                
                sampleList.forEach(sample => {{
                    tableHtml += '<tr>';
                    tableHtml += '<td>' + sample['比较基因组'] + '</td>';
                    tableHtml += '<td>' + sample['染色体'] + '</td>';
                    tableHtml += '<td>' + sample['区块数量'] + '</td>';
                    tableHtml += '<td>' + sample['总长度（bp）'].toLocaleString() + '</td>';
                    tableHtml += '<td>' + sample['平均长度（bp）'].toFixed(2) + '</td>';
                    tableHtml += '</tr>';
                }});
                
                tableHtml += '</tbody></table>';
                tableContainer.innerHTML = tableHtml;
            }});
        }});
    </script>
</body>
</html>'''
    
    try:
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        logger.info("HTML报告已生成: %s", html_path)
        return html_path
    except Exception as e:
        logger.error("生成HTML报告失败: %s", e)
        return ""


def generate_all_visualizations(output_dir: str, chr_list: List[str], ref_genome: str, resume: bool = False) -> bool:
    """
    生成可视化报告（HTML交互报告）。
    
    Args:
        output_dir: 输出根目录
        chr_list: 染色体列表
        ref_genome: 参考基因组名称
        resume: 是否启用断点续跑
    
    Returns:
        成功返回True，失败返回False
    """
    logger = logging.getLogger("panCB.report")
    
    # 生成HTML交互报告
    logger.info("生成HTML交互报告...")
    html_path = generate_summary_html_report(output_dir, chr_list, ref_genome, resume)
    if not html_path:
        logger.warning("HTML交互报告生成失败")
        return False
    
    return True