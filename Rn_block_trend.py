#!/usr/bin/env python3
"""
Rn_block_trend.py - 绘制核心block迭代变化趋势折线图

功能：读取TSV文件，绘制从R1到Rn迭代过程中核心block数量的变化趋势
（block数量可能上升或下降）

输入TSV格式：
iteration\tblock_count\tsample_name
R1\t1234\tSample1
R2\t1100\tSample2
...

用法：
    python Rn_block_trend.py -i input.tsv -o output.png [-t "Title"]
"""

import argparse
import logging
import os
import sys
from typing import List, Tuple, Optional, Dict

import matplotlib
matplotlib.use('Agg')  # 非交互式后端，适用于服务器环境
import matplotlib.pyplot as plt

# 获取模块级别的logger
logger = logging.getLogger(__name__)


def read_tsv(tsv_path: str) -> Tuple[List[str], List[int], List[str]]:
    """
    读取TSV文件，返回迭代轮次列表和对应的block数量列表
    
    Args:
        tsv_path: TSV文件路径
    
    Returns:
        (iterations, block_counts, sample_names): 三个列表
        如果文件不存在或读取失败，返回三个空列表
    
    Raises:
        无异常抛出，所有错误通过日志记录并返回空列表
    """
    iterations: List[str] = []
    block_counts: List[int] = []
    sample_names: List[str] = []
    
    # 检查文件是否存在
    if not os.path.isfile(tsv_path):
        logger.warning("TSV文件不存在: %s", tsv_path)
        return iterations, block_counts, sample_names
    
    try:
        with open(tsv_path, 'r', encoding='utf-8') as f:
            # 读取并验证表头
            header = f.readline()
            if not header:
                logger.warning("TSV文件为空: %s", tsv_path)
                return iterations, block_counts, sample_names

            header = header.strip()
            # 验证表头格式（可选，但有助于检测格式问题）
            expected_header_parts = ['iteration', 'block_count', 'sample_name']
            header_parts = header.lower().split('\t')
            if len(header_parts) < 2:
                logger.warning("TSV表头格式不正确，期望至少2列: %s", tsv_path)
            
            line_num = 1
            for line in f:
                line_num += 1
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) < 2:
                    logger.warning("第%d行格式不正确，跳过: %s", line_num, line)
                    continue
                
                # 解析迭代轮次
                iteration = parts[0].strip()
                if not iteration:
                    logger.warning("第%d行迭代轮次为空，跳过", line_num)
                    continue
                
                # 解析block数量
                try:
                    block_count = int(parts[1].strip())
                    if block_count < 0:
                        logger.warning("第%d行block数量为负数(%d)，跳过", line_num, block_count)
                        continue
                except ValueError:
                    logger.warning("第%d行block数量无法解析为整数: %s", line_num, parts[1])
                    continue
                
                # 解析样本名称（可选）
                sample_name = parts[2].strip() if len(parts) >= 3 else ''
                
                iterations.append(iteration)
                block_counts.append(block_count)
                sample_names.append(sample_name)
        
        logger.debug("成功读取TSV文件: %s, 共%d条记录", tsv_path, len(iterations))
        
    except IOError as e:
        logger.error("读取TSV文件失败: %s, 错误: %s", tsv_path, str(e))
    except Exception as e:
        logger.error("解析TSV文件时发生未知错误: %s, 错误: %s", tsv_path, str(e))
    
    return iterations, block_counts, sample_names


def plot_trend(
    iterations: List[str],
    block_counts: List[int],
    sample_names: List[str],
    output_path: str,
    title: Optional[str] = None,
    chr_name: Optional[str] = None
) -> bool:
    """
    绘制核心block变化趋势折线图
    
    Args:
        iterations: 迭代轮次列表 ['R1', 'R2', ...]
        block_counts: block数量列表
        sample_names: 样本名称列表
        output_path: 输出图片路径
        title: 图表标题（可选）
        chr_name: 染色体名称（可选）
    
    Returns:
        成功返回True，失败返回False
    """
    # 验证输入数据
    if not iterations or not block_counts:
        logger.warning("没有数据可绘制")
        return False
    
    if len(iterations) != len(block_counts):
        logger.error("迭代轮次数量(%d)与block数量(%d)不匹配", 
                    len(iterations), len(block_counts))
        return False
    
    # 确保输出目录存在
    out_dir = os.path.dirname(output_path)
    if out_dir:
        try:
            os.makedirs(out_dir, exist_ok=True)
        except OSError as e:
            logger.error("创建输出目录失败: %s, 错误: %s", out_dir, str(e))
            return False
    
    try:
        # 设置图表大小
        fig_width = max(10, len(iterations) * 0.3)
        fig, ax = plt.subplots(figsize=(fig_width, 6))
        
        # 绘制折线图
        x = range(len(iterations))
        ax.plot(x, block_counts, marker='o', linewidth=2, markersize=6, color='#2E86AB')
        ax.fill_between(x, block_counts, alpha=0.3, color='#2E86AB')
        
        # 设置x轴
        ax.set_xticks(list(x))
        # 如果迭代次数太多，只显示部分标签
        if len(iterations) > 30:
            step = max(1, len(iterations) // 20)
            labels = [iterations[i] if i % step == 0 else '' for i in range(len(iterations))]
            ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        else:
            ax.set_xticklabels(iterations, rotation=45, ha='right', fontsize=9)
        
        # 设置标签
        ax.set_xlabel('Iteration Round', fontsize=12)
        ax.set_ylabel('Core Block Count', fontsize=12)
        
        # 设置标题
        if title:
            ax.set_title(title, fontsize=14, fontweight='bold')
        elif chr_name:
            ax.set_title(f'Core Block Trend - {chr_name}', fontsize=14, fontweight='bold')
        else:
            ax.set_title('Core Block Trend During Iteration', fontsize=14, fontweight='bold')
        
        # 添加网格
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # 添加起始和结束点的数值标注
        if len(block_counts) > 0:
            ax.annotate(f'{block_counts[0]}', xy=(0, block_counts[0]), 
                       xytext=(5, 10), textcoords='offset points', fontsize=10, color='green')
        if len(block_counts) > 1:
            ax.annotate(f'{block_counts[-1]}', xy=(len(block_counts)-1, block_counts[-1]), 
                       xytext=(5, 10), textcoords='offset points', fontsize=10, color='red')
        
        # 添加变化率信息
        if len(block_counts) >= 2 and block_counts[0] > 0:
            change_rate = ((block_counts[-1] - block_counts[0]) / block_counts[0]) * 100
            change_label = 'Change'
            info_text = f'Initial: {block_counts[0]}\nFinal: {block_counts[-1]}\n{change_label}: {change_rate:+.1f}%'
            ax.text(0.98, 0.98, info_text, transform=ax.transAxes, fontsize=10,
                   verticalalignment='top', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        
        # 验证文件是否成功创建
        if os.path.isfile(output_path) and os.path.getsize(output_path) > 0:
            logger.info("图表已保存: %s", output_path)
            return True
        else:
            logger.error("图表文件创建失败或为空: %s", output_path)
            return False
            
    except Exception as e:
        logger.error("绘制图表时发生错误: %s", str(e))
        # 确保关闭任何打开的图形
        plt.close('all')
        return False


def collect_all_chr_data(output_dir: str, chr_list: List[str]) -> Dict[str, Dict]:
    """
    收集所有染色体的TSV数据
    
    扫描输出目录中的染色体TSV文件，解析每个文件并提取数据。
    对于缺失或无效的文件，记录警告并跳过该染色体。
    
    Args:
        output_dir: 输出根目录（如 'output'）
        chr_list: 染色体名称列表（如 ['Chr01', 'Chr02', ...]）
    
    Returns:
        字典，键为染色体名称，值为包含以下字段的字典：
        - iterations: 迭代轮次列表 ['R1', 'R2', ...]
        - block_counts: block数量列表 [1234, 1100, ...]
        - sample_names: 样本名称列表 ['Sample1', 'Sample2', ...]
        - final_count: 最终block数量（最后一轮的值）
        - initial_count: 初始block数量（第一轮的值）
        - change_rate: 变化率 ((final - initial) / initial * 100)
        
        如果某染色体数据无效或文件不存在，该染色体不会出现在返回字典中。
    
    Example:
        >>> chr_data = collect_all_chr_data('output', ['Chr01', 'Chr02'])
        >>> chr_data['Chr01']['final_count']
        1678
    """
    chr_data: Dict[str, Dict] = {}
    
    for chr_name in chr_list:
        # 构建TSV文件路径: output/genomeFilter_{chr_name}/coreB/core_block_trend/{chr_name}_trend.tsv
        tsv_path = os.path.join(
            output_dir,
            f"genomeFilter_{chr_name}",
            "coreB",
            "core_block_trend",
            f"{chr_name}_trend.tsv"
        )
        
        # 读取TSV数据
        iterations, block_counts, sample_names = read_tsv(tsv_path)
        
        # 跳过无效数据
        if not iterations or not block_counts:
            logger.warning("染色体 %s 没有有效数据，跳过", chr_name)
            continue
        
        # 计算统计信息
        initial_count = block_counts[0]
        final_count = block_counts[-1]
        
        # 计算变化率，避免除零错误
        if initial_count > 0:
            change_rate = ((final_count - initial_count) / initial_count) * 100
        else:
            change_rate = 0.0
        
        # 存储染色体数据
        chr_data[chr_name] = {
            "iterations": iterations,
            "block_counts": block_counts,
            "sample_names": sample_names,
            "final_count": final_count,
            "initial_count": initial_count,
            "change_rate": change_rate
        }
        
        logger.debug("成功收集染色体 %s 数据: %d 轮迭代, 最终block数=%d", 
                    chr_name, len(iterations), final_count)
    
    logger.info("共收集 %d 条染色体数据 (总共 %d 条)", len(chr_data), len(chr_list))
    return chr_data


def generate_html_report(chr_data: Dict[str, Dict], output_path: str) -> bool:
    """
    生成基因组级别的HTML交互式报告
    
    生成一个自包含的HTML文件，包含：
    - 上部：柱状图显示各染色体最终block数量
    - 中部：可点击的染色体指示器
    - 下部：数据表格显示选中染色体的详细迭代数据
    
    Args:
        chr_data: 所有染色体的数据字典，由collect_all_chr_data返回
        output_path: HTML输出路径
    
    Returns:
        成功返回True，失败返回False
    """
    if not chr_data:
        logger.warning("没有染色体数据可生成报告")
        return False
    
    # 确保输出目录存在
    out_dir = os.path.dirname(output_path)
    if out_dir:
        try:
            os.makedirs(out_dir, exist_ok=True)
        except OSError as e:
            logger.error("创建输出目录失败: %s, 错误: %s", out_dir, str(e))
            return False
    
    # 准备染色体数据为JavaScript格式
    import json
    
    # 转换数据为JSON格式
    js_chr_data = {}
    for chr_name, data in chr_data.items():
        js_chr_data[chr_name] = {
            "iterations": data["iterations"],
            "block_counts": data["block_counts"],
            "sample_names": data["sample_names"],
            "final_count": data["final_count"],
            "initial_count": data["initial_count"],
            "change_rate": round(data["change_rate"], 2)
        }
    
    chr_data_json = json.dumps(js_chr_data, ensure_ascii=False, indent=2)
    
    # 获取排序后的染色体列表（用于柱状图和指示器）
    chr_names = sorted(chr_data.keys())
    final_counts = [chr_data[c]["final_count"] for c in chr_names]
    
    # 生成HTML内容
    html_content = f'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Core Block Trend Report - Genome Level</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background-color: #f5f7fa;
            color: #333;
            line-height: 1.6;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        
        h1 {{
            text-align: center;
            color: #2c3e50;
            margin-bottom: 30px;
            font-size: 28px;
        }}
        
        .section {{
            background: white;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            padding: 20px;
        }}
        
        .section-title {{
            font-size: 18px;
            color: #34495e;
            margin-bottom: 15px;
            padding-bottom: 10px;
            border-bottom: 2px solid #3498db;
        }}
        
        /* Bar Chart Section */
        .bar-chart-container {{
            width: 100%;
            overflow-x: auto;
        }}
        
        .bar-chart {{
            display: flex;
            align-items: flex-end;
            justify-content: center;
            height: 300px;
            padding: 20px 10px;
            gap: 8px;
            min-width: fit-content;
        }}
        
        .bar-wrapper {{
            display: flex;
            flex-direction: column;
            align-items: center;
            min-width: 40px;
        }}
        
        .bar {{
            width: 35px;
            background: linear-gradient(180deg, #3498db 0%, #2980b9 100%);
            border-radius: 4px 4px 0 0;
            transition: all 0.3s ease;
            cursor: pointer;
            position: relative;
        }}
        
        .bar:hover {{
            background: linear-gradient(180deg, #e74c3c 0%, #c0392b 100%);
            transform: scaleX(1.1);
        }}
        
        .bar-value {{
            font-size: 11px;
            color: #666;
            margin-bottom: 5px;
            white-space: nowrap;
        }}
        
        .bar-label {{
            font-size: 11px;
            color: #666;
            margin-top: 8px;
            transform: rotate(-45deg);
            white-space: nowrap;
        }}
        
        /* Chromosome Indicators Section */
        .chr-indicators {{
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            justify-content: center;
            padding: 10px;
        }}
        
        .chr-indicator {{
            padding: 10px 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 8px;
            cursor: pointer;
            transition: all 0.3s ease;
            font-weight: 500;
            box-shadow: 0 2px 5px rgba(0,0,0,0.2);
        }}
        
        .chr-indicator:hover {{
            transform: translateY(-3px);
            box-shadow: 0 4px 15px rgba(0,0,0,0.3);
        }}
        
        .chr-indicator.active {{
            background: linear-gradient(135deg, #e74c3c 0%, #c0392b 100%);
            transform: scale(1.05);
        }}
        
        .chr-indicator .chr-name {{
            font-size: 14px;
        }}
        
        .chr-indicator .chr-count {{
            font-size: 12px;
            opacity: 0.9;
        }}
        
        /* Data Table Section */
        .table-container {{
            overflow-x: auto;
        }}
        
        .placeholder {{
            text-align: center;
            padding: 40px;
            color: #95a5a6;
            font-size: 16px;
        }}
        
        .data-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
        }}
        
        .data-table th {{
            background: #34495e;
            color: white;
            padding: 12px 15px;
            text-align: left;
            font-weight: 600;
        }}
        
        .data-table td {{
            padding: 10px 15px;
            border-bottom: 1px solid #ecf0f1;
        }}
        
        .data-table tr:hover {{
            background-color: #f8f9fa;
        }}
        
        .data-table tr:nth-child(even) {{
            background-color: #fafbfc;
        }}
        
        .selected-chr-title {{
            font-size: 16px;
            color: #2c3e50;
            margin-bottom: 10px;
        }}
        
        .stats-row {{
            display: flex;
            gap: 20px;
            margin-bottom: 15px;
            flex-wrap: wrap;
        }}
        
        .stat-item {{
            background: #ecf0f1;
            padding: 10px 20px;
            border-radius: 6px;
        }}
        
        .stat-label {{
            font-size: 12px;
            color: #7f8c8d;
        }}
        
        .stat-value {{
            font-size: 18px;
            font-weight: bold;
            color: #2c3e50;
        }}
        
        .stat-value.positive {{
            color: #27ae60;
        }}
        
        .stat-value.negative {{
            color: #e74c3c;
        }}
        
        /* Footer */
        .footer {{
            text-align: center;
            padding: 20px;
            color: #95a5a6;
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Core Block Trend Report</h1>
        
        <!-- Bar Chart Section -->
        <div class="section">
            <div class="section-title">📊 Final Block Count by Chromosome</div>
            <div class="bar-chart-container">
                <div class="bar-chart" id="barChart"></div>
            </div>
        </div>
        
        <!-- Chromosome Indicators Section -->
        <div class="section">
            <div class="section-title">🔬 Select Chromosome to View Details</div>
            <div class="chr-indicators" id="chrIndicators"></div>
        </div>
        
        <!-- Trend Chart Section (for selected chromosome) -->
        <div class="section">
            <div class="section-title">📈 Iteration Trend Chart</div>
            <div id="trendChartContainer">
                <div class="placeholder">👆 Click on a chromosome above to view the trend chart</div>
            </div>
        </div>
        
        <!-- Data Table Section -->
        <div class="section">
            <div class="section-title">📋 Iteration Data</div>
            <div id="dataTableContainer">
                <div class="placeholder">👆 Click on a chromosome above to view detailed iteration data</div>
            </div>
        </div>
        
        <div class="footer">
            Generated by panCB Pipeline - Core Block Trend Reporter
        </div>
    </div>
    
    <script>
        // Embedded chromosome data
        const chromosomeData = {chr_data_json};
        
        // Get sorted chromosome names
        const chrNames = {json.dumps(chr_names)};
        
        // Calculate max value for bar chart scaling
        const maxCount = Math.max(...chrNames.map(c => chromosomeData[c].final_count));
        
        // Initialize bar chart
        function initBarChart() {{
            const container = document.getElementById('barChart');
            container.innerHTML = '';
            
            chrNames.forEach(chrName => {{
                const data = chromosomeData[chrName];
                const barHeight = (data.final_count / maxCount) * 250;
                
                const wrapper = document.createElement('div');
                wrapper.className = 'bar-wrapper';
                
                const value = document.createElement('div');
                value.className = 'bar-value';
                value.textContent = data.final_count.toLocaleString();
                
                const bar = document.createElement('div');
                bar.className = 'bar';
                bar.style.height = barHeight + 'px';
                bar.onclick = () => selectChromosome(chrName);
                bar.title = chrName + ': ' + data.final_count.toLocaleString() + ' blocks';
                
                const label = document.createElement('div');
                label.className = 'bar-label';
                label.textContent = chrName;
                
                wrapper.appendChild(value);
                wrapper.appendChild(bar);
                wrapper.appendChild(label);
                container.appendChild(wrapper);
            }});
        }}
        
        // Initialize chromosome indicators
        function initChrIndicators() {{
            const container = document.getElementById('chrIndicators');
            container.innerHTML = '';
            
            chrNames.forEach(chrName => {{
                const data = chromosomeData[chrName];
                
                const indicator = document.createElement('div');
                indicator.className = 'chr-indicator';
                indicator.id = 'indicator-' + chrName;
                indicator.onclick = () => selectChromosome(chrName);
                
                const nameDiv = document.createElement('div');
                nameDiv.className = 'chr-name';
                nameDiv.textContent = chrName;
                
                const countDiv = document.createElement('div');
                countDiv.className = 'chr-count';
                countDiv.textContent = data.final_count.toLocaleString() + ' blocks';
                
                indicator.appendChild(nameDiv);
                indicator.appendChild(countDiv);
                container.appendChild(indicator);
            }});
        }}
        
        // Chart instance for trend chart
        let trendChart = null;
        
        // Select chromosome and display data
        function selectChromosome(chrName) {{
            // Update active indicator
            document.querySelectorAll('.chr-indicator').forEach(el => {{
                el.classList.remove('active');
            }});
            document.getElementById('indicator-' + chrName).classList.add('active');
            
            // Get chromosome data
            const data = chromosomeData[chrName];
            
            // Build trend chart
            buildTrendChart(chrName, data);
            
            // Build stats row
            const changeClass = data.change_rate >= 0 ? 'positive' : 'negative';
            const changeSign = data.change_rate >= 0 ? '+' : '';
            
            let html = '<div class="selected-chr-title">📍 ' + chrName + ' - Iteration Details</div>';
            html += '<div class="stats-row">';
            html += '<div class="stat-item"><div class="stat-label">Initial Count</div><div class="stat-value">' + data.initial_count.toLocaleString() + '</div></div>';
            html += '<div class="stat-item"><div class="stat-label">Final Count</div><div class="stat-value">' + data.final_count.toLocaleString() + '</div></div>';
            html += '<div class="stat-item"><div class="stat-label">Change Rate</div><div class="stat-value ' + changeClass + '">' + changeSign + data.change_rate.toFixed(2) + '%</div></div>';
            html += '<div class="stat-item"><div class="stat-label">Iterations</div><div class="stat-value">' + data.iterations.length + '</div></div>';
            html += '</div>';
            
            // Build data table
            html += '<table class="data-table">';
            html += '<thead><tr><th>iteration</th><th>block_count</th><th>sample_name</th></tr></thead>';
            html += '<tbody>';
            
            for (let i = 0; i < data.iterations.length; i++) {{
                html += '<tr>';
                html += '<td>' + data.iterations[i] + '</td>';
                html += '<td>' + data.block_counts[i].toLocaleString() + '</td>';
                html += '<td>' + (data.sample_names[i] || '-') + '</td>';
                html += '</tr>';
            }}
            
            html += '</tbody></table>';
            
            document.getElementById('dataTableContainer').innerHTML = html;
        }}
        
        // Build trend chart for selected chromosome
        function buildTrendChart(chrName, data) {{
            const container = document.getElementById('trendChartContainer');
            
            // Create canvas for chart
            container.innerHTML = '<div class="selected-chr-title">📈 ' + chrName + ' - Core Block Trend</div>' +
                '<div style="height: 350px; position: relative;"><canvas id="trendCanvas"></canvas></div>';
            
            const ctx = document.getElementById('trendCanvas').getContext('2d');
            
            // Destroy previous chart if exists
            if (trendChart) {{
                trendChart.destroy();
            }}
            
            // Prepare labels (show every Nth label if too many)
            const labels = data.iterations;
            const step = labels.length > 30 ? Math.ceil(labels.length / 20) : 1;
            
            // Create new chart
            trendChart = new Chart(ctx, {{
                type: 'line',
                data: {{
                    labels: labels,
                    datasets: [{{
                        label: 'Core Block Count',
                        data: data.block_counts,
                        borderColor: '#2E86AB',
                        backgroundColor: 'rgba(46, 134, 171, 0.2)',
                        borderWidth: 2,
                        pointRadius: labels.length > 50 ? 1 : 3,
                        pointHoverRadius: 5,
                        fill: true,
                        tension: 0.1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            display: false
                        }},
                        tooltip: {{
                            callbacks: {{
                                title: function(context) {{
                                    const idx = context[0].dataIndex;
                                    return data.iterations[idx] + (data.sample_names[idx] ? ' (' + data.sample_names[idx] + ')' : '');
                                }},
                                label: function(context) {{
                                    return 'Block Count: ' + context.raw.toLocaleString();
                                }}
                            }}
                        }}
                    }},
                    scales: {{
                        x: {{
                            title: {{
                                display: true,
                                text: 'Iteration Round'
                            }},
                            ticks: {{
                                maxRotation: 45,
                                minRotation: 45,
                                callback: function(value, index) {{
                                    return index % step === 0 ? this.getLabelForValue(value) : '';
                                }}
                            }}
                        }},
                        y: {{
                            title: {{
                                display: true,
                                text: 'Core Block Count'
                            }},
                            beginAtZero: false
                        }}
                    }}
                }}
            }});
        }}
        
        // Initialize on page load
        document.addEventListener('DOMContentLoaded', function() {{
            initBarChart();
            initChrIndicators();
        }});
    </script>
</body>
</html>'''
    
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        # 验证文件是否成功创建
        if os.path.isfile(output_path) and os.path.getsize(output_path) > 0:
            logger.info("HTML报告已生成: %s", output_path)
            return True
        else:
            logger.error("HTML报告文件创建失败或为空: %s", output_path)
            return False
            
    except IOError as e:
        logger.error("写入HTML报告失败: %s, 错误: %s", output_path, str(e))
        return False
    except Exception as e:
        logger.error("生成HTML报告时发生未知错误: %s", str(e))
        return False


def generate_genome_report(output_dir: str, chr_list: List[str]) -> bool:
    """
    主入口函数：收集数据并生成HTML报告
    
    协调数据收集和HTML生成过程，处理错误并返回适当的状态。
    
    Args:
        output_dir: 输出根目录（如 'output'）
        chr_list: 染色体名称列表（如 ['Chr01', 'Chr02', ...]）
    
    Returns:
        成功返回True，失败返回False
    
    Example:
        >>> success = generate_genome_report('output', ['Chr01', 'Chr02', 'Chr03'])
        >>> print(success)
        True
    """
    logger.info("开始生成基因组级别HTML报告...")
    logger.info("输出目录: %s", output_dir)
    logger.info("染色体列表: %s", chr_list)
    
    # 验证输入参数
    if not output_dir:
        logger.error("输出目录不能为空")
        return False
    
    if not chr_list:
        logger.error("染色体列表不能为空")
        return False
    
    # 步骤1: 收集所有染色体数据
    logger.info("步骤1: 收集染色体数据...")
    try:
        chr_data = collect_all_chr_data(output_dir, chr_list)
    except Exception as e:
        logger.error("收集染色体数据时发生错误: %s", str(e))
        return False
    
    if not chr_data:
        logger.warning("没有收集到任何有效的染色体数据")
        return False
    
    logger.info("成功收集 %d 条染色体数据", len(chr_data))
    
    # 步骤2: 构建HTML报告输出路径
    # 输出路径: output/pan_core_block/report/core_block_trend_report.html
    report_dir = os.path.join(output_dir, "pan_core_block", "report")
    report_path = os.path.join(report_dir, "core_block_trend_report.html")
    
    logger.info("步骤2: 生成HTML报告...")
    logger.info("报告输出路径: %s", report_path)
    
    # 步骤3: 生成HTML报告
    try:
        success = generate_html_report(chr_data, report_path)
    except Exception as e:
        logger.error("生成HTML报告时发生错误: %s", str(e))
        return False
    
    if success:
        logger.info("基因组级别HTML报告生成完成!")
        logger.info("报告文件: %s", report_path)
    else:
        logger.error("HTML报告生成失败")
    
    return success


def main():
    parser = argparse.ArgumentParser(
        description='绘制核心block迭代变化趋势折线图',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
    python Rn_block_trend.py -i Chr01_trend.tsv -o Chr01_trend.png
    python Rn_block_trend.py -i Chr01_trend.tsv -o Chr01_trend.png -t "Chr01 Core Block Trend"
        '''
    )
    parser.add_argument('-i', '--input', required=True, help='输入TSV文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出图片路径（PNG格式）')
    parser.add_argument('-t', '--title', default=None, help='图表标题（可选）')
    parser.add_argument('-c', '--chr', default=None, help='染色体名称（可选，用于自动生成标题）')
    parser.add_argument('-v', '--verbose', action='store_true', help='显示详细日志')
    
    args = parser.parse_args()
    
    # 配置日志
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    if not os.path.isfile(args.input):
        logger.error("输入文件不存在: %s", args.input)
        sys.exit(1)
    
    # 读取数据
    iterations, block_counts, sample_names = read_tsv(args.input)
    
    if not iterations:
        logger.error("TSV文件为空或格式不正确: %s", args.input)
        sys.exit(1)
    
    # 绘制图表
    success = plot_trend(
        iterations, block_counts, sample_names,
        args.output, args.title, args.chr
    )
    
    if not success:
        sys.exit(1)


if __name__ == '__main__':
    main()
