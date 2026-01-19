"""
core_blocks.py - 向后兼容层

本文件作为向后兼容层，从 core_blocks 包中重新导出所有函数。
现有代码可以继续使用 `from core_blocks import xxx` 或 `import core_blocks` 的方式导入。

实际实现已拆分到 core_blocks/ 包中的各个模块：
- core_blocks/coord_utils.py: 坐标计算工具函数
- core_blocks/block_filters.py: Block后处理过滤逻辑
- core_blocks/fast_iteration.py: 快速迭代模式实现
- core_blocks/standard_iteration.py: 标准迭代模式实现
- core_blocks/iterate_control.py: 迭代控制与进度管理
"""

# 从 core_blocks 包导入所有公共接口
from core_blocks import (
    # 坐标工具函数
    calculate_coordinate_update,
    calculate_flanking_coords,
    nt_stat,
    
    # Block过滤函数
    clean_transposition,
    process_overlaps,
    filter_by_length,
    apply_all_filters,
    FilterStats,
    
    # 快速迭代函数
    extract_core_blocks_fast,
    validate_final_coordinates,
    
    # 标准迭代函数
    extract_core_blocks_once,
    
    # 迭代控制函数
    run_iterative_core_blocks_for_chr,
    generate_alignment_priority_for_chr,
    _parse_query_name,
    _sum_core_length_by_chr,
    
    # 进度检查函数
    get_expected_iter_count,
    get_last_completed_iter,
    has_final_iter_result,
    is_iteration_complete,
    
    # 完整性检查函数
    has_chr_level_coords,
    has_chr_trend_files,
    check_chr_completeness,
    
    # 重生成函数
    regenerate_chr_trend_from_iter_files,
    regenerate_chr_trend_png_from_tsv,
    ensure_chr_outputs_complete,
    
    # 报告生成函数
    generate_genome_trend_report,
)

# 保持原有的导入兼容性
import os
import re
import copy
import logging
from typing import List, Dict, Tuple, Optional
from Bio import SeqIO
import utils
import debug_log
from Rn_block_trend import generate_genome_report

__all__ = [
    # 坐标工具
    'calculate_coordinate_update',
    'calculate_flanking_coords',
    'nt_stat',
    # Block过滤
    'clean_transposition',
    'process_overlaps',
    'filter_by_length',
    'apply_all_filters',
    'FilterStats',
    # 迭代函数
    'extract_core_blocks_fast',
    'validate_final_coordinates',
    'extract_core_blocks_once',
    # 迭代控制
    'run_iterative_core_blocks_for_chr',
    'generate_alignment_priority_for_chr',
    '_parse_query_name',
    '_sum_core_length_by_chr',
    # 进度检查
    'get_expected_iter_count',
    'get_last_completed_iter',
    'has_final_iter_result',
    'is_iteration_complete',
    # 完整性检查
    'has_chr_level_coords',
    'has_chr_trend_files',
    'check_chr_completeness',
    # 重生成
    'regenerate_chr_trend_from_iter_files',
    'regenerate_chr_trend_png_from_tsv',
    'ensure_chr_outputs_complete',
    # 报告生成
    'generate_genome_trend_report',
]
