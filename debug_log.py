"""
debug_log.py - 调试日志模块

提供统一的调试日志管理功能，用于输出详细的 block 过滤日志。
当启用 --debug 参数时，会将详细的过滤过程记录到 output/logs/ChrXX_debug.log 文件中。
"""

import os
from typing import Optional
import utils


class DebugLogger:
    """
    调试日志管理器
    
    用于管理染色体级别的调试日志文件，提供统一的日志写入接口。
    """
    
    def __init__(self, output_dir: str, chr_name: str, enabled: bool = False):
        """
        初始化调试日志管理器
        
        Args:
            output_dir: 输出根目录
            chr_name: 染色体名称
            enabled: 是否启用调试模式
        """
        self.output_dir = output_dir
        self.chr_name = chr_name
        self.enabled = enabled
        self._file = None
        self._log_path = None
        
        if enabled:
            logs_dir = os.path.join(output_dir, "logs")
            utils.ensure_dir(logs_dir)
            self._log_path = os.path.join(logs_dir, f"{chr_name}_debug.log")
    
    def __enter__(self):
        """上下文管理器入口"""
        if self.enabled and self._log_path:
            self._file = open(self._log_path, 'a', encoding='utf-8')
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """上下文管理器退出"""
        if self._file:
            self._file.close()
            self._file = None
        return False
    
    @property
    def file(self) -> Optional[object]:
        """获取文件句柄（用于兼容现有代码）"""
        return self._file
    
    def write(self, message: str) -> None:
        """
        写入调试日志
        
        Args:
            message: 日志消息
        """
        if self._file:
            self._file.write(message)
    
    def write_section_header(self, phase: str, name: str) -> None:
        """
        写入阶段标题
        
        Args:
            phase: 阶段名称（如"阶段1"）
            name: 处理对象名称
        """
        if self._file:
            self._file.write(f"\n{'='*80}\n")
            self._file.write(f"【{phase}】{name}\n")
            self._file.write(f"{'='*80}\n")
    
    def write_section_footer(self) -> None:
        """写入阶段结束分隔符"""
        if self._file:
            self._file.write(f"{'='*80}\n\n")
    
    def write_subsection(self, title: str) -> None:
        """
        写入子节标题
        
        Args:
            title: 子节标题
        """
        if self._file:
            self._file.write(f"--- {title} ---\n")
    
    def write_param(self, name: str, value) -> None:
        """
        写入参数信息
        
        Args:
            name: 参数名称
            value: 参数值
        """
        if self._file:
            self._file.write(f"{name}: {value}\n")
    
    def write_filtered_block(self, block_num: int, reason: str, details: str) -> None:
        """
        写入被过滤的 block 信息
        
        Args:
            block_num: block 编号
            reason: 过滤原因
            details: 详细信息
        """
        if self._file:
            self._file.write(f"Block #{block_num} 被过滤 [{reason}]: {details}\n")
    
    def write_overlap_info(self, overlap_num: int, chr_name: str, overlap_len: int, 
                           block1_info: str, block2_info: str) -> None:
        """
        写入重叠处理信息
        
        Args:
            overlap_num: 重叠编号
            chr_name: 染色体名称
            overlap_len: 重叠长度
            block1_info: block1 信息
            block2_info: block2 信息
        """
        if self._file:
            self._file.write(f"  重叠 #{overlap_num}: chr={chr_name}, overlap_len={overlap_len}bp, "
                           f"{block1_info}, {block2_info}\n")
    
    def write_stats(self, stats: dict) -> None:
        """
        写入统计信息
        
        Args:
            stats: 统计字典
        """
        if self._file:
            for key, value in stats.items():
                self._file.write(f"  - {key}: {value}\n")
    
    def write_summary(self, title: str, items: dict) -> None:
        """
        写入汇总信息
        
        Args:
            title: 汇总标题
            items: 汇总项目字典
        """
        if self._file:
            self._file.write(f"\n--- {title} ---\n")
            for key, value in items.items():
                self._file.write(f"{key}: {value}\n")


def create_debug_logger(output_dir: str, chr_name: str, debug: bool = False) -> DebugLogger:
    """
    创建调试日志管理器的工厂函数
    
    Args:
        output_dir: 输出根目录
        chr_name: 染色体名称
        debug: 是否启用调试模式
    
    Returns:
        DebugLogger 实例
    """
    return DebugLogger(output_dir, chr_name, debug)
