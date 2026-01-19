import json
import os
import subprocess
import logging
from datetime import datetime
from typing import Any, Dict, List, Optional
from Bio import SeqIO


def ensure_dir(path: str) -> None:
    """
    确保目录存在，如果不存在则创建。
    
    Args:
        path: 目录路径
    """
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


def read_accs(accs_file: str) -> List[str]:
    """
    从样本列表文件中读取样本名称。
    
    Args:
        accs_file: 样本列表文件路径（每行一个样本名称，第一列为名称）
    
    Returns:
        样本名称列表
    """
    accs: List[str] = []
    with open(accs_file) as f:
        for line in f:
            name = line.strip().split()[0]
            if name:
                accs.append(name)
    return accs


def relation(interval1, interval2):
    """
    判断两个区间的空间关系。
    
    Args:
        interval1: 区间1 [start, end]
        interval2: 区间2 [start, end]
    
    Returns:
        0: 完全相同
        1: 不相交
        2: 部分重叠
        3: 一个包含另一个
    """
    min1, max1 = sorted(interval1)
    min2, max2 = sorted(interval2)
    if min1 == min2 and max1 == max2:
        return 0
    if max1 < min2 or max2 < min1:
        return 1
    if (min1 < min2 <= max1 < max2) or (min2 < min1 <= max2 < max1):
        return 2
    if (min1 <= min2 <= max2 <= max1) or (min2 <= min1 <= max1 <= max2):
        return 3


def intersection(interval1, interval2):
    """
    计算两个区间的交集。
    
    Args:
        interval1: 区间1 [start, end]
        interval2: 区间2 [start, end]
    
    Returns:
        交集区间 [start, end]，如果无交集返回空列表
    """
    nums = sorted(interval1 + interval2)
    if relation(interval1, interval2) != 1:
        return [nums[1], nums[2]]
    else:
        return []


def merge(intervals):
    """
    合并重叠的区间列表。
    
    Args:
        intervals: 区间列表，每个区间为[start, end]
    
    Returns:
        合并后的区间列表
    """
    intervals.sort(key=lambda x: x[0])
    merged = []
    for interval in intervals:
        if not merged or merged[-1][-1] < interval[0]:
            merged.append(interval)
        else:
            merged[-1][-1] = max(merged[-1][-1], interval[-1])
    return merged


def fasta_has_bases(fasta_path: str) -> bool:
    """
    检查FASTA文件是否包含非空序列。
    
    Args:
        fasta_path: FASTA文件路径
    
    Returns:
        如果文件包含至少一条非空序列，返回True；否则返回False
    """
    try:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            if len(rec.seq) > 0:
                return True
        return False
    except Exception:
        return False


def write_minimal_placeholder_fasta(fasta_path: str, seq_id: str) -> None:
    """
    写入最小占位FASTA文件（1bp的N序列）。
    用于处理缺失染色体的占位。
    
    Args:
        fasta_path: 输出FASTA文件路径
        seq_id: 序列ID
    """
    with open(fasta_path, "w") as fo:
        fo.write(f">{seq_id}\nN\n")


def run_cmd(
    cmd: str,
    logger: logging.Logger,
    critical: bool = False,
    allowed_exit_codes: tuple = (0,),
) -> int:
    """
    执行系统命令
    
    Args:
        cmd: 要执行的命令
        logger: 日志记录器
        critical: 是否为关键命令（失败时记录为ERROR）
        allowed_exit_codes: 允许的退出码集合（默认只允许0）
    
    Returns:
        命令的退出码
    """
    ret = subprocess.call(cmd, shell=True)
    # 仅当返回码不在允许集合内时，才输出 WARNING/ERROR，
    # 例如：某些外部工具（如 TRF）在“无结果/短序列”时可能返回 1，但并非真正失败。
    if ret not in allowed_exit_codes:
        if critical:
            logger.error("命令失败(ret=%s): %s", ret, cmd)
        else:
            logger.warning("命令失败(ret=%s): %s", ret, cmd)
    return ret


def ensure_fai_exists(fasta_path: str, logger: logging.Logger = None) -> None:
    """
    确保FASTA文件的索引文件(.fai)存在，如果不存在则生成。
    
    Args:
        fasta_path: FASTA文件路径
        logger: 可选的日志记录器
    """
    fai = fasta_path + ".fai"
    if not os.path.isfile(fai):
        if logger:
            logger.info("samtools faidx 生成索引: %s", fasta_path)
        subprocess.call(f"samtools faidx {fasta_path}", shell=True)


def count_lines(path: str) -> int:
    """
    统计文件行数。
    
    Args:
        path: 文件路径
    
    Returns:
        文件行数，如果文件不存在或读取失败返回0
    """
    try:
        n = 0
        with open(path, 'r') as f:
            for _ in f:
                n += 1
        return n
    except Exception:
        return 0


def file_size(path: str) -> int:
    """
    获取文件大小（字节）。
    
    Args:
        path: 文件路径
    
    Returns:
        文件大小（字节），如果文件不存在或读取失败返回0
    """
    try:
        return os.path.getsize(path)
    except Exception:
        return 0


class ResumeManager:
    """
    维护流程的断点续跑状态，记录每个步骤的执行情况。
    """

    STEP_SEQUENCE = [
        ("split", "步骤1：按染色体拆分"),
        ("align_core", "步骤2-4：核心比对与迭代"),
        ("seq_extract", "步骤5：提取核心序列"),
        ("alternative", "步骤5.5：alternative block"),
        ("report_summary", "步骤6：核心区块统计"),
        ("report_detail", "步骤7：详细报告"),
    ]

    def __init__(self, output_dir: str) -> None:
        self.output_dir = output_dir
        self.checkpoint_dir = os.path.join(output_dir, "checkpoints")
        ensure_dir(self.checkpoint_dir)
        self.state_path = os.path.join(self.checkpoint_dir, "resume_state.json")
        self.state: Dict[str, Any] = {"steps": {}}
        self._load()

    def _load(self) -> None:
        if os.path.isfile(self.state_path):
            try:
                with open(self.state_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                    if isinstance(data, dict):
                        self.state = data
            except Exception:
                logging.getLogger("panCB.resume").warning("无法读取历史断点信息，忽略旧状态: %s", self.state_path)

    def reset(self) -> None:
        """新运行时清空状态文件。"""
        self.state = {"steps": {}}
        if os.path.isfile(self.state_path):
            try:
                os.remove(self.state_path)
            except Exception:
                logging.getLogger("panCB.resume").warning("清理旧断点信息失败: %s", self.state_path)

    def _save(self) -> None:
        try:
            with open(self.state_path, "w", encoding="utf-8") as f:
                json.dump(self.state, f, ensure_ascii=False, indent=2)
        except Exception:
            logging.getLogger("panCB.resume").warning("写入断点信息失败: %s", self.state_path)

    def mark_step_started(self, step: str, meta: Optional[Dict[str, Any]] = None) -> None:
        entry = self.state.setdefault("steps", {}).setdefault(step, {})
        entry["status"] = "running"
        entry["updated_at"] = datetime.utcnow().isoformat()
        if meta is not None:
            entry["meta"] = meta
        self.state["last_started"] = step
        if self.state.get("last_completed") == step:
            self.state.pop("last_completed", None)
        self._save()

    def mark_step_done(self, step: str, meta: Optional[Dict[str, Any]] = None) -> None:
        entry = self.state.setdefault("steps", {}).setdefault(step, {})
        entry["status"] = "done"
        entry["updated_at"] = datetime.utcnow().isoformat()
        if meta is not None:
            entry["meta"] = meta
        self.state["last_completed"] = step
        self._save()

    def mark_step_skipped(self, step: str, reason: str) -> None:
        entry = self.state.setdefault("steps", {}).setdefault(step, {})
        entry["status"] = "skipped"
        entry["reason"] = reason
        entry["updated_at"] = datetime.utcnow().isoformat()
        self.state["last_completed"] = step
        self._save()

    def is_step_done(self, step: str) -> bool:
        status = self.state.get("steps", {}).get(step, {}).get("status")
        return status == "done" or status == "skipped"

    def get_step_status(self, step: str) -> Optional[str]:
        status = self.state.get("steps", {}).get(step, {}).get("status")
        return status if isinstance(status, str) else None

    def get_step_meta(self, step: str) -> Dict[str, Any]:
        entry = self.state.get("steps", {}).get(step, {})
        meta = entry.get("meta")
        if isinstance(meta, dict):
            return meta
        return {}

    def update_step_meta(self, step: str, meta: Dict[str, Any]) -> None:
        if not isinstance(meta, dict):
            return
        entry = self.state.setdefault("steps", {}).setdefault(step, {})
        current_meta = entry.get("meta")
        if not isinstance(current_meta, dict):
            current_meta = {}
        current_meta.update(meta)
        entry["meta"] = current_meta
        entry["updated_at"] = datetime.utcnow().isoformat()
        self._save()

    def reopen_step(self, step: str, meta: Optional[Dict[str, Any]] = None) -> None:
        """
        将指定步骤重新标记为 running。用于检测到产物缺失时自动回退。
        """
        entry = self.state.setdefault("steps", {}).setdefault(step, {})
        entry["status"] = "running"
        entry["updated_at"] = datetime.utcnow().isoformat()
        if meta is not None:
            entry["meta"] = meta
        if self.state.get("last_completed") == step:
            self.state.pop("last_completed", None)
        self.state["last_started"] = step
        self._save()

    def get_last_completed_step(self) -> Optional[str]:
        step = self.state.get("last_completed")
        return step if isinstance(step, str) else None

    def get_next_pending_step(self) -> Optional[str]:
        steps = self.state.get("steps", {})
        for step, _ in self.STEP_SEQUENCE:
            status = steps.get(step, {}).get("status")
            if status != "done" and status != "skipped":
                return step
        return None

    def describe_progress(self) -> str:
        last_completed = self.get_last_completed_step() or "无"
        next_pending = self.get_next_pending_step() or "全部完成"
        label_map = {k: v for k, v in self.STEP_SEQUENCE}
        last_label = label_map.get(last_completed, last_completed)
        next_label = label_map.get(next_pending, next_pending)
        return f"上次已完成：{last_label}；下一步：{next_label}"