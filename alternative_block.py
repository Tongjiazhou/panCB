import os
import logging
from typing import List, Tuple
from Bio import SeqIO

import utils


def _read_chr_length_from_fai(fai_path: str, chr_name: str) -> int:
    """
    从FASTA索引文件(.fai)中读取指定染色体的长度。
    
    Args:
        fai_path: FASTA索引文件路径
        chr_name: 染色体名称
    
    Returns:
        染色体长度，如果未找到或读取失败返回0
    """
    try:
        with open(fai_path) as f:
            for line in f:
                cols = line.strip().split()
                if cols and cols[0] == chr_name:
                    return int(cols[1])
    except Exception:
        return 0
    return 0


def _read_core_intervals(core_bed_path: str) -> List[Tuple[int, int]]:
    """
    从BED文件中读取核心区间列表。
    
    Args:
        core_bed_path: 核心block的BED文件路径（3列：chr start end）
    
    Returns:
        排序后的区间列表，每个区间为(start, end)元组
    """
    intervals: List[Tuple[int, int]] = []
    with open(core_bed_path) as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) >= 3:
                try:
                    s = int(cols[1])
                    e = int(cols[2])
                    intervals.append((s, e))
                except ValueError:
                    pass
    intervals.sort(key=lambda x: x[0])
    return intervals


def _build_core_and_alternative(chr_len: int, core_intervals: List[Tuple[int, int]]) -> List[Tuple[int, int, str]]:
    """
    基于核心区间构建核心和alternative区间的完整列表。
    
    Args:
        chr_len: 染色体总长度
        core_intervals: 核心区间列表（已排序）
    
    Returns:
        包含所有区间的列表，每个元素为(start, end, label)元组，label为"core"或"alternative"
    """
    merged: List[Tuple[int, int, str]] = []
    if chr_len <= 0:
        return merged

    if not core_intervals:
        merged.append((0, chr_len, "alternative"))
        return merged

    # Leading alternative
    first_s = max(0, core_intervals[0][0])
    if first_s > 0:
        merged.append((0, first_s, "alternative"))

    # Core and gaps
    for i, (s, e) in enumerate(core_intervals):
        if e > s:
            merged.append((s, e, "core"))
        if i + 1 < len(core_intervals):
            ns, _ = core_intervals[i + 1]
            if ns > e:
                merged.append((e, ns, "alternative"))

    # Trailing alternative
    last_e = core_intervals[-1][1]
    if chr_len > last_e:
        merged.append((last_e, chr_len, "alternative"))

    # Sort by start
    merged.sort(key=lambda x: x[0])
    return merged


def _write_sorted_bed_and_len(records: List[Tuple[int, int, str]], chr_name: str, sorted_bed: str, len_path: str) -> None:
    """
    将核心和alternative区间写入排序后的BED文件和长度文件。
    
    Args:
        records: 区间记录列表，每个元素为(start, end, label)
        chr_name: 染色体名称
        sorted_bed: 输出BED文件路径
        len_path: 输出长度文件路径
    """
    with open(sorted_bed, "w") as fo_bed, open(len_path, "w") as fo_len:
        idx = 1
        for s, e, label in records:
            fo_bed.write(f"{chr_name}\t{s}\t{e}\t{label}\n")
            length_value = (e - s) if label == "core" else (s - e)  # mimic original: alternative uses negative length
            fo_len.write(f"{idx}\t{length_value}\t{label}\n")
            idx += 1


def _extract_sorted_fasta(genome_fa: str, sorted_bed: str, out_fasta: str, logger: logging.Logger) -> None:
    """
    从排序后的BED文件中提取序列（包含core和alternative）。
    
    Args:
        genome_fa: 基因组FASTA文件路径
        sorted_bed: 排序后的BED文件路径
        out_fasta: 输出FASTA文件路径
        logger: 日志记录器
    """
    cmd = f"bedtools getfasta -name -fi {genome_fa} -bed {sorted_bed} -fo {out_fasta}"
    utils.run_cmd(cmd, logger)

def _write_alternative_bed(records: List[Tuple[int, int, str]], chr_name: str, alt_bed: str) -> None:
    """
    将alternative区间写入BED文件。
    
    Args:
        records: 区间记录列表
        chr_name: 染色体名称
        alt_bed: 输出alternative BED文件路径
    """
    with open(alt_bed, "w") as fo:
        for s, e, label in records:
            if label == "alternative" and e > s:
                fo.write(f"{chr_name}\t{s}\t{e}\talternative\n")


def _extract_alt_fasta(genome_fa: str, alt_bed: str, alt_fasta: str, logger: logging.Logger) -> None:
    """
    从alternative BED文件中提取序列。
    
    Args:
        genome_fa: 基因组FASTA文件路径
        alt_bed: alternative BED文件路径
        alt_fasta: 输出alternative FASTA文件路径
        logger: 日志记录器
    """
    cmd = f"bedtools getfasta -fi {genome_fa} -bed {alt_bed} -fo {alt_fasta}"
    utils.run_cmd(cmd, logger)


def _name_alternative_blocks(sorted_bed: str, out_named_bed: str) -> None:
    """
    为alternative区间添加编号（ChrXX-AB0000格式）。
    
    Args:
        sorted_bed: 排序后的BED文件路径（输入）
        out_named_bed: 命名后的BED文件路径（输出）
    """
    n = 0
    with open(sorted_bed) as fi, open(out_named_bed, "w") as fo:
        for line in fi:
            cols = line.strip().split()
            if len(cols) >= 4 and cols[3] == "alternative":
                ab = f"{cols[0]}-AB{n:04d}"
                fo.write("\t".join([cols[0], cols[1], cols[2], ab]) + "\n")
                n += 1
            else:
                fo.write(line)


def run_alternative_for_all_chr(
    genome_dir: str,
    output_dir: str,
    chr_list: List[str],
    ref_genome: str,
    resume: bool = False,
) -> None:
    """
    为所有染色体计算alternative block。
    基于core block的补集计算alternative区间，生成BED和FASTA文件。
    
    Args:
        genome_dir: 基因组目录（包含{name}.fasta）
        output_dir: 输出根目录
        chr_list: 染色体列表
        ref_genome: 参考基因组名称（未使用，保留兼容性）
        resume: 是否启用断点续跑
    """
    logger = logging.getLogger("panCB.alt")
    for chr_name in chr_list:
        coreB_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}", "coreB")
        if not os.path.isdir(coreB_dir):
            logger.warning("[%s] 缺少核心目录，跳过 alternative: %s", chr_name, coreB_dir)
            continue

        # 枚举样本（按 coreB 的 per-sample 坐标）
        names_set = set()
        for fname in os.listdir(coreB_dir):
            if fname.startswith(f"{chr_name}.Os.RP_qryRecord."):
                names_set.add(fname.split(".Os.RP_qryRecord.")[-1])

        if not names_set:
            logger.warning("[%s] 未找到任何 ChrXX.Os.RP_qryRecord.* 文件，跳过 alternative 计算。请确保步骤4.5已正确执行。", chr_name)
            continue

        logger.info("[%s] 找到 %d 个样本的 core 坐标文件，开始计算 alternative", chr_name, len(names_set))

        alterB_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}", "alterB")
        fasta_dir = os.path.join(alterB_dir, "alter-FASTA")
        bed_dir = os.path.join(alterB_dir, "alter-BED")

        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(bed_dir, exist_ok=True)

        for name in sorted(list(names_set)):
            genome_fa = os.path.join(genome_dir, f"{name}.fasta")
            utils.ensure_fai_exists(genome_fa, logger)
            fai_path = genome_fa + ".fai"

            core_bed = os.path.join(coreB_dir, f"{chr_name}.Os.RP_qryRecord.{name}")
            if not os.path.isfile(core_bed):
                logger.warning("[%s] 缺少 core 坐标：%s", chr_name, core_bed)
                continue

            sorted_bed = os.path.join(bed_dir, f"{chr_name}.RP_qryRecord.{name}.sorted.bed")
            len_path = os.path.join(bed_dir, f"{chr_name}.RP_qryRecord.{name}.sorted.len")
            alt_bed = os.path.join(bed_dir, f"{chr_name}.RP_qryRecord.{name}.alternative.bed")
            named_bed = os.path.join(bed_dir, f"{chr_name}.RP_qryRecord.{name}.sorted.bed.named")
            sorted_fasta = os.path.join(fasta_dir, f"{chr_name}.RP_qryRecord.{name}.sorted.fasta")
            alt_fasta = os.path.join(fasta_dir, f"{chr_name}.RP_qryRecord.{name}.alternative.fasta")

            if resume and os.path.isfile(named_bed) and os.path.isfile(alt_fasta):
                logger.info("[resume][%s] 跳过已有 alternative：%s", chr_name, name)
                continue

            chr_len = _read_chr_length_from_fai(fai_path, chr_name)
            if chr_len <= 0:
                logger.warning("[%s] fai 无法获取染色体长度：%s", chr_name, fai_path)
                continue

            core_intervals = _read_core_intervals(core_bed)
            records = _build_core_and_alternative(chr_len, core_intervals)
            _write_sorted_bed_and_len(records, chr_name, sorted_bed, len_path)

            # alternative-only BED and FASTA
            _write_alternative_bed(records, chr_name, alt_bed)
            _extract_alt_fasta(genome_fa, alt_bed, alt_fasta, logger)

            # sorted FASTA (both core & alternative) for compatibility
            _extract_sorted_fasta(genome_fa, sorted_bed, sorted_fasta, logger)

            # named BED from sorted BED
            _name_alternative_blocks(sorted_bed, named_bed)

            logger.info("[%s] alternative 生成完成：%s", chr_name, name)