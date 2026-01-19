import os
import re
import subprocess
import logging
import sys
from multiprocessing import Pool
from typing import List, Dict
import utils
import debug_log

# 移除Python 3.11+对整数转字符串的位数限制（用于debug日志输出大坐标值）
sys.set_int_max_str_digits(0)


# Merged from split_chromosomes.py
from Bio import SeqIO


def detect_chromosomes(ref_fasta: str) -> List[str]:
    """
    从参考基因组FASTA文件中检测所有染色体ID。
    
    Args:
        ref_fasta: 参考基因组FASTA文件路径
    
    Returns:
        排序后的染色体ID列表
    """
    chr_ids: List[str] = []
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        chr_ids.append(rec.id)
    chr_ids = sorted(chr_ids)
    return chr_ids


def _list_existing_chr_dirs(chr_split_root: str) -> List[str]:
    """
    列出已存在的染色体拆分目录。
    
    Args:
        chr_split_root: 染色体拆分根目录路径
    
    Returns:
        已存在的染色体目录名称列表（排序后）
    """
    if not os.path.isdir(chr_split_root):
        return []
    chr_ids = []
    for name in os.listdir(chr_split_root):
        path = os.path.join(chr_split_root, name)
        if os.path.isdir(path):
            chr_ids.append(name)
    return sorted(chr_ids)


def split_all_genomes_by_chr(genome_dir: str, accs_file: str, output_dir: str, ref_name: str, resume: bool = False) -> List[str]:
    """
    拆分为 {output_dir}/Chr_split/ChrXX/{genome}_ChrXX.fasta
    若某个体缺失该染色体，则写入1bp占位序列（N）。
    """
    logger = logging.getLogger("panCB.split")
    accs = utils.read_accs(accs_file)
    chr_split_root = os.path.join(output_dir, "Chr_split")

    if resume:
        existing_chr_ids = _list_existing_chr_dirs(chr_split_root)
        if existing_chr_ids:
            logger.info("[resume] 检测到已有拆分目录，跳过拆分: %s", chr_split_root)
            return existing_chr_ids

    utils.ensure_dir(chr_split_root)

    ref_fasta = os.path.join(genome_dir, f"{ref_name}.fasta")
    chr_ids = detect_chromosomes(ref_fasta)

    for chr_id in chr_ids:
        utils.ensure_dir(os.path.join(chr_split_root, chr_id))

    for name in accs:
        fasta_path = os.path.join(genome_dir, f"{name}.fasta")
        recs = {rec.id: rec for rec in SeqIO.parse(fasta_path, "fasta")}
        for chr_id in chr_ids:
            out_dir = os.path.join(chr_split_root, chr_id)
            out_fa = os.path.join(out_dir, f"{name}_{chr_id}.fasta")
            if resume and os.path.isfile(out_fa) and utils.fasta_has_bases(out_fa):
                logger.info("[resume] 跳过已有切片: %s", out_fa)
            elif chr_id in recs and len(recs[chr_id].seq) > 0:
                SeqIO.write(recs[chr_id], out_fa, "fasta")
            else:
                logger.warning("样本 %s 缺失 %s，写入占位序列", name, chr_id)
                utils.write_minimal_placeholder_fasta(out_fa, chr_id)

    return chr_ids


def core_block_filter(filter_coords: str, ilen: int, out_core_dir: str, chr_name: str = "", debug_log=None, core_identity: float = 98.0) -> str:
    """
    对过滤后的比对坐标进行核心块初筛。
    过滤规则：身份度>=core_identity、非反向、长度>ilen、去除重叠区间。
    
    Args:
        filter_coords: 过滤后的coords文件路径（show-coords输出）
        ilen: 最小长度阈值
        out_core_dir: 输出目录（coreB目录）
        chr_name: 染色体名称（用于debug日志）
        debug_log: debug日志文件句柄
        core_identity: 身份度阈值（默认98.0%）
    
    Returns:
        初筛结果文件路径（*.filtered.1227）
    """
    name = os.path.basename(filter_coords)
    
    if debug_log:
        debug_log.write(f"\n{'='*80}\n")
        debug_log.write(f"【阶段1】初筛过滤 - {name}\n")
        debug_log.write(f"{'='*80}\n")
        debug_log.write(f"输入文件: {filter_coords}\n")
        debug_log.write(f"最小长度阈值: {ilen} bp\n")
        debug_log.write(f"身份度阈值: {core_identity}%\n\n")

    qryd = {}
    total_blocks = 0
    filtered_stats = {
        'chr_mismatch': 0,
        'identity_low': 0,
        'reverse': 0,
        'length_short': 0,
        'passed': 0
    }
    
    with open(filter_coords) as f:
        for i in f:
            ii = re.split('\t', i.strip())
            total_blocks += 1
            
            # 条件1: 染色体ID一致
            if ii[9] != ii[10]:
                filtered_stats['chr_mismatch'] += 1
                if debug_log:
                    debug_log.write(f"Block #{total_blocks} 被过滤 [染色体不一致]: ref={ii[9]}, qry={ii[10]}, "
                                   f"ref_pos={ii[0]}-{ii[1]}, qry_pos={ii[2]}-{ii[3]}\n")
                continue
            
            # 条件2: 一致性≥core_identity
            if float(ii[6]) < core_identity:
                filtered_stats['identity_low'] += 1
                if debug_log:
                    debug_log.write(f"Block #{total_blocks} 被过滤 [一致性<{core_identity}%]: identity={ii[6]}%, "
                                   f"chr={ii[9]}, ref_pos={ii[0]}-{ii[1]}, qry_pos={ii[2]}-{ii[3]}\n")
                continue
            
            # 条件3: 正向比对
            if int(ii[2]) >= int(ii[3]):
                filtered_stats['reverse'] += 1
                if debug_log:
                    debug_log.write(f"Block #{total_blocks} 被过滤 [反向比对]: chr={ii[9]}, "
                                   f"ref_pos={ii[0]}-{ii[1]}, qry_pos={ii[2]}-{ii[3]}\n")
                continue
            
            # 条件4: 长度过滤
            if int(ii[4]) <= int(ilen) and int(ii[5]) <= int(ilen):
                filtered_stats['length_short'] += 1
                if debug_log:
                    debug_log.write(f"Block #{total_blocks} 被过滤 [长度不足]: chr={ii[9]}, "
                                   f"ref_len={ii[4]}bp, qry_len={ii[5]}bp, threshold={ilen}bp, "
                                   f"ref_pos={ii[0]}-{ii[1]}, qry_pos={ii[2]}-{ii[3]}\n")
                continue
            
            # 通过所有条件
            filtered_stats['passed'] += 1
            if ii[9] in qryd:
                qryd[ii[9]].append([int(ii[0]), int(ii[1]), int(ii[2]), int(ii[3])])
            else:
                qryd[ii[9]] = [[int(ii[0]), int(ii[1]), int(ii[2]), int(ii[3])]]
    
    if debug_log:
        debug_log.write(f"\n--- 初筛统计 ---\n")
        debug_log.write(f"总输入blocks: {total_blocks}\n")
        debug_log.write(f"  - 染色体不一致: {filtered_stats['chr_mismatch']}\n")
        debug_log.write(f"  - 一致性<{core_identity}%: {filtered_stats['identity_low']}\n")
        debug_log.write(f"  - 反向比对: {filtered_stats['reverse']}\n")
        debug_log.write(f"  - 长度不足: {filtered_stats['length_short']}\n")
        debug_log.write(f"  - 通过初筛: {filtered_stats['passed']}\n\n")

    # 预处理：过滤掉无效区间（start >= end）
    invalid_interval_count = 0
    for chr_key in qryd.keys():
        valid_blocks = []
        for block in qryd[chr_key]:
            # 检查ref和qry区间是否有效（start < end）
            if block[0] < block[1] and block[2] < block[3]:
                valid_blocks.append(block)
            else:
                invalid_interval_count += 1
                if debug_log:
                    debug_log.write(f"  过滤无效区间: chr={chr_key}, ref={block[0]}-{block[1]}, qry={block[2]}-{block[3]}\n")
        qryd[chr_key] = valid_blocks
    
    if invalid_interval_count > 0:
        if debug_log:
            debug_log.write(f"预处理过滤无效区间: {invalid_interval_count} 个\n\n")
        logging.getLogger("panCB.align").info("预处理过滤无效区间: %d 个", invalid_interval_count)

    # 参考序列坐标去重叠
    ref_overlap_count = 0
    if debug_log:
        debug_log.write(f"--- 参考序列坐标去重叠 ---\n")
    
    for i in qryd.keys():
        endn = 1
        max_iterations = len(qryd[i]) * 100  # 防止无限循环的安全阈值
        iteration_count = 0
        while endn != 0:
            endn = 0
            iteration_count += 1
            if iteration_count > max_iterations:
                if debug_log:
                    debug_log.write(f"  警告: 参考序列去重叠迭代次数超过阈值({max_iterations})，强制退出\n")
                logging.getLogger("panCB.align").warning("[%s] 参考序列去重叠迭代次数超过阈值，强制退出", i)
                break
            for c in range(len(qryd[i][:]) - 1):
                intetmp = utils.intersection(qryd[i][c][0:2], qryd[i][c + 1][0:2])
                if len(intetmp) != 0:
                    intetmpl = intetmp[1] - intetmp[0]
                    # 跳过无效重叠（overlap_len <= 0 或区间无效）
                    if intetmpl <= 0:
                        continue
                    # 检查区间有效性（start < end）
                    if qryd[i][c][0] >= qryd[i][c][1] or qryd[i][c + 1][0] >= qryd[i][c + 1][1]:
                        if debug_log:
                            debug_log.write(f"  跳过无效区间: block1_ref={qryd[i][c][0]}-{qryd[i][c][1]}, "
                                           f"block2_ref={qryd[i][c+1][0]}-{qryd[i][c+1][1]}\n")
                        continue
                    ref_overlap_count += 1
                    if debug_log:
                        debug_log.write(f"  重叠 #{ref_overlap_count}: chr={i}, overlap_len={intetmpl}bp, "
                                       f"block1_ref={qryd[i][c][0]}-{qryd[i][c][1]}, "
                                       f"block2_ref={qryd[i][c+1][0]}-{qryd[i][c+1][1]}\n")
                    qryd[i].insert(c, [qryd[i][c][0], qryd[i][c][1] - intetmpl - 1, qryd[i][c][2], qryd[i][c][3] - intetmpl - 1])
                    del qryd[i][c + 1]
                    qryd[i].insert(c + 1, [qryd[i][c + 1][0] + intetmpl, qryd[i][c + 1][1], qryd[i][c + 1][2] + intetmpl, qryd[i][c + 1][3]])
                    del qryd[i][c + 2]
                    endn += 1
    
    if debug_log:
        debug_log.write(f"参考序列重叠处理: {ref_overlap_count} 次\n\n")

    # 查询序列坐标去重叠
    qry_overlap_count = 0
    if debug_log:
        debug_log.write(f"--- 查询序列坐标去重叠 ---\n")
    
    for i in qryd.keys():
        endn = 1
        max_iterations = len(qryd[i]) * 100  # 防止无限循环的安全阈值
        iteration_count = 0
        while endn != 0:
            endn = 0
            iteration_count += 1
            if iteration_count > max_iterations:
                if debug_log:
                    debug_log.write(f"  警告: 查询序列去重叠迭代次数超过阈值({max_iterations})，强制退出\n")
                logging.getLogger("panCB.align").warning("[%s] 查询序列去重叠迭代次数超过阈值，强制退出", i)
                break
            for c in range(len(qryd[i]) - 1):
                intetmp = utils.intersection(qryd[i][c][2:4], qryd[i][c + 1][2:4])
                if len(intetmp) != 0:
                    intetmpl = intetmp[1] - intetmp[0]
                    # 跳过无效重叠（overlap_len <= 0 或区间无效）
                    if intetmpl <= 0:
                        continue
                    # 检查区间有效性（start < end）
                    if qryd[i][c][2] >= qryd[i][c][3] or qryd[i][c + 1][2] >= qryd[i][c + 1][3]:
                        if debug_log:
                            debug_log.write(f"  跳过无效区间: block1_qry={qryd[i][c][2]}-{qryd[i][c][3]}, "
                                           f"block2_qry={qryd[i][c+1][2]}-{qryd[i][c+1][3]}\n")
                        continue
                    qry_overlap_count += 1
                    if debug_log:
                        debug_log.write(f"  重叠 #{qry_overlap_count}: chr={i}, overlap_len={intetmpl}bp, "
                                       f"block1_qry={qryd[i][c][2]}-{qryd[i][c][3]}, "
                                       f"block2_qry={qryd[i][c+1][2]}-{qryd[i][c+1][3]}\n")
                    qryd[i].insert(c, [qryd[i][c][0], qryd[i][c][1] - intetmpl - 1, qryd[i][c][2], qryd[i][c][3] - intetmpl - 1])
                    del qryd[i][c + 1]
                    qryd[i].insert(c + 1, [qryd[i][c + 1][0] + intetmpl, qryd[i][c + 1][1], qryd[i][c + 1][2] + intetmpl, qryd[i][c + 1][3]])
                    del qryd[i][c + 2]
                    endn += 1
    
    if debug_log:
        debug_log.write(f"查询序列重叠处理: {qry_overlap_count} 次\n\n")

    utils.ensure_dir(out_core_dir)
    out_path = os.path.join(out_core_dir, f"{name}.{ilen}.filtered.1227")
    
    final_output_count = 0
    final_filtered_count = 0
    
    if debug_log:
        debug_log.write(f"--- 最终长度过滤 ---\n")
    
    with open(out_path, "w") as rout:
        for i in qryd.keys():
            for j in qryd[i]:
                if j[1] - j[0] > int(ilen) or j[3] - j[2] > int(ilen):
                    rout.write(i + "\t" + str(j[0]) + "\t" + str(j[1]) + "\t" + str(j[2]) + "\t" + str(j[3]) + "\n")
                    final_output_count += 1
                else:
                    final_filtered_count += 1
                    if debug_log:
                        debug_log.write(f"  Block 被过滤 [最终长度不足]: chr={i}, "
                                       f"ref_len={j[1]-j[0]}bp, qry_len={j[3]-j[2]}bp, threshold={ilen}bp\n")
    
    if debug_log:
        debug_log.write(f"\n--- 初筛最终结果 ---\n")
        debug_log.write(f"最终输出blocks: {final_output_count}\n")
        debug_log.write(f"最终过滤blocks: {final_filtered_count}\n")
        debug_log.write(f"输出文件: {out_path}\n")
        debug_log.write(f"{'='*80}\n\n")

    logging.getLogger("panCB.align").info("初筛完成: %s", out_path)
    return out_path


def _run_single_alignment_chr(args_tuple) -> str:
    """
    对单个查询样本与参考基因组在指定染色体上进行比对和初筛。
    
    Args:
        args_tuple: 参数元组，包含参考名称、查询名称、染色体名称、拆分根目录、
                   输出目录、nucmer参数、身份度、最小长度、resume标志、debug标志、core_identity
    
    Returns:
        初筛结果文件路径（*.filtered.1227），失败返回空字符串
    """
    # 支持13或14个参数（向后兼容）
    if len(args_tuple) == 14:
        (
            ref_name,
            query_name,
            chr_name,
            split_root,
            output_dir,
            min_match,
            min_cluster,
            max_gap,
            nucmer_threads,
            identity,
            min_len,
            resume,
            debug,
            core_identity,
        ) = args_tuple
    else:
        (
            ref_name,
            query_name,
            chr_name,
            split_root,
            output_dir,
            min_match,
            min_cluster,
            max_gap,
            nucmer_threads,
            identity,
            min_len,
            resume,
            debug,
        ) = args_tuple
        core_identity = 98.0  # 默认值

    logger = logging.getLogger("panCB.align")

    ref_fa = os.path.join(split_root, chr_name, f"{ref_name}_{chr_name}.fasta")
    qry_fa = os.path.join(split_root, chr_name, f"{query_name}_{chr_name}.fasta")

    if not os.path.isfile(ref_fa) or not utils.fasta_has_bases(ref_fa):
        logger.warning("[%s] 跳过：参考切片缺失或为空: %s", chr_name, ref_fa)
        return ""
    if not os.path.isfile(qry_fa) or not utils.fasta_has_bases(qry_fa):
        logger.warning("[%s] 跳过：查询切片缺失或为空: %s", chr_name, qry_fa)
        return ""

    gf_chr_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}")
    utils.ensure_dir(gf_chr_dir)
    coreB_dir = os.path.join(gf_chr_dir, "coreB")
    utils.ensure_dir(coreB_dir)

    prefix = os.path.join(gf_chr_dir, f"{ref_name}_RP_{query_name}")

    delta_path = f"{prefix}.delta"
    filtered_delta = f"{prefix}.filtered.delta"
    coords_path = f"{prefix}.filtered.coords"
    expected_1227 = os.path.join(coreB_dir, os.path.basename(coords_path) + f".{min_len}.filtered.1227")

    # 细粒度断点续跑：检查各阶段产物
    if resume and os.path.isfile(expected_1227) and utils.file_size(expected_1227) > 0:
        logger.info("[resume][%s] 跳过全部步骤，.filtered.1227 已存在: %s", chr_name, expected_1227)
        return expected_1227

    # 检查是否可以从 coords 文件继续（跳过 nucmer + delta-filter + show-coords）
    coords_exists = os.path.isfile(coords_path) and os.path.getsize(coords_path) > 0
    filtered_delta_exists = os.path.isfile(filtered_delta) and os.path.getsize(filtered_delta) > 0
    delta_exists = os.path.isfile(delta_path) and os.path.getsize(delta_path) > 0
    
    if resume:
        logger.info("[resume][%s][%s] 中间文件检测: delta=%s, filtered_delta=%s, coords=%s", 
                   chr_name, query_name, delta_exists, filtered_delta_exists, coords_exists)
    
    if resume and coords_exists:
        logger.info("[resume][%s] 检测到 coords 文件，跳过比对步骤，直接进行初筛: %s", chr_name, coords_path)
    else:
        # 检查是否可以从 filtered.delta 继续（跳过 nucmer + delta-filter）
        if resume and filtered_delta_exists:
            logger.info("[resume][%s] 检测到 filtered.delta，跳过 nucmer 和 delta-filter: %s", chr_name, filtered_delta)
        else:
            # 检查是否可以从 delta 继续（跳过 nucmer）
            if resume and delta_exists:
                logger.info("[resume][%s] 检测到 delta 文件，跳过 nucmer: %s", chr_name, delta_path)
            else:
                # 从头开始：运行 nucmer
                nucmer_cmd = (
                    f"nucmer -l {min_match} -c {min_cluster} -g {max_gap} -t {nucmer_threads} "
                    f"-p {prefix} {ref_fa} {qry_fa}"
                )
                ret = utils.run_cmd(nucmer_cmd, logger)
                if ret != 0:
                    logger.warning("[%s] nucmer 非零退出，跳过: %s", chr_name, prefix)
                    return ""

            # 运行 delta-filter
            df_cmd = f"delta-filter -g -r -q -i {identity} -l {min_match} {delta_path} > {filtered_delta}"
            utils.run_cmd(df_cmd, logger)

        # 运行 show-coords
        sc_cmd = f"show-coords -THrd {filtered_delta} > {coords_path}"
        utils.run_cmd(sc_cmd, logger)

    if not os.path.isfile(coords_path) or os.path.getsize(coords_path) == 0:
        logger.warning("[%s] coords 为空，跳过初筛: %s", chr_name, coords_path)
        return ""

    # 使用 debug_log 模块管理调试日志
    with debug_log.create_debug_logger(output_dir, chr_name, debug) as dbg_logger:
        filtered_1227 = core_block_filter(coords_path, min_len, coreB_dir, chr_name, dbg_logger.file, core_identity)
        if filtered_1227:
            lines = utils.count_lines(filtered_1227)
            sizeb = utils.file_size(filtered_1227)
            logger.info("[%s] 初筛文件: %s 行数=%d 大小=%dB (identity_threshold=%.1f%%)", 
                       chr_name, filtered_1227, lines, sizeb, core_identity)
        return filtered_1227


def run_all_alignments_for_chr(
    ref_name: str,
    chr_name: str,
    split_root: str,
    output_dir: str,
    accs_file: str,
    processes: int,  # 保留用于兼容性，本函数内部不使用
    resume: bool = False,
    debug: bool = False,
    core_identity: float = 98.0,
) -> List[str]:
    """
    对指定染色体的所有查询样本进行比对和初筛。

    注意：此函数在染色体级别的并行进程池中被调用，因此内部采用顺序执行，
    避免嵌套多进程导致的"daemonic processes are not allowed to have children"错误。

    Args:
        ref_name: 参考基因组名称
        chr_name: 染色体名称
        split_root: 染色体拆分根目录（Chr_split）
        output_dir: 输出根目录
        accs_file: 样本列表文件路径
        processes: 并行进程数（此参数在本函数中不使用，保留用于兼容性）
        resume: 是否启用断点续跑
        debug: 是否启用调试模式
        core_identity: 核心block身份度阈值（默认98.0%）

    Returns:
        所有查询样本的初筛结果文件路径列表（*.filtered.1227）
    """
    # 硬编码参数值（需求 2.2-2.7）
    min_match = 50          # nucmer 最小匹配长度
    min_cluster = 100       # nucmer 聚簇阈值
    max_gap = 30            # nucmer 允许缺口
    nucmer_threads = 20     # nucmer 线程数
    identity = 90           # delta-filter 身份度阈值（%）
    min_len = 50            # 核心块最小长度阈值
    
    logger = logging.getLogger("panCB.align")
    accs = utils.read_accs(accs_file)
    queries = [x for x in accs if x != ref_name]
    if not queries:
        logger.warning("[%s] 无需比对：仅参考基因组", chr_name)
        return []

    chr_root = os.path.join(output_dir, f"genomeFilter_{chr_name}")
    coreB_dir = os.path.join(chr_root, "coreB")
    if resume and os.path.isdir(coreB_dir):
        existing_filtered = [
            os.path.join(coreB_dir, f)
            for f in os.listdir(coreB_dir)
            if f.endswith(".filtered.1227")
        ]
        if existing_filtered:
            logger.info("[resume][%s] 检测到已有比对目录，跳过 nucmer: %s", chr_name, coreB_dir)
            return sorted(existing_filtered)

    tasks = [
        (
            ref_name,
            q,
            chr_name,
            split_root,
            output_dir,
            min_match,
            min_cluster,
            max_gap,
            nucmer_threads,
            identity,
            min_len,
            resume,
            debug,
            core_identity,
        )
        for q in queries
    ]

    # 顺序执行所有样本的比对任务，避免嵌套多进程
    # 并行已经在染色体级别实现（panCB.py中的Pool处理多条染色体）
    results: List[str] = []
    logger.info("[%s] 开始顺序处理 %d 个样本的比对任务 (core_identity=%.1f%%)", chr_name, len(tasks), core_identity)
    for i, t in enumerate(tasks, 1):
        query_name = t[1]  # 查询样本名称
        logger.info("[%s] 正在处理样本 %d/%d: %s", chr_name, i, len(tasks), query_name)
        out_path = _run_single_alignment_chr(t)
        if out_path:
            results.append(out_path)
    
    logger.info("[%s] 完成所有样本比对，成功 %d/%d", chr_name, len(results), len(tasks))
    return results


def run_all_alignments_parallel(
    ref_name: str,
    chr_list: List[str],
    split_root: str,
    output_dir: str,
    accs_file: str,
    processes: int,
    resume: bool = False,
    debug: bool = False,
    core_identity: float = 98.0,
) -> Dict[str, List[str]]:
    """
    并行处理所有染色体的所有样本比对（任务级并行）。
    
    这是优化后的比对策略：将所有 (样本, 染色体) 对作为独立任务并行处理，
    而不是按染色体分组处理。这样可以最大化CPU利用率，避免长染色体导致的资源闲置。
    
    Args:
        ref_name: 参考基因组名称
        chr_list: 染色体列表
        split_root: 染色体拆分根目录（Chr_split）
        output_dir: 输出根目录
        accs_file: 样本列表文件路径
        processes: 并行进程数
        resume: 是否启用断点续跑
        debug: 是否启用调试模式
        core_identity: 核心block身份度阈值（默认98.0%）
    
    Returns:
        字典，key为染色体名称，value为该染色体的初筛结果文件路径列表
        例如: {'Chr01': ['/path/to/file1.filtered.1227', ...], 'Chr02': [...], ...}
    """
    # 硬编码参数值
    min_match = 50
    min_cluster = 100
    max_gap = 30
    nucmer_threads = 20
    identity = 90
    min_len = 50
    
    logger = logging.getLogger("panCB.align")
    accs = utils.read_accs(accs_file)
    queries = [x for x in accs if x != ref_name]
    
    if not queries:
        logger.warning("无需比对：仅参考基因组")
        return {chr_name: [] for chr_name in chr_list}
    
    # 检查resume模式，收集已完成的任务
    completed_tasks = set()
    chr_to_results = {chr_name: [] for chr_name in chr_list}
    
    if resume:
        for chr_name in chr_list:
            chr_root = os.path.join(output_dir, f"genomeFilter_{chr_name}")
            coreB_dir = os.path.join(chr_root, "coreB")
            if os.path.isdir(coreB_dir):
                for f in os.listdir(coreB_dir):
                    if f.endswith(".filtered.1227"):
                        # 从文件名提取样本名称
                        # 格式: {ref_name}_RP_{query_name}.filtered.coords.{min_len}.filtered.1227
                        match = re.search(rf'{ref_name}_RP_(.+?)\.filtered\.coords', f)
                        if match:
                            query_name = match.group(1)
                            completed_tasks.add((query_name, chr_name))
                            chr_to_results[chr_name].append(os.path.join(coreB_dir, f))
        
        if completed_tasks:
            logger.info("[resume] 检测到 %d 个已完成的比对任务，将跳过", len(completed_tasks))
    
    # 生成所有待处理的任务：(样本, 染色体) 对
    all_tasks = []
    for chr_name in chr_list:
        for query_name in queries:
            if (query_name, chr_name) not in completed_tasks:
                all_tasks.append(
                    (
                        ref_name,
                        query_name,
                        chr_name,
                        split_root,
                        output_dir,
                        min_match,
                        min_cluster,
                        max_gap,
                        nucmer_threads,
                        identity,
                        min_len,
                        resume,
                        debug,
                        core_identity,
                    )
                )
    
    total_tasks = len(queries) * len(chr_list)
    pending_tasks = len(all_tasks)
    completed_count = total_tasks - pending_tasks
    
    logger.info("=" * 80)
    logger.info("阶段1：并行比对所有样本×染色体")
    logger.info("  总任务数: %d (样本数: %d × 染色体数: %d)", total_tasks, len(queries), len(chr_list))
    logger.info("  已完成: %d", completed_count)
    logger.info("  待处理: %d", pending_tasks)
    logger.info("  进程池大小: %d", processes)
    logger.info("  核心block身份度阈值: %.1f%%", core_identity)
    logger.info("  策略: 任务级并行（自动负载均衡）")
    logger.info("=" * 80)
    
    if not all_tasks:
        logger.info("[resume] 所有比对任务已完成，跳过比对阶段")
        return chr_to_results
    
    # 使用进程池并行处理所有任务
    processed_count = 0
    if processes and processes > 1:
        with Pool(processes) as pool:
            for out_path in pool.imap_unordered(_run_single_alignment_chr, all_tasks):
                if out_path:
                    # 从路径中提取染色体名称
                    # 路径格式: .../genomeFilter_ChrXX/coreB/...
                    match = re.search(r'genomeFilter_(Chr\d+)', out_path)
                    if match:
                        chr_name = match.group(1)
                        chr_to_results[chr_name].append(out_path)
                
                processed_count += 1
                if processed_count % 50 == 0 or processed_count == pending_tasks:
                    progress = (completed_count + processed_count) / total_tasks * 100
                    logger.info("比对进度: %d/%d (%.1f%%) - 已完成 %d 个任务",
                               completed_count + processed_count, total_tasks, progress, processed_count)
    else:
        # 顺序执行（用于调试）
        for i, task in enumerate(all_tasks, 1):
            out_path = _run_single_alignment_chr(task)
            if out_path:
                match = re.search(r'genomeFilter_(Chr\d+)', out_path)
                if match:
                    chr_name = match.group(1)
                    chr_to_results[chr_name].append(out_path)
            
            processed_count += 1
            if processed_count % 50 == 0 or processed_count == pending_tasks:
                progress = (completed_count + processed_count) / total_tasks * 100
                logger.info("比对进度: %d/%d (%.1f%%)",
                           completed_count + processed_count, total_tasks, progress)
    
    logger.info("=" * 80)
    logger.info("阶段1完成：所有比对任务已完成")
    logger.info("  总计: %d/%d 个任务成功", completed_count + processed_count, total_tasks)
    
    # 统计每条染色体的结果
    for chr_name in chr_list:
        result_count = len(chr_to_results[chr_name])
        logger.info("  [%s] %d 个样本比对完成", chr_name, result_count)
    
    logger.info("=" * 80)
    
    # 对每条染色体的结果进行排序
    for chr_name in chr_list:
        chr_to_results[chr_name] = sorted(chr_to_results[chr_name])
    
    return chr_to_results
