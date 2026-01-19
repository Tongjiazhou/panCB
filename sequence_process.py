import os
import re
import logging
from typing import List
import utils


def can_regenerate_chr_coords(output_dir: str, chr_name: str, min_len: int = 50) -> bool:
    """
    检查是否可以从迭代结果重新生成染色体级别坐标文件。
    
    判断条件：
    1. 最终迭代结果文件存在
    2. 文件非空且有效（至少有一行数据）
    
    Args:
        output_dir: 输出根目录
        chr_name: 染色体名称
        min_len: 最小长度阈值（默认50）
    
    Returns:
        如果迭代结果存在且可用于重新生成，返回True；否则返回False
    
    _Requirements: 2.1, 2.2_
    """
    logger = logging.getLogger("panCB.seq")
    
    coreB_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}", "coreB")
    if not os.path.isdir(coreB_dir):
        logger.debug("[%s] coreB目录不存在: %s", chr_name, coreB_dir)
        return False
    
    # 查找最终迭代结果文件
    candidates = [
        f for f in os.listdir(coreB_dir)
        if f.startswith(f"{chr_name}.{min_len}.R") and f.endswith("_RP_qryRecord.coords")
    ]
    
    if not candidates:
        logger.debug("[%s] 未找到迭代结果文件", chr_name)
        return False
    
    # 按迭代轮次数字排序
    def extract_iter_num(filename):
        match = re.search(r'\.R(\d+)_RP_qryRecord\.coords$', filename)
        return int(match.group(1)) if match else 0
    
    candidates.sort(key=extract_iter_num)
    final_file = os.path.join(coreB_dir, candidates[-1])
    
    # 检查文件是否非空
    if utils.file_size(final_file) <= 0:
        logger.debug("[%s] 最终迭代结果文件为空: %s", chr_name, final_file)
        return False
    
    # 检查文件是否有效（至少有一行有效数据）
    try:
        with open(final_file) as f:
            for line in f:
                cols = line.strip().split()
                # 有效行至少需要3列：chr ref_start ref_end
                if len(cols) >= 3:
                    logger.debug("[%s] 可以从迭代结果重新生成坐标文件: %s", chr_name, final_file)
                    return True
        logger.debug("[%s] 最终迭代结果文件无有效数据行: %s", chr_name, final_file)
        return False
    except Exception as e:
        logger.debug("[%s] 读取迭代结果文件失败: %s, 错误: %s", chr_name, final_file, e)
        return False


def extract_chr_level_coords(
    output_dir: str,
    chr_name: str,
    ref_genome: str,
    resume: bool = False,
    force_regenerate: bool = False,
) -> bool:
    """
    为单个染色体生成染色体级别的核心block坐标结果。
    输出到: output/genomeFilter_ChrXX/coreB/ChrXX.Os.RP_qryRecord.{sample}
    
    关键逻辑说明：
    - 最终迭代结果文件 Rn_RP_qryRecord.coords 的列结构：
      第1列=染色体, 第2-3列=参考坐标, 第4-5列=第1个查询样本坐标, ...
    - 根据最终文件的实际列数来确定有多少样本参与了完整迭代
    - 只为实际参与迭代的样本生成坐标文件
    
    Args:
        output_dir: 输出根目录
        chr_name: 染色体名称
        ref_genome: 参考基因组名称
        resume: 是否启用断点续跑模式（跳过已存在的文件）
        force_regenerate: 是否强制重新生成（即使文件已存在）
    
    Returns:
        True表示成功，False表示失败
    
    _Requirements: 2.2, 2.3_
    """
    logger = logging.getLogger("panCB.seq")
    
    # 硬编码参数：核心块最小长度阈值
    min_len = 50
    
    coreB_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}", "coreB")
    if not os.path.isdir(coreB_dir):
        logger.warning("[%s] 缺少核心目录: %s", chr_name, coreB_dir)
        return False
    
    # 找到最终的迭代结果文件
    candidates = [
        f for f in os.listdir(coreB_dir) 
        if f.startswith(f"{chr_name}.{min_len}.R") and f.endswith("_RP_qryRecord.coords")
    ]
    if not candidates:
        logger.warning("[%s] 未找到最终 R 结果，跳过。", chr_name)
        return False
    
    # 按迭代轮次数字排序（而不是字符串排序）
    # 文件名格式：ChrXX.50.R{n}_RP_qryRecord.coords
    def extract_iter_num(filename):
        """从文件名中提取迭代轮次数字"""
        import re
        match = re.search(r'\.R(\d+)_RP_qryRecord\.coords$', filename)
        return int(match.group(1)) if match else 0
    
    candidates.sort(key=extract_iter_num)
    final_file = os.path.join(coreB_dir, candidates[-1])
    logger.info("[%s] 选择最终迭代文件: %s (共%d个候选)", chr_name, candidates[-1], len(candidates))
    
    # 读取迭代顺序文件
    iter_file = os.path.join(coreB_dir, f"{ref_genome}.{chr_name}.iter")
    if not os.path.isfile(iter_file):
        logger.warning("[%s] 缺少迭代顺序文件: %s", chr_name, iter_file)
        return False
    
    # 首先检测最终文件的实际列数，确定有多少样本参与了完整迭代
    actual_col_count = 0
    with open(final_file) as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) > actual_col_count:
                actual_col_count = len(cols)
    
    # 根据列数计算实际参与迭代的样本数
    # 列结构：1(chr) + 2(ref) + 2*(query_count) = total_cols
    # 所以 query_count = (total_cols - 1) / 2 - 1 = (total_cols - 3) / 2
    if actual_col_count < 3:
        logger.warning("[%s] 最终迭代文件列数异常: %d", chr_name, actual_col_count)
        return False
    
    actual_sample_count = (actual_col_count - 1) // 2  # 包括参考基因组
    actual_query_count = actual_sample_count - 1  # 不包括参考基因组
    
    logger.info("[%s] 最终迭代文件列数=%d, 实际样本数=%d (参考+%d个查询)", 
                chr_name, actual_col_count, actual_sample_count, actual_query_count)
    
    # 读取.iter文件中的样本顺序
    iter_samples: List[str] = []
    for line in open(iter_file):
        name = line.strip().split()[0]
        if name:
            iter_samples.append(name)
    
    # 构建样本名称到列索引的映射，只包含实际参与迭代的样本
    # 参考基因组固定使用第1-2列（0-indexed: cols[1], cols[2]）
    iterorder = {ref_genome: [1, 2]}
    name_order: List[str] = [ref_genome]
    
    # 只处理实际参与迭代的查询样本
    for i, name in enumerate(iter_samples):
        if i >= actual_query_count:
            # 这个样本没有参与完整迭代，跳过
            logger.debug("[%s] 样本 %s 未参与完整迭代（索引=%d, 实际查询数=%d）", 
                        chr_name, name, i, actual_query_count)
            continue
        # 第i个查询样本的列索引：3+2*i 和 4+2*i（0-indexed）
        sidx = 3 + 2 * i
        eidx = 4 + 2 * i
        iterorder[name] = [sidx, eidx]
        name_order.append(name)
    
    logger.info("[%s] 将为 %d 个样本生成坐标文件（参考=%s, 查询=%d个）", 
                chr_name, len(name_order), ref_genome, len(name_order) - 1)
    
    # 如果是强制重新生成模式，记录日志
    if force_regenerate:
        logger.info("[regenerate][%s] 强制重新生成染色体级别坐标文件", chr_name)
    
    # 为每个样本生成坐标文件到coreB目录
    success_count = 0
    empty_count = 0
    regenerated_count = 0
    for name in name_order:
        # 输出到: genomeFilter_ChrXX/coreB/ChrXX.Os.RP_qryRecord.{sample}
        coreB_save_path = os.path.join(coreB_dir, f"{chr_name}.Os.RP_qryRecord.{name}")
        
        # 检查文件是否已存在
        file_exists = os.path.isfile(coreB_save_path) and utils.file_size(coreB_save_path) > 0
        
        # resume模式下跳过已存在的文件（除非force_regenerate为True）
        if resume and file_exists and not force_regenerate:
            logger.info("[resume][%s] 跳过已有染色体级别坐标: %s", chr_name, coreB_save_path)
            success_count += 1
        else:
            # 如果是强制重新生成且文件已存在，记录重新生成日志
            if force_regenerate and file_exists:
                regenerated_count += 1
            try:
                # 读取并处理数据，添加block编号
                coord_lines = []
                sidx, eidx = iterorder[name]
                block_num = 1  # 初始化block编号
                for line in open(final_file):
                    cols = line.strip().split()
                    if len(cols) < 3:
                        continue
                    if len(cols) <= max(sidx, eidx):
                        # 这一行的列数不够，记录警告但继续处理其他行
                        logger.debug("[%s] 样本 %s 在某行列数不足: 需要索引%d, 实际列数%d", 
                                    chr_name, name, max(sidx, eidx), len(cols))
                        continue
                    # 添加block编号到染色体名称：Chr01_block_1
                    chr_with_block = f"{cols[0]}_block_{block_num}"
                    new_cols = [chr_with_block, cols[sidx], cols[eidx]]
                    coord_lines.append("\t".join(new_cols) + "\n")
                    block_num += 1  # 递增block编号
                
                # 写入到coreB目录
                with open(coreB_save_path, "w") as save:
                    save.writelines(coord_lines)
                
                if utils.file_size(coreB_save_path) > 0:
                    success_count += 1
                    logger.info("[%s] 染色体级别坐标生成: %s 行=%d 大小=%dB", 
                              chr_name, coreB_save_path, 
                              utils.count_lines(coreB_save_path), utils.file_size(coreB_save_path))
                else:
                    empty_count += 1
                    logger.warning("[%s] 生成的坐标文件为空: %s (索引=[%d,%d])", 
                                  chr_name, coreB_save_path, sidx, eidx)
            except Exception as e:
                logger.error("[%s] 生成染色体级别坐标失败 %s: %s", chr_name, coreB_save_path, e)
    
    if success_count == len(name_order):
        if regenerated_count > 0:
            logger.info("[regenerate][%s] 染色体级别核心block坐标重新生成完成，共%d个样本（重新生成%d个）", 
                       chr_name, success_count, regenerated_count)
        else:
            logger.info("[%s] 染色体级别核心block坐标生成完成，共%d个样本", chr_name, success_count)
        return True
    else:
        logger.warning("[%s] 染色体级别坐标生成不完整: 成功%d/%d, 空文件%d", 
                      chr_name, success_count, len(name_order), empty_count)
        return success_count > 0  # 只要有成功的就返回True


def merge_chr_coords_and_extract_sequences(
    genome_dir: str,
    output_dir: str,
    chr_list: List[str],
    resume: bool = False,
) -> None:
    """
    从染色体级别结果合并生成全基因组级别的核心block坐标，并提取序列。
    染色体级别结果位于: output/genomeFilter_ChrXX/coreB/ChrXX.Os.RP_qryRecord.{sample}
    基因组级别结果输出到: output/pan_core_block/RP_qryRecord.{sample} 和 fasta/
    """
    logger = logging.getLogger("panCB.seq")

    # 收集所有样本名称（从染色体级别结果中）
    names_set = set()
    for chr_name in chr_list:
        coreB_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}", "coreB")
        if not os.path.isdir(coreB_dir):
            logger.warning("[%s] 核心目录不存在: %s", chr_name, coreB_dir)
            continue
        for fname in os.listdir(coreB_dir):
            if fname.startswith(f"{chr_name}.Os.RP_qryRecord."):
                names_set.add(fname.split(".Os.RP_qryRecord.")[-1])

    if not names_set:
        logger.warning("未找到任何染色体级别结果，无法合并")
        return

    pan_core_root = os.path.join(output_dir, "pan_core_block")
    os.makedirs(pan_core_root, exist_ok=True)

    # 合并全基因组坐标
    for name in sorted(list(names_set)):
        merged = os.path.join(pan_core_root, f"RP_qryRecord.{name}")
        if resume and os.path.isfile(merged) and utils.file_size(merged) > 0:
            logger.info("[resume] 跳过已有合并坐标: %s", merged)
        else:
            with open(merged, "w") as out:
                for chr_name in chr_list:
                    # 从coreB目录读取: genomeFilter_ChrXX/coreB/ChrXX.Os.RP_qryRecord.{sample}
                    path = os.path.join(output_dir, f"genomeFilter_{chr_name}", "coreB", f"{chr_name}.Os.RP_qryRecord.{name}")
                    if os.path.isfile(path):
                        with open(path) as fi:
                            for line in fi:
                                out.write(line)
                    else:
                        logger.debug("[%s] 样本 %s 的染色体级别结果不存在: %s", chr_name, name, path)
            logger.info("合并全基因组坐标: %s 行=%d 大小=%dB", 
                       merged, utils.count_lines(merged), utils.file_size(merged))

    # 提取序列
    fasta_dir = os.path.join(pan_core_root, "fasta")
    os.makedirs(fasta_dir, exist_ok=True)
    for name in sorted(list(names_set)):
        genome_fa = os.path.join(genome_dir, f"{name}.fasta")
        bed_path = os.path.join(pan_core_root, f"RP_qryRecord.{name}")
        out_fa = os.path.join(fasta_dir, f"RP_qryRecord.{name}.fasta")
        # 使用 -name 参数保留BED文件第一列的block编号信息
        cmd = f"bedtools getfasta -fi {genome_fa} -bed {bed_path} -fo {out_fa} -name"
        if resume and os.path.isfile(out_fa) and utils.file_size(out_fa) > 0:
            logger.info("[resume] 跳过已有 fasta: %s", out_fa)
        else:
            utils.run_cmd(cmd, logger)
        if utils.file_size(out_fa) == 0:
            logger.warning("bedtools 输出为空: %s", out_fa)
        else:
            logger.info("序列提取完成: %s 大小=%dB", out_fa, utils.file_size(out_fa))


def has_sequence_outputs(output_dir: str, chr_list: List[str]) -> bool:
    """
    判断序列提取步骤是否完成：
    1. 检查染色体级别结果是否存在（genomeFilter_ChrXX/coreB/ChrXX.Os.RP_qryRecord.*）
    2. 检查基因组级别合并坐标是否存在（pan_core_block/RP_qryRecord.*）
    3. 检查序列文件是否存在（pan_core_block/fasta/RP_qryRecord.*.fasta）
    """
    # 检查染色体级别结果（从coreB目录）
    has_chr_outputs = False
    for chr_name in chr_list:
        coreB_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}", "coreB")
        if os.path.isdir(coreB_dir):
            chr_files = [
                f for f in os.listdir(coreB_dir) if f.startswith(f"{chr_name}.Os.RP_qryRecord.")
            ]
            if chr_files:
                has_chr_outputs = True
                break
    
    if not has_chr_outputs:
        return False

    # 检查基因组级别合并坐标
    pan_core_root = os.path.join(output_dir, "pan_core_block")
    if not os.path.isdir(pan_core_root):
        return False
    merged_files = [
        f
        for f in os.listdir(pan_core_root)
        if f.startswith("RP_qryRecord.") and not os.path.isdir(os.path.join(pan_core_root, f))
    ]
    if not merged_files:
        return False
    
    # 检查序列文件
    fasta_dir = os.path.join(pan_core_root, "fasta")
    if not os.path.isdir(fasta_dir):
        return False
    fasta_files = [
        f for f in os.listdir(fasta_dir) if f.endswith(".fasta")
    ]
    return len(fasta_files) > 0 and any(utils.file_size(os.path.join(fasta_dir, f)) > 0 for f in fasta_files)