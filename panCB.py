import argparse
import logging
import os
import shutil
import sys
from datetime import datetime
from multiprocessing import Pool
from typing import Any, Dict, List, Tuple

# Ensure this panCB directory is prioritized for imports
sys.path.insert(0, os.path.dirname(__file__))

import genome_alignment
import core_blocks
import sequence_process
import utils
import report
import alternative_block


def ensure_dir(path: str) -> None:
    """
    确保目录存在，如果不存在则创建。
    
    Args:
        path: 目录路径
    """
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


def cleanup_tmp_directories(output_dir: str, chr_list: List[str]) -> None:
    """
    清理所有染色体的临时目录：output/genomeFilter_ChrXX/tmp
    """
    logger = logging.getLogger("panCB.cleanup")
    cleaned_count = 0
    failed_count = 0
    
    for chr_name in chr_list:
        tmp_dir = os.path.join(output_dir, f"genomeFilter_{chr_name}", "tmp")
        if os.path.isdir(tmp_dir):
            try:
                shutil.rmtree(tmp_dir)
                logger.info("[%s] 已删除临时目录: %s", chr_name, tmp_dir)
                cleaned_count += 1
            except Exception as e:
                logger.warning("[%s] 删除临时目录失败: %s, 错误: %s", chr_name, tmp_dir, e)
                failed_count += 1
    
    if cleaned_count > 0:
        logger.info("临时目录清理完成: 成功删除 %d 个目录", cleaned_count)
    if failed_count > 0:
        logger.warning("临时目录清理部分失败: %d 个目录删除失败", failed_count)


def setup_logging(output_dir: str) -> None:
    """
    设置日志系统
    
    Args:
        output_dir: 输出目录
    """
    ensure_dir(output_dir)
    logs_dir = os.path.join(output_dir, "logs")
    ensure_dir(logs_dir)
    log_file = os.path.join(logs_dir, f"run_panCB_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")

    root = logging.getLogger()
    root.setLevel(logging.INFO)

    # 文件日志：固定为INFO级别，减少日志文件大小
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    root.addHandler(fh)

    # 控制台日志：固定为INFO级别
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    root.addHandler(sh)


def _detect_chr_iter_status(output_dir: str, chr_list: List[str], min_len: int, ref_genome: str) -> Tuple[List[str], List[str]]:
    """
    检测每条染色体的完整状态，返回待处理和已完成的染色体列表。
    
    判断条件（层级检测）：
    1. 迭代必须完全完成（达到预期轮次）
    2. 趋势文件必须存在（TSV和PNG）
    3. 染色体级别坐标文件必须存在
    
    如果某个层级的文件缺失但上一层级存在，会尝试重新生成。
    
    Args:
        output_dir: 输出目录
        chr_list: 染色体列表
        min_len: 最小长度阈值
        ref_genome: 参考基因组名称
    
    Returns:
        (pending_chr, done_chr): 待处理染色体列表和已完成染色体列表
    
    _Requirements: 1.1, 2.1, 4.2, 4.3_
    """
    logger = logging.getLogger("panCB.resume")
    pending: List[str] = []
    done: List[str] = []
    
    for chr_name in chr_list:
        logger.debug("[%s] 开始层级检测...", chr_name)
        
        # Level 1: 检查迭代是否完全完成（达到预期轮次）
        iter_complete = core_blocks.is_iteration_complete(output_dir, chr_name, min_len, ref_genome)
        
        if not iter_complete:
            # 迭代未完成，记录详细状态
            expected = core_blocks.get_expected_iter_count(output_dir, chr_name, ref_genome)
            last_completed = core_blocks.get_last_completed_iter(output_dir, chr_name, min_len)
            
            # 区分三种状态：迭代未开始、迭代未完成、迭代状态异常
            if last_completed == 0:
                logger.info("[%s] 迭代未开始（仅有比对文件），需要执行 [已完成: R%d, 预期: R%d]", 
                           chr_name, last_completed, expected if expected > 0 else 0)
            elif expected > 0 and last_completed > 0:
                logger.info("[%s] 迭代未完成: 已完成R%d, 预期R%d，需要继续执行", 
                           chr_name, last_completed, expected)
            else:
                logger.info("[%s] 迭代状态异常（无法确定预期轮次），需要重新执行 [已完成: R%d, 预期: 未知]", 
                           chr_name, last_completed)
            pending.append(chr_name)
            continue
        
        # Level 2: 检查趋势文件（TSV和PNG）
        has_trend = core_blocks.has_chr_trend_files(output_dir, chr_name)
        
        if not has_trend:
            # 趋势文件缺失，尝试从迭代结果重新生成
            logger.info("[%s] 趋势文件缺失，尝试从迭代结果重新生成...", chr_name)
            if core_blocks.ensure_chr_outputs_complete(output_dir, chr_name, min_len, ref_genome):
                logger.info("[%s] 趋势文件重新生成成功", chr_name)
                has_trend = True
            else:
                logger.warning("[%s] 趋势文件重新生成失败", chr_name)
        
        # Level 3: 检查染色体级别坐标文件
        has_coords = core_blocks.has_chr_level_coords(output_dir, chr_name)
        
        if not has_coords:
            # 坐标文件缺失，检查是否可以重新生成
            logger.info("[%s] 染色体级别坐标文件缺失，检查是否可以重新生成...", chr_name)
            if sequence_process.can_regenerate_chr_coords(output_dir, chr_name, min_len):
                # 尝试重新生成坐标文件
                logger.info("[%s] 尝试从迭代结果重新生成坐标文件...", chr_name)
                if sequence_process.extract_chr_level_coords(
                    output_dir=output_dir,
                    chr_name=chr_name,
                    ref_genome=ref_genome,
                    resume=False,
                    force_regenerate=True,
                ):
                    logger.info("[%s] 坐标文件重新生成成功", chr_name)
                    has_coords = True
                else:
                    logger.warning("[%s] 坐标文件重新生成失败", chr_name)
            else:
                logger.warning("[%s] 无法重新生成坐标文件（迭代结果不可用）", chr_name)
        
        # 综合判断：迭代完成且趋势文件和坐标文件都存在
        if iter_complete and has_trend and has_coords:
            logger.info("[%s] 检测完成: 迭代=完成, 趋势=存在, 坐标=存在 -> 已完成", chr_name)
            done.append(chr_name)
        else:
            # 记录缺失的内容
            missing = []
            if not has_trend:
                missing.append("趋势文件")
            if not has_coords:
                missing.append("坐标文件")
            logger.info("[%s] 检测完成: 迭代=完成, 缺失=[%s] -> 需要重新处理", 
                       chr_name, ", ".join(missing))
            pending.append(chr_name)
    
    # 汇总日志
    logger.info("层级检测完成: 已完成=%d条染色体, 待处理=%d条染色体", len(done), len(pending))
    if done:
        logger.debug("已完成染色体: %s", ", ".join(done))
    if pending:
        logger.debug("待处理染色体: %s", ", ".join(pending))
    
    return pending, done


def _run_post_alignment_task(task: Dict[str, Any]) -> Dict[str, Any]:
    """
    执行单个染色体的后续处理任务（步骤3-4.5）。
    
    此函数在染色体级别的进程池中并行执行。每个染色体的处理包括：
    - 步骤3：排序与优先级
    - 步骤4：迭代识别核心block
    - 步骤4.5：生成染色体级别核心block坐标
    
    注意：步骤2（比对）已在阶段1完成，此函数处理后续步骤。
    
    Args:
        task: 任务字典，包含染色体名称、filtered_files等
    
    Returns:
        包含chr_name、success、message的字典
    """
    chr_name = task["chr_name"]
    filtered_files = task["filtered_files"]
    logger = logging.getLogger(f"panCB.core.{chr_name}")
    
    try:
        logger.info("[%s] 步骤3：排序与优先级", chr_name)
        acc_order = core_blocks.generate_alignment_priority_for_chr(
            chr_name=chr_name,
            filtered_1227_files=filtered_files,
            output_dir=task["output_dir"],
            ref_genome=task["ref_name"],
        )

        logger.info("[%s] 步骤4：迭代识别核心block", chr_name)
        final_iter_file = core_blocks.run_iterative_core_blocks_for_chr(
            genome_dir=task["genome_dir"],
            output_dir=task["output_dir"],
            chr_name=chr_name,
            ref_genome=task["ref_name"],
            acc_order=acc_order,
            resume=task["resume"],
            blastn_threads=task["blastn_threads"],
            debug=task["debug"],
            core_coverage=task.get("core_coverage", 80.0),
            fast_iteration=task.get("fast_iteration", False),
        )

        # 硬编码参数：核心块最小长度阈值
        min_len = 50
        success = core_blocks.has_final_iter_result(task["output_dir"], chr_name, min_len)
        if success:
            # 步骤4.5：立即生成染色体级别的核心block坐标结果
            logger.info("[%s] 步骤4.5：生成染色体级别核心block坐标", chr_name)
            chr_coords_success = sequence_process.extract_chr_level_coords(
                output_dir=task["output_dir"],
                chr_name=chr_name,
                ref_genome=task["ref_name"],
                resume=task["resume"],
            )
            if chr_coords_success:
                logger.info("[%s] 染色体级别核心block坐标生成完成", chr_name)
            else:
                logger.warning("[%s] 染色体级别核心block坐标生成失败", chr_name)
        
        message = "" if success else "未检测到迭代结果"
        return {"chr_name": chr_name, "success": success, "message": message}
    except Exception as e:
        logger.exception("[%s] 后续处理异常: %s", chr_name, e)
        return {"chr_name": chr_name, "success": False, "message": str(e)}


def parse_args():
    """
    解析命令行参数。
    
    Returns:
        解析后的参数对象
    """
    p = argparse.ArgumentParser(description="panCB: Pan-genome core block identification pipeline with chromosome-based alignment and merging")
    p.add_argument("-r", "--ref", required=True, help="Reference genome name (without file extension)")
    p.add_argument("-g", "--genome_dir", required=True, help="Genome directory containing {name}.fasta files")
    p.add_argument("-o", "--output", required=True, help="Output directory")
    p.add_argument("-accs", required=True, help="Genome list file, one genome name per line")

    p.add_argument("--processes-nucmer", type=int, default=5, 
                   help="Number of parallel processes for nucmer alignment phase (default: 5). "
                        "Controls how many nucmer alignment tasks run simultaneously.")
    p.add_argument("--processes-core", type=int, default=5, 
                   help="Number of parallel processes for core block iteration phase (default: 5). "
                        "Controls how many chromosomes are processed in parallel during iteration.")
    p.add_argument("--resume", action="store_true", help="Enable resume mode, skip existing outputs")
    p.add_argument("--step", choices=["all", "core"], default="all", 
                   help="Execution steps: all=core+alternative, core=core blocks only (default: all)")
    
    # Tool configuration parameters
    p.add_argument("--blastn-threads", type=int, default=1, help="Number of threads for blastn (default: 1)")
    p.add_argument("--debug", action="store_true", 
                   help="Enable debug mode, output detailed block filtering logs to output/logs/ChrXX_debug.log")
    
    # Core block filtering parameters
    p.add_argument("--core-identity", type=float, default=98.0,
                   help="Identity threshold for core block filtering (default: 98.0, range: 0-100)")
    p.add_argument("--core-coverage", type=float, default=80.0,
                   help="Coverage threshold for core block filtering (default: 80.0, range: 0-100)")
    
    # Iteration optimization parameters
    p.add_argument("--fast", action="store_true",
                   help="Enable fast iteration mode with delayed BLAST validation. "
                        "Only the final iteration performs BLAST verification with flanking sequences. "
                        "When not specified, standard per-iteration BLAST validation is used.")
    p.add_argument("--round-verify", type=int, default=0,
                   help="[DEPRECATED] This parameter is deprecated and will be ignored. "
                        "Fast mode now uses a fixed verification interval of every 2 rounds. "
                        "Kept for backward compatibility only.")
    p.add_argument("--flanking-size", type=int, default=0,
                   help="[DEPRECATED] This parameter is deprecated and will be ignored in fast mode. "
                        "Fast mode now uses dynamic flanking (0.5%% of block length). "
                        "Kept for backward compatibility only.")
    
    return p.parse_args()


def main():
    """
    panCB主流程入口函数。
    执行完整的泛基因组核心block鉴定流程。
    """
    args = parse_args()
    setup_logging(args.output)
    logging.info("参数: %s", vars(args))
    logging.info("核心block过滤参数: 身份度阈值=%.1f%%, 覆盖度阈值=%.1f%%", 
                 args.core_identity, args.core_coverage)
    
    # 处理弃用参数警告
    if args.round_verify > 0:
        logging.warning("[DEPRECATED] --round-verify参数已弃用，快速模式现在固定每2轮验证一次。"
                       "提供的值 %d 将被忽略。", args.round_verify)
    if args.flanking_size > 0:
        logging.warning("[DEPRECATED] --flanking-size参数在快速模式下已弃用，"
                       "现在使用动态flanking（block长度的0.5%%）。提供的值 %d 将被忽略。", 
                       args.flanking_size)
    
    if args.fast:
        logging.info("使用快速迭代模式 + 增量验证（固定每2轮）+ 动态flanking")
    else:
        logging.info("使用标准迭代模式（每轮BLAST验证）")

    resume_mgr = utils.ResumeManager(args.output)
    if args.resume:
        logging.info("[resume] %s", resume_mgr.describe_progress())
    else:
        resume_mgr.reset()

    # 预生成 .fai 索引
    try:
        names = utils.read_accs(args.accs)
        if args.ref not in names:
            names.insert(0, args.ref)
        for n in names:
            utils.ensure_fai_exists(os.path.join(args.genome_dir, f"{n}.fasta"), logging.getLogger("panCB.index"))
    except Exception as e:
        logging.warning("预生成 .fai 索引时出现异常（忽略继续）: %s", e)

    # 1) 按染色体拆分
    logging.info("步骤1：按染色体拆分开始")
    if args.resume and resume_mgr.is_step_done("split"):
        logging.info("[resume] 步骤1 已完成，直接加载历史染色体信息")
        chr_list_meta = resume_mgr.get_step_meta("split")
        chr_list: List[str] = chr_list_meta.get("chr_list", [])
        if not chr_list:
            chr_list = genome_alignment.split_all_genomes_by_chr(
                genome_dir=args.genome_dir,
                accs_file=args.accs,
                output_dir=args.output,
                ref_name=args.ref,
                resume=True,
            )
    else:
        resume_mgr.mark_step_started("split")
        chr_list = genome_alignment.split_all_genomes_by_chr(
            genome_dir=args.genome_dir,
            accs_file=args.accs,
            output_dir=args.output,
            ref_name=args.ref,
            resume=args.resume,
        )
        resume_mgr.mark_step_done("split", {"chr_list": chr_list})
    logging.info("步骤1完成：检测到染色体 %s", ",".join(chr_list))

    split_root = os.path.join(args.output, "Chr_split")

    # 2) 逐染色体：比对→初筛→排序→迭代（仅核心block）
    pending_chr: List[str] = []
    done_chr: List[str] = []
    need_align = True
    if args.resume:
        meta = resume_mgr.get_step_meta("align_core")
        pending_chr = list(meta.get("pending_chr", []))
        done_chr = list(meta.get("done_chr", []))
        step_status = resume_mgr.get_step_status("align_core")

        if resume_mgr.is_step_done("align_core"):
            if not pending_chr:
                # 硬编码参数：核心块最小长度阈值
                min_len = 50
                pending_chr, done_chr = _detect_chr_iter_status(args.output, chr_list, min_len, args.ref)
                if not pending_chr:
                    logging.info("[resume] 步骤2-4 已完成，跳过核心比对与迭代")
                    need_align = False
                else:
                    logging.info(
                        "[resume] 检测到 %d 条染色体缺少核心结果，重新执行：%s",
                        len(pending_chr),
                        ",".join(pending_chr),
                    )
                    resume_mgr.reopen_step("align_core", {"pending_chr": pending_chr, "done_chr": done_chr})
            else:
                resume_mgr.reopen_step("align_core", {"pending_chr": pending_chr, "done_chr": done_chr})
        else:
            if not pending_chr:
                # 硬编码参数：核心块最小长度阈值
                min_len = 50
                pending_chr, detected_done = _detect_chr_iter_status(args.output, chr_list, min_len, args.ref)
                done_chr = detected_done
            if not pending_chr:
                logging.info("[resume] 检测到所有染色体已有核心迭代结果，标记 align_core 已完成")
                resume_mgr.mark_step_done("align_core", {"pending_chr": [], "done_chr": done_chr})
                need_align = False
            elif step_status != "running":
                resume_mgr.mark_step_started("align_core", {"pending_chr": pending_chr, "done_chr": done_chr})
    else:
        pending_chr = list(chr_list)
        done_chr = []
        resume_mgr.mark_step_started("align_core", {"pending_chr": pending_chr, "done_chr": done_chr})

    if need_align:
        target_chrs = set(pending_chr) if pending_chr else set(chr_list)
        
        # ========== 阶段1：并行比对所有样本×染色体（任务级并行）==========
        logging.info("\n" + "=" * 80)
        logging.info("开始两阶段并行处理策略")
        logging.info("=" * 80)
        
        # 调用新的任务级并行函数，一次性完成所有比对
        chr_to_filtered_files = genome_alignment.run_all_alignments_parallel(
            ref_name=args.ref,
            chr_list=chr_list,
            split_root=split_root,
            output_dir=args.output,
            accs_file=args.accs,
            processes=args.processes_nucmer,
            resume=args.resume,
            debug=args.debug,
            core_identity=args.core_identity,
        )
        
        # ========== 阶段2：并行处理后续步骤（染色体级并行）==========
        logging.info("\n" + "=" * 80)
        logging.info("阶段2：并行处理排序、迭代、提取（染色体级并行）")
        logging.info("=" * 80)
        
        # 构建后续处理任务
        post_tasks: List[Dict[str, Any]] = []
        for chr_name in chr_list:
            if target_chrs and chr_name not in target_chrs:
                logging.info("[%s] 已完成，跳过后续处理", chr_name)
                continue
            
            filtered_files = chr_to_filtered_files.get(chr_name, [])
            if not filtered_files:
                logging.warning("[%s] 没有比对结果，跳过后续处理", chr_name)
                continue
            
            post_tasks.append(
                {
                    "chr_name": chr_name,
                    "ref_name": args.ref,
                    "output_dir": args.output,
                    "genome_dir": args.genome_dir,
                    "filtered_files": filtered_files,
                    "resume": args.resume,
                    "blastn_threads": args.blastn_threads,
                    "debug": args.debug,
                    "core_identity": args.core_identity,
                    "core_coverage": args.core_coverage,
                    "fast_iteration": args.fast,
                }
            )
        
        if not post_tasks:
            logging.info("没有需要后续处理的染色体，跳过步骤3-4")
        else:
            def handle_post_result(res: Dict[str, Any]) -> None:
                chr_name = res.get("chr_name")
                if res.get("success"):
                    if chr_name in pending_chr:
                        pending_chr.remove(chr_name)
                    if chr_name not in done_chr:
                        done_chr.append(chr_name)
                    resume_mgr.update_step_meta(
                        "align_core",
                        {
                            "pending_chr": pending_chr,
                            "done_chr": done_chr,
                        },
                    )
                    logging.info("[%s] 后续处理完成", chr_name)
                else:
                    logging.warning("[%s] 后续处理失败：%s", chr_name, res.get("message"))
            
            # 染色体级并行处理后续步骤
            requested_pool_size = args.processes_core if args.processes_core and args.processes_core > 1 else 1
            actual_pool_size = min(requested_pool_size, len(post_tasks)) if requested_pool_size > 1 else 1
            
            if actual_pool_size > 1 and len(post_tasks) > 1:
                if actual_pool_size < requested_pool_size:
                    logging.info("开始并行处理 %d 条染色体的后续步骤，进程池大小: %d (请求: %d, 限制: 任务数)", 
                                 len(post_tasks), actual_pool_size, requested_pool_size)
                else:
                    logging.info("开始并行处理 %d 条染色体的后续步骤，进程池大小: %d", 
                                 len(post_tasks), actual_pool_size)
                with Pool(actual_pool_size) as pool:
                    for res in pool.imap_unordered(_run_post_alignment_task, post_tasks):
                        handle_post_result(res)
            else:
                logging.info("顺序处理 %d 条染色体的后续步骤", len(post_tasks))
                for task in post_tasks:
                    res = _run_post_alignment_task(task)
                    handle_post_result(res)
        
        logging.info("\n" + "=" * 80)
        logging.info("两阶段处理完成")
        logging.info("=" * 80 + "\n")

        if not pending_chr:
            resume_mgr.mark_step_done("align_core", {"pending_chr": [], "done_chr": done_chr or chr_list})

    # 3) 合并所有染色体坐标并提取序列（输出到 pan_core_block）
    logging.info("步骤6：合并全基因组核心坐标并提取序列")
    seq_ready = sequence_process.has_sequence_outputs(args.output, chr_list)
    if args.resume and resume_mgr.is_step_done("seq_extract") and seq_ready:
        logging.info("[resume] 步骤5 已完成，跳过序列提取")
    else:
        if args.resume and resume_mgr.is_step_done("seq_extract") and not seq_ready:
            logging.info("[resume] 检测到序列输出缺失，重新执行步骤5")
        resume_mgr.mark_step_started("seq_extract")
        sequence_process.merge_chr_coords_and_extract_sequences(
            genome_dir=args.genome_dir,
            output_dir=args.output,
            chr_list=chr_list,
            resume=args.resume,
        )
        if sequence_process.has_sequence_outputs(args.output, chr_list):
            resume_mgr.mark_step_done("seq_extract")
        else:
            logging.warning("步骤5执行后仍未检测到核心序列输出，请检查日志")

    # 3.5) 可选：鉴定 alternative_block（按染色体输出到 alterB）
    if args.step == "all":
        logging.info("步骤6.5：按染色体计算 alternative_block")
        if args.resume and resume_mgr.is_step_done("alternative"):
            logging.info("[resume] 步骤6.5 已完成，跳过 alternative 计算")
        else:
            resume_mgr.mark_step_started("alternative")
            alternative_block.run_alternative_for_all_chr(
                genome_dir=args.genome_dir,
                output_dir=args.output,
                chr_list=chr_list,
                ref_genome=args.ref,
                resume=args.resume,
            )
            resume_mgr.mark_step_done("alternative")
    else:
        resume_mgr.mark_step_skipped("alternative", "参数 step=core，未执行 alternative")

    # 4) 统计报告（TSV格式）
    summary_ready = report.has_summary_outputs(args.output)
    if args.resume and resume_mgr.is_step_done("report_summary") and summary_ready:
        logging.info("[resume] 步骤7 核心长度统计已完成，跳过")
    else:
        if args.resume and resume_mgr.is_step_done("report_summary") and not summary_ready:
            logging.info("[resume] 检测到统计报告缺失，重新生成")
        logging.info("步骤7：核心长度统计")
        resume_mgr.mark_step_started("report_summary")
        report.summarize_core_lengths(
            output_dir=args.output,
            chr_list=chr_list,
            ref_genome=args.ref,
            resume=args.resume,
        )
        if report.has_summary_outputs(args.output):
            resume_mgr.mark_step_done("report_summary")
        else:
            logging.warning("步骤7执行后仍未检测到报告汇总文件，请检查")

    # 5) 详细文本报告和可视化报告
    detail_ready = report.has_detail_report(args.output)
    html_report_path = os.path.join(args.output, "pan_core_block", "report", "core_summary.html")
    html_ready = os.path.isfile(html_report_path) and utils.file_size(html_report_path) > 0
    
    if args.resume and resume_mgr.is_step_done("report_detail") and detail_ready and html_ready:
        logging.info("[resume] 步骤8 详细报告和可视化报告已完成，跳过")
    else:
        if args.resume and resume_mgr.is_step_done("report_detail") and (not detail_ready or not html_ready):
            logging.info("[resume] 检测到详细报告或可视化报告缺失，重新生成")
        logging.info("步骤8：生成详细文本报告和可视化报告")
        resume_mgr.mark_step_started("report_detail")
        
        # 生成TXT详细报告
        report.write_detailed_summary(
            output_dir=args.output,
            chr_list=chr_list,
            ref_genome=args.ref,
            resume=args.resume,
        )
        
        # 生成HTML可视化报告
        report.generate_all_visualizations(
            output_dir=args.output,
            chr_list=chr_list,
            ref_genome=args.ref,
            resume=args.resume,
        )
        
        if report.has_detail_report(args.output):
            resume_mgr.mark_step_done("report_detail")
        else:
            logging.warning("步骤8执行后仍未检测到详细报告，请检查")

    # 6) 生成基因组级别核心block变化趋势HTML报告
    # _Requirements: 3.1, 3.2, 3.3_
    trend_report_path = os.path.join(args.output, "pan_core_block", "report", "core_block_trend_report.html")
    
    # 检查HTML报告是否存在
    html_report_exists = os.path.isfile(trend_report_path) and utils.file_size(trend_report_path) > 0
    
    if args.resume and html_report_exists:
        logging.info("[resume] 步骤9 迭代趋势报告已存在，跳过")
    else:
        # 检查是否有足够的趋势文件来生成报告
        trend_files_available = 0
        for chr_name in chr_list:
            if core_blocks.has_chr_trend_files(args.output, chr_name):
                trend_files_available += 1
        
        if trend_files_available == 0:
            logging.warning("步骤9：没有可用的染色体趋势文件，无法生成基因组级别报告")
        else:
            if args.resume and not html_report_exists:
                logging.info("[resume] 步骤9：检测到HTML报告缺失但趋势文件存在（%d/%d条染色体），重新生成报告...", 
                           trend_files_available, len(chr_list))
            else:
                logging.info("步骤9：生成基因组级别核心block变化趋势报告（使用%d条染色体的趋势数据）", 
                           trend_files_available)
            
            # 生成报告
            if core_blocks.generate_genome_trend_report(args.output, chr_list):
                logging.info("步骤9：基因组级别HTML报告生成成功: %s", trend_report_path)
            else:
                logging.warning("步骤9：基因组级别HTML报告生成失败")

    logging.info("panCB流程完成。输出根目录：%s", args.output)
    
    # 清理临时目录
    cleanup_tmp_directories(args.output, chr_list)


if __name__ == "__main__":
    main()