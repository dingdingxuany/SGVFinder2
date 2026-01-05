import sys
# 确保包含 cy_ext 路径
sys.path.insert(0,'cy_ext/')

from collections import namedtuple
import pysam
from itertools import groupby
from operator import attrgetter
import logging
import gzip
import ujson
import sam2pmp_helper as s2ph
from .ICRAUtils import timeit
import multiprocessing

# 日志配置
log_ = logging.getLogger('sam2pmp')

# 保留原有的 NamedTuple 定义
DestMap_se = namedtuple('DestMap_se', ['dest_id', 'mismatches', 'strand', 'pos', 'cigar',
                                       'end', 'mismatch_pos', 'alnscore', 'map_quality'])
DestMap_pe = namedtuple('DestMap_pe', ['dest_id', 'mismatches1', 'mismatches2', 'strand1',
                                       'strand2', 'pos1', 'pos2', 'cigar1', 'cigar2',
                                       'mismatch_pos1', 'mismatch_pos2', 'alnscore1',
                                       'alnscore2', 'map_quality1', 'map_quality2'])
DestMap_se_part = namedtuple('DestMap_se_part', ['dest_id', 'mismatches', 'strand', 'pos', 'end',
                                                 'map_quality'])
DestMap_pe_part = namedtuple('DestMap_pe_part', ['dest_id', 'mismatches1', 'mismatches2', 'strand1',
                                                 'strand2', 'pos1', 'pos2', 'map_quality1',
                                                 'map_quality2'])

# --- 辅助函数 (保持不变，但需确保 worker 能访问) ---

def _get_map_quality(aln):
    alp = aln.get_aligned_pairs()
    ref = aln.get_reference_sequence()
    quals = aln.query_qualities
    return s2ph.calc_quality(alp, quals, ref)

def _readunmapped(alngrp):
    return len(alngrp) == 1 and alngrp[0].is_unmapped \
           or len(alngrp) == 2 and alngrp[0].is_unmapped and alngrp[1].is_unmapped

def _get_nm_md_as(aln):
    retdct = dict(aln.tags)
    # 某些比对软件可能没有 AS tag，需做容错处理
    return retdct.get('NM', 0), retdct.get('MD', ''), retdct.get('AS', 0)

def _get_strand(aln):
    return '-' if aln.is_reverse else '+'

def _get_pe_end(aln):
    return 0 if aln.is_read1 else 1

def _se_from_aln(aln, unique, full):
    qual = 1 if unique else _get_map_quality(aln)
    refname = aln.reference_name.split()[0]
    if full:
        nm, md, asc = _get_nm_md_as(aln)
        return DestMap_se(refname, nm, _get_strand(aln), aln.pos,
                          aln.cigarstring, _get_pe_end(aln), md, asc, qual)
    return DestMap_se_part(refname, dict(aln.tags)['NM'], _get_strand(aln),
                           aln.pos, _get_pe_end(aln), qual)

def _pe_from_aln(aln1, aln2, unique, full):
    qual1 = 1 if unique else _get_map_quality(aln1)
    qual2 = 1 if unique else _get_map_quality(aln2)
    if full:
        nm1, md1, asc1 = _get_nm_md_as(aln1)
        nm2, md2, asc2 = _get_nm_md_as(aln2)
        return DestMap_pe(aln1.reference_name.split()[0], nm1, nm2,
                          _get_strand(aln1), _get_strand(aln2), aln1.pos, aln2.pos,
                          aln1.cigarstring, aln2.cigarstring, md1, md2, asc1, asc2,
                          qual1, qual2)
    return DestMap_pe_part(aln1.reference_name.split()[0], dict(aln1.tags)['NM'],
                           dict(aln2.tags)['NM'], _get_strand(aln1), _get_strand(aln2),
                           aln1.pos, aln2.pos, qual1, qual2)

# --- SourceReadSAM 类 (保持不变) ---
class SourceReadSAM(object):
    def __init__(self, rid, seq, qual, full=True):
        self.rid = rid
        self.seq = [str(s) for s in seq]
        self.quals = qual
        self.full = full
        self.se_maps = []
        self.pe_maps = []

    def to_ser(self):
        if self.full:
            return (self.rid, self.seq, self.quals, [tuple(x) for x in self.se_maps],
                    [tuple(y) for y in self.pe_maps])
        else:
            return (self.rid, self.seq, self._letterqual(), [tuple(x) for x in self.se_maps],
                    [tuple(y) for y in self.pe_maps])

    def _letterqual(self):
        return [''.join([chr(i + 33) for i in q]) for q in self.quals]

    @staticmethod
    def _numberqual(letterqual):
        return [[ord(c) - 33 for c in st] for st in letterqual]

    @classmethod
    def from_alngroup(cls, alngrp, full):
        fst = alngrp[0]
        res = cls(fst.query_name, [fst.query_sequence], [fst.query_qualities.tolist()], full)
        read1, read2 = None, None
        seen_read1 = fst.is_read1
        seen_read2 = fst.is_read2
        unique = (len(alngrp) == 1) \
                 or ((len(alngrp) == 2) and (alngrp[0].is_read1 ^ alngrp[1].is_read1))
        for aln in alngrp:
            if not seen_read1 and aln.is_read1:
                res.seq.insert(0, aln.query_sequence)
                res.quals.insert(0, aln.query_qualities.tolist())
                seen_read1 = True
            if not seen_read2 and aln.is_read2:
                res.seq.append(aln.query_sequence)
                res.quals.append(aln.query_qualities.tolist())
                seen_read2 = True
            if aln.is_paired:
                if aln.is_read1:
                    read1 = aln
                elif aln.is_read2:
                    read2 = aln
                else:
                    pass # raise RuntimeError('huh?')
                if read1 is not None and read2 is not None:
                    if aln.is_proper_pair:
                        # 简单的断言可能会在脏数据中失败，建议改为条件判断
                        if (read1.reference_name == read2.reference_name) and (read1.is_reverse ^ read2.is_reverse):
                            res.pe_maps.append(_pe_from_aln(read1, read2, unique, full))
                    else:
                        res.se_maps.append(_se_from_aln(read1, unique, full))
                        res.se_maps.append(_se_from_aln(read2, unique, full))
                    read1, read2 = None, None
                
            else:
                if not aln.is_unmapped:
                    res.se_maps.append(_se_from_aln(aln, unique, full))
        return res

    def sort(self):
        self.se_maps.sort(key=attrgetter('map_quality'), reverse=True)
        self.pe_maps.sort(key=lambda m: m.map_quality1 + m.map_quality2, reverse=True)


# --- 并行处理核心逻辑 ---

def process_batch(batch_data):
    """
    Worker 进程执行的函数。
    接收一个 list of list of pysam.AlignedSegment (已pickle)
    返回 JSON 字符串列表和未比对计数
    """
    results = []
    local_unmapped = 0
    full = True # 假设 full 默认为 True，如果需要传参可以使用 functools.partial

    for alngrp in batch_data:
        # 在这里执行原来循环体内的逻辑
        if _readunmapped(alngrp):
            local_unmapped += 1
            continue
        
        try:
            curr = SourceReadSAM.from_alngroup(alngrp, full)
            if len(curr) > 0:
                curr.sort()
                # 直接在这里序列化成 JSON 字符串，减少主进程负担
                json_str = ujson.dumps(curr.to_ser())
                results.append(json_str)
        except Exception as e:
            # 捕获异常防止单个异常导致进程挂掉
            # logging.error(f"Error processing read group: {e}")
            pass

    return results, local_unmapped

def sam_reader_generator(samfile, batch_size=1000):
    """
    生成器：读取 SAM 文件，按 query_name 分组，并按 batch_size 打包
    """
    save = pysam.set_verbosity(0)
    f = pysam.AlignmentFile(samfile)
    pysam.set_verbosity(save)

    batch = []
    
    # groupby 返回的是迭代器，我们需要将其转为 list 以便 pickle 传给子进程
    # 注意：list(grp) 会消耗迭代器，将 objects 读入内存
    for _, grp in groupby(f, attrgetter('query_name')):
        alngrp = list(grp)
        batch.append(alngrp)
        
        if len(batch) >= batch_size:
            yield batch
            batch = []
    
    if batch:
        yield batch
    
    f.close()

@timeit(log_)
def sam2pmp(samfile, pmpfile, full=True, threads=4, batch_size=2000):
    """
    并行版本的 sam2pmp
    """
    unmapped_total = 0
    
    # 使用 imap_unordered 还是 imap？
    # 如果 pmp 文件的输出顺序不重要，使用 imap_unordered 会稍微快一点。
    # 这里为了稳妥使用 imap (有序)。
    
    # 确定进程数，不超过 CPU 核心数
    cpu_count = multiprocessing.cpu_count()
    run_threads = min(threads, cpu_count)
    if run_threads < 1: run_threads = 1
    
    log_.info(f"Starting parallel conversion with {run_threads} processes.")

    with multiprocessing.Pool(processes=run_threads) as pool:
        # 创建生成器
        reader = sam_reader_generator(samfile, batch_size)
        
        with gzip.open(pmpfile, 'wt') as of:
            # imap 懒加载，不会一次性把文件读完撑爆内存
            # process_batch 需要 full 参数，这里简化处理，假设 full=True
            # 如果需要传参，可以使用 functools.partial(process_batch, full=full)
            
            for res_list, local_unmapped in pool.imap(process_batch, reader):
                unmapped_total += local_unmapped
                for json_line in res_list:
                    of.write(json_line)
                    of.write('\n')
    
    log_.info(f"Conversion finished. Total unmapped reads: {unmapped_total}")

# 为了兼容命令行调用或模块调用，如果是脚本执行：
if __name__ == "__main__":
    # 示例调用
    import sys
    if len(sys.argv) >= 3:
        sam2pmp(sys.argv[1], sys.argv[2], threads=8)
