from Bio.Seq import Seq
import re
from Bio.SeqUtils import gc_fraction
"""DNA序列分析器（增强版）"""
class DNAAnalyzer:
    def __init__(self, sequence):
        self.seq = Seq(sequence.upper())
        self.length = len(self.seq)
    """获取基本统计信息（长度、GC、碱基数目）"""
    def get_info(self):
        return {
            'length': self.length,
            'GC值': round(gc_fraction(self.seq) * 100, 2),
            'A_count': self.seq.count('A'),
            'T_count': self.seq.count('T'),
            'G_count': self.seq.count('G'),
            'C_count': self.seq.count('C')
        }
    """获取反向互补序列"""
    def get_reverse_complement(self):
        rev_seq = self.seq.reverse_complement()
        return str(rev_seq)
    """翻译为蛋白质（修复警告：自动修剪不完整密码子）"""
    def translate(self, table=1):
        # ⚠️ 修复警告：将序列修剪为3的整数倍，去除末尾不完整密码子
        trimmed_seq = self.seq[:len(self.seq) // 3 * 3]
        return str(trimmed_seq.translate(table=table, to_stop=True))
    """六框翻译（修复警告：自动修剪不完整密码子）"""
    def six_frame_translation(self):
        results = {}
        # 正链3个阅读框
        for frame in range(3):
            seq_fragment = self.seq[frame:]
            # ⚠️ 修复警告：修剪片段为3的整数倍
            trimmed_frag = seq_fragment[:len(seq_fragment) // 3 * 3]
            results[f'+{frame + 1}'] = str(trimmed_frag.translate(to_stop=True))

        rev_seq = self.seq.reverse_complement()
        # 反链3个阅读框
        for frame in range(3):
            seq_fragment = rev_seq[frame:]
            # ⚠️ 修复警告：修剪片段为3的整数倍
            trimmed_frag = seq_fragment[:len(seq_fragment) // 3 * 3]
            results[f'-{frame + 1}'] = str(trimmed_frag.translate(to_stop=True))
        return results

    """寻找最长阅读框（修复警告：自动修剪不完整密码子）"""

    def find_orfs(self, min_length):
        orfs = []
        # 正链3个阅读框
        for frame in range(3):
            seq_fragment = self.seq[frame:]
            # ⚠️ 修复警告：修剪片段为3的整数倍
            trimmed_frag = seq_fragment[:len(seq_fragment) // 3 * 3]
            trans = str(trimmed_frag.translate(table=1, to_stop=True))
            for match in re.finditer(r'M[^*]*', trans):
                if len(match.group()) >= min_length:
                    orfs.append({
                        "strand": "+",
                        "frame": frame + 1,
                        "protein": match.group(),
                        'length_aa': len(match.group())
                    })
        rev_seq = self.seq.reverse_complement()
        # 反链3个阅读框
        for frame in range(3):
            seq_fragment = rev_seq[frame:]
            # ⚠️ 修复警告：修剪片段为3的整数倍
            trimmed_frag = seq_fragment[:len(seq_fragment) // 3 * 3]
            trans = trimmed_frag.translate(table=1, to_stop=True)
            for match in re.finditer(r'M[^*]*', str(trans)):
                if len(match.group()) >= min_length:
                    orfs.append({
                        "strand": "-",
                        "frame": frame + 1,
                        "protein": match.group(),
                        "length_aa": len(match.group())
                    })
        return orfs
# 测试代码
analyzer = DNAAnalyzer("ATGGCCCTGTGGATGCGCCTCCTGCCCGTGCTGTAAGCCGCCGCC")
orfs1 = analyzer.find_orfs(10)
for orf in orfs1:
    print(f"链: {orf['strand']}, 阅读框: {orf['frame']}, 长度: {orf['length_aa']}aa")
    print(f"蛋白质: {orf['protein'][:50]}\n")
print(analyzer.get_info())
print(f"反向互补: {analyzer.get_reverse_complement()}")
print(f"标准翻译: {analyzer.translate()}")
print(f"六框翻译: {analyzer.six_frame_translation()}")