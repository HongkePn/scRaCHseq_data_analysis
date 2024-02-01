# detect mutations in bam files.
import os
import sys
import pysam
import gzip
import numpy as np
from scipy.stats import hypergeom
from collections import Counter, namedtuple
os.chdir("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/script/lr/")
from parse_gene_anno import parse_gff_tree
from sc_longread import blocks_to_junctions, get_gene_blocks
import copy

def get_gene_flat(gene_to_transcript, transcript_to_exon):
    """
    Return a dictionary with gene_id as key and list of exons as value.
    Overlapping exons are merged, e.g. (1,4),(2,5) -> (1,5)
    Exons are either a list of 2 int or a tuple of 2 int.
    """
    gene_dict = {}
    for g in gene_to_transcript:
        exons = []
        for tr in gene_to_transcript[g]:
            for exon in transcript_to_exon[tr]:
                exons.append(copy.deepcopy(exon))
        if len(exons) < 2:
            gene_dict[g] = exons
        else:
            exons.sort(key=lambda exon: exon[0])
            merged_exons = [exons[0]]
            for higher in exons[1:]:  # start of higher >= start of lower
                lower = merged_exons[-1]
                if higher[0] <= lower[1]:
                    end = max(lower[1], higher[1])
                    merged_exons[-1] = (lower[0], end)  # Tuple
                else:
                    merged_exons.append(higher)
            gene_dict[g] = merged_exons
    return gene_dict

def get_fa(fn):
    ch = ""
    seq = []
    for line in open(fn):
        if line[0] == ">":
            if ch != "":
                yield ch, "".join(seq)
            ch = line[1:].strip().split()[0]
            seq = []
        else:
            seq.append(line.strip().upper())
    yield ch, "".join(seq)


def seq_entropy(seq):
    res = 0.
    for st in list(set(seq)):
        p = float(seq.count(st))/len(seq)
        res += -p*np.log(p)
    return res


def find_homo_regions(fa_seq, chr_bl, min_len=3, min_gap=1):
    """
    find regions with at least `min_len` homopolymers and
    call +/- `min_gap` in surrounding regions as homo-regions as well.
    """
    homo_dict = {}
    for bl in chr_bl:
        i=bl.s
        while i < bl.e-min_len-1:
            if fa_seq[i:(i+min_len)] == "A" * min_len:
                j= i+min_len
                while j<bl.e-1 and fa_seq[j]=="A":
                    j += 1
                for ix in range(max(0,i-min_gap),min(j+min_gap,bl.e-1)):
                    homo_dict[ix] = "A"
                i = j
            elif fa_seq[i:(i+min_len)] == "T" * min_len:
                j= i+min_len
                while j<bl.e-1 and fa_seq[j]=="T":
                    j += 1
                for ix in range(max(0,i-min_gap),min(j+min_gap,bl.e-1)):
                    homo_dict[ix] = "T"
                i = j
            elif fa_seq[i:(i+min_len)] == "G" * min_len:
                j= i+min_len
                while j<bl.e-1 and fa_seq[j]=="G":
                    j += 1
                for ix in range(max(0,i-min_gap),min(j+min_gap,bl.e-1)):
                    homo_dict[ix] = "G"
                i = j
            elif fa_seq[i:(i+min_len)] == "C" * min_len:
                j= i+min_len
                while j<bl.e-1 and fa_seq[j]=="C":
                    j += 1
                for ix in range(max(0,i-min_gap),min(j+min_gap,bl.e-1)):
                    homo_dict[ix] = "C"
                i = j
            else:
                i += 1
    return homo_dict


def update_corr_cnt(int_l, cb_corr_cnt):
    for i in range(len(int_l)-1):
        for j in range(i+1, len(int_l)):
            cb_corr_cnt[(int_l[i],int_l[j])] += 1
            cb_corr_cnt[(int_l[j],int_l[i])] += 1


def get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, out_dir, cb_seq_dict, bam_short, known_position_dict, min_cov=100, report_pct=(0.15,0.85),known_range_dict=None):
    def _ir(ch,st,en):
        for ix in known_range_dict[ch]:
            if ix[0]<=st<=ix[1] or ix[0]<=en<=ix[1] or st<=ix[0]<=en:
                return True
        return False
    c2i = {"A":0, "C":1, "G":2, "T":3}  # four array.arrays of the same length in order A C G T
    fa_dict={}
    acc_pct = []
    REF_cnt_dict = {}
    ALT_cnt_dict = {}
    cb_seq_set = set(cb_seq_dict.keys())
    reporting_summary = []
    for c in get_fa(fa_f):
        fa_dict[c[0]] = c[1]
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    if bam_short is not None:
        bam_s = pysam.AlignmentFile(bam_short, "rb")
    cb_corr_cnt = Counter()
    for ch in chr_to_blocks:
        if ch not in known_range_dict:
            continue
        #print ch
        #if ch !="chr18":
        #    continue
        #homo_dict = find_homo_regions(fa_dict[ch], chr_to_blocks[ch])
        homo_dict = {}
        for ith, bl in enumerate(chr_to_blocks[ch]):
            if not _ir(ch,bl.s,bl.e):
                continue
            print ch, bl.s, bl.e
            tmp_bl_flat = get_gene_flat({"NNN":bl.transcript_list}, transcript_to_exon)
            print(tmp_bl_flat["NNN"])
            for ex in tmp_bl_flat["NNN"]:
                #63123346-63320128
                #if ex[0]>63320128 or ex[1]<63123346:
                #    continue  # only look at BCL2
                cnt = bamfile.count(ch, ex[0], ex[1])
                if cnt < min_cov:
                    continue
                cov = bamfile.count_coverage(ch, ex[0], ex[1],
                quality_threshold=0)  # four array.arrays of the same length in order A C G T
                if len(cov[0])<20:
                    continue  # ignore tiny exons
                for i in range(2, len(cov[0])-2):  # ignore the bases at the beginning and the end (close to splicing site)
                    tot =  float(cov[0][i]+cov[1][i]+cov[2][i]+cov[3][i])
                    v_pos = ex[0]+i
                    if tot>min_cov and (fa_dict[ch][v_pos]!="N"):
                        freq = cov[c2i[fa_dict[ch][v_pos]]][i]/tot
                        acc_pct.append(freq)
                        base_freq = [("A",cov[0][i]),("C",cov[1][i]),("G",cov[2][i]),("T",cov[3][i])]
                        base_freq.sort(key=lambda x:x[1],reverse=True)
                        if v_pos == 63318364:
                            print base_freq
                        ALT = [it[0] for it in base_freq if it[0] != fa_dict[ch][v_pos]][0] # the most enriched ALT allele
                        alt_freq = cov[c2i[ALT]][i]/tot
                        if (report_pct[0]< alt_freq < report_pct[1]) or ((ch,v_pos) in known_position_dict):
                            tmp_atcg_set = {}
                            if bam_short is not None:
                                try:
                                    cov_s = bam_s.count_coverage(ch, v_pos, v_pos+1, quality_threshold=20)
                                    s_tot = cov_s[0][0]+cov_s[1][0]+cov_s[2][0]+cov_s[3][0]
                                    if s_tot> (min_cov/2):
                                        s_freq = cov_s[c2i[fa_dict[ch][v_pos]]][0]/float(s_tot)
                                    else:
                                        s_freq = -1
                                except:
                                    s_freq = -1
                            else:
                                s_freq = -1
                            seq_ent = seq_entropy(fa_dict[ch][(v_pos-10):(v_pos+10)])
                            indel_freq = -1
                            if ((ch,v_pos) in known_position_dict) or ((v_pos not in homo_dict) and (seq_ent > 1) and (s_freq==-1 or (0.05<s_freq<0.95))):
                                for pileupcolumn in bamfile.pileup(ch, v_pos, v_pos+1,truncate=True, min_base_quality=0,ignore_overlaps=False,max_depth=30000):
                                    c_keep = 0
                                    c_del = 0
                                    for pileupread in pileupcolumn.pileups:
                                        if not pileupread.is_del:
                                            if not pileupread.is_refskip:
                                                c_keep += 1
                                                cb_seq, umi_seq = pileupread.alignment.query_name.split("#")[0].split("_")
                                                if cb_seq in cb_seq_set:
                                                    tmp_atcg_set.setdefault(pileupread.alignment.query_sequence[pileupread.query_position],Counter())[cb_seq] += 1
                                                    #tmp_set[cb_seq] += 1
                                                    if pileupread.alignment.query_sequence[pileupread.query_position] == fa_dict[ch][v_pos]:
                                                        REF_cnt_dict.setdefault((ch, v_pos),[]).append(cb_seq)
                                                    if pileupread.alignment.query_sequence[pileupread.query_position] == ALT:
                                                        ALT_cnt_dict.setdefault((ch, v_pos),[]).append(cb_seq)
                                        else:
                                            if not pileupread.is_refskip:
                                                c_del += 1
                                indel_freq = c_del/float(c_keep+c_del)
                                #tmp_set = set()
                                #for b in tmp_atcg_set:
                                #    tmp_atcg_set[b] = set(it for it in tmp_atcg_set[b] if tmp_atcg_set[b][it]<=2)
                                #if (base_freq[0][0] in tmp_atcg_set) and (base_freq[1][0] in tmp_atcg_set):
                                #    tmp_set.update(tmp_atcg_set[base_freq[0][0]])
                                #    tmp_set.update(tmp_atcg_set[base_freq[1][0]])
                                #    rv = hypergeom(len(tmp_set), len(tmp_atcg_set[base_freq[0][0]]), len(tmp_atcg_set[base_freq[1][0]]))
                                #    hpg_prob = rv.pmf(len(tmp_atcg_set[base_freq[0][0]].intersection(tmp_atcg_set[base_freq[1][0]])))
                                #else:
                                hpg_prob = 1
                                reporting_summary.append((ch, v_pos, fa_dict[ch][v_pos], ALT, freq, s_freq, hpg_prob, seq_ent, indel_freq))
    print "number:", len(reporting_summary)
    subfolder_name = "{}_mutation".format(smp_id)
    if not os.path.exists(os.path.join(out_dir,subfolder_name)):
        os.makedirs(os.path.join(out_dir,subfolder_name))
    with gzip.open(os.path.join(out_dir,subfolder_name,"ref_cnt.csv.gz"),"wb") as ref_cnt_f:
        ref_cnt_f.write("chr,position,"+",".join(cb_seq_dict.keys())+"\n")  # write header
        for p in REF_cnt_dict:
            tmp_c = Counter(REF_cnt_dict[p])
            ref_cnt_f.write("{},{},".format(p[0],p[1])+",".join( str(tmp_c[it]) for it in cb_seq_dict.keys() )+"\n" )
    with gzip.open(os.path.join(out_dir,subfolder_name,"alt_cnt.csv.gz"),"wb") as alt_cnt_f:
        alt_cnt_f.write("chr,position,"+",".join(cb_seq_dict.keys())+"\n")  # write header
        for p in ALT_cnt_dict:
            tmp_c = Counter(ALT_cnt_dict[p])
            alt_cnt_f.write("{},{},".format(p[0],p[1])+",".join( str(tmp_c[it]) for it in cb_seq_dict.keys() )+"\n" )
    with gzip.open(os.path.join(out_dir,subfolder_name,"allele_stat.csv.gz"),"wb") as al_stat:
        al_stat.write("chr,position,REF,ALT,REF_frequency,REF_frequency_in_short_reads,hypergeom_test_p_value,sequence_entrophy,INDEL_frequency\n")  # write header
        for rec in reporting_summary:
            al_stat.write(",".join( str(it) for it in rec )+"\n" )
    pct_bin, pt = np.histogram(acc_pct, bins=500, range=(0, 1))
    with open(os.path.join(out_dir,subfolder_name,"freq_summary.csv"),"w") as cov_bin_out:
        for ix in range(500):
            cov_bin_out.write("{},{}\n".format(pt[ix],pct_bin[ix]))


def get_err_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, out_dir, min_cov=100):
    c2i = {"A":0, "C":1, "G":2, "T":3}  # four array.arrays of the same length in order A C G T
    fa_dict={}
    acc_pct = []
    REF_cnt_dict = {}
    ALT_cnt_dict = {}
    reporting_summary = []
    for c in get_fa(fa_f):
        fa_dict[c[0]] = c[1]
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    subfolder_name = "mutation"
    if not os.path.exists(os.path.join(out_dir,subfolder_name)):
        os.makedirs(os.path.join(out_dir,subfolder_name))
    mut_freq_f = gzip.open(os.path.join(out_dir,subfolder_name,"base_frequency.csv.gz"),"wb")
    for ch in chr_to_blocks:
        print ch
        homo_dict = find_homo_regions(fa_dict[ch], chr_to_blocks[ch])
        for ith, bl in enumerate(chr_to_blocks[ch]):
            tmp_bl_flat = get_gene_flat({"NNN":bl.transcript_list}, transcript_to_exon)
            for ex in tmp_bl_flat["NNN"]:
                cnt = bamfile.count(ch, ex[0], ex[1])
                if cnt < min_cov:
                    continue
                cov = bamfile.count_coverage(ch, ex[0], ex[1],
                quality_threshold=0)  # four array.arrays of the same length in order A C G T
                if len(cov[0])<20:
                    continue  # ignore tiny exons
                for i in range(5, len(cov[0])-5):  # ignore the bases at the beginning and the end (close to splicing site)
                    tot =  float(cov[0][i]+cov[1][i]+cov[2][i]+cov[3][i])
                    v_pos = ex[0]+i
                    if tot>min_cov and (fa_dict[ch][v_pos]!="N"):
                        freq = cov[c2i[fa_dict[ch][v_pos]]][i]/tot
                        seq_ent = seq_entropy(fa_dict[ch][(v_pos-10):(v_pos+10)])
                        for pileupcolumn in bamfile.pileup(ch, v_pos, v_pos+1,truncate=True, min_base_quality=0,ignore_overlaps=False,max_depth=1000):
                            indel_l = [1 if pileupread.is_refskip else 0 for pileupread in pileupcolumn.pileups]
                        indel_freq = indel_l.count(1)/float(len(indel_l))
                        in_homo= 1 if v_pos in homo_dict else 0
                        mut_freq_f.write("{},{},{},{},{},{}\n".format(ch,v_pos,freq,indel_freq,in_homo,seq_ent))
    mut_freq_f.close()

if __name__ == "__main__":
    smp_id = sys.argv[1]
    iso_dir = "/vast/scratch/users/peng.h/rachseq/fastq_polyA/filtered_bam/"
    if(os.path.isdir(iso_dir)):
        print(iso_dir)
    else:
        print(iso_dir,"not exist")
        exit(1)
    known_range_dict = {}
    #make SF3B1 position csv file
    for l in open("/vast/scratch/users/peng.h/rachseq/fastq_polyA/flames_out/SF3B1_position.csv").readlines()[1:]:
        ch ="chr" + l.strip().split(",")[-1]
        st = int(l.strip().split(",")[1])
        en = int(l.strip().split(",")[2])
        known_range_dict.setdefault(ch,[]).append((st,en))
    print known_range_dict
    known_position_dict = {("chr2",197402109):0}
    fa_f="/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/GRCh38.primary_assembly.genome.fa"
    gff_f="/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/gencode.v33.annotation.gff3"
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    bam_in=os.path.join(iso_dir, "{}_align2genome.bam".format(smp_id))
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open(os.path.join("/stornext/Genomics/data/CLL_venetoclax/RaCHseq_paper/barcode/", "{}.tsv".format(smp_id))))
    #os.chdir("/stornext/Genomics/data/CLL_venetoclax/RaCHseq_paper/filtered_polyA/flames_out/")
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict,
                        bam_short=None,known_position_dict=known_position_dict,
                        min_cov=90, report_pct=(0.1,0.8),known_range_dict=known_range_dict)

    """
    #RaCH-seq CLL141
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/Thijssen_count80/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short="/stornext/General/data/user_managed/grpu_mritchie_1/hongkePeng/Rachel/all_fastq/CLL141-CLL-cells_S8_Rsubread.sorted.bam"
    iso_dir = "/stornext/Genomics/data/CLL_venetoclax/RaCHseq/cll141/isoform_outs"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    """
    """
    # CLL141 capture
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/Thijssen_count80/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short="/stornext/General/data/user_managed/grpu_mritchie_1/hongkePeng/Rachel/all_fastq/CLL141-CLL-cells_S8_Rsubread.sorted.bam"
    iso_dir = "/stornext/Genomics/data/CLL_venetoclax/single_cell_data/capture_test/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)

    ### CLL141
    cb_seq_dict = dict( (it.strip().split(",")[1], it.strip().split(",")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Rachel_scRNA_Aug19/cluster_barcode_anno_lib20.csv"))
    bam_short="/stornext/General/data/user_managed/grpu_mritchie_1/hongkePeng/Rachel/all_fastq/CLL141-CLL-cells_S8_Rsubread.sorted.bam"
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Rachel_scRNA_Aug19/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)

    ### CLL267
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/patient2/Thijssen2_count20/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short="/stornext/General/data/user_managed/grpu_mritchie_1/hongkePeng/Rachel/all_fastq/CLL267_S4_Rsubread.sorted.bam"
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL267/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)

    ### CLL318
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/CLL318/CLL318_count20/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL318"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    ### CLL152

    cb_seq_dict = dict( (it.strip().split(",")[1], it.strip().split(",")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/CLL152/CLL152_count20/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL152/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)

    ### CLL153
    cb_seq_dict = dict( (it.strip().split(",")[0], it.strip().split(",")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/CLL153/cellranger_code/CLL153_count20/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL153/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)

    ### scmix1
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/scbench_5cellline_10x/10percent_cellranger/filtered_gene_bc_matrices/hg38/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_out_8"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    #get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    get_err_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir)

    ### scmix2
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixologyV3/luyiT_10X_260319/cellmix_Lib10/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION_April19/v3_long_read/isoform_all"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    """

    """
    #mouse
    known_position_dict = {}
    fa_f="/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/GRCm38.primary_assembly.genome.fa"
    gff_f="/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/gencode.vM24.annotation.gff3"
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)

    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/JamesRyall/10X/AGRF_CAGRF18671_CD27KANXX_cellranger/MuSC_10P_cellranger/filtered_gene_bc_matrices/mm10/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/JamesRyall/PromethION/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    """
