import json
import statistics as stats
from collections import defaultdict, Counter

from il_logging import Ilogger
from il_chk_conv import informed_round


ilog = Ilogger()
ilog.module = __name__

# def initialize_classes_logger(args):
#     ilog.verbosity = args.verbosity
class BlastCandidateRead:
    """Class for creation of BLAST candidate read objects"""
    def __init__(
        self, ID, RNAME, slen, length, constr_cov_range, read_cov_range, pident, mismatch, 
        gapopen, evalue, bitscore, qcovhsp, qcovus, constr_len, seq, orig_alignment, strand,
        args, aligner="BLAST"
        ):
        self.aligner = aligner
        self.args = args
        self.ID = ID
        self.RNAME = RNAME
        self.constr_cov_range = constr_cov_range
        self.slen = slen
        self.al_length = length
        self.read_cov_range = read_cov_range
        self.pident = pident
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.evalue = evalue
        self.bitscore = bitscore
        self.qcovhsp = qcovhsp
        self.qcovus = qcovus
        self.constr_len = constr_len
        self.trunc_constr = False
        # self.seq = seq
        # self.seq_wo_constr = self.gen_seq_wo_constr()
        # self.len_seq_wo_constr = len(self.seq_wo_constr)
        self.orig_alignment = orig_alignment
        self.strand = strand
        self.read_cov_len = abs(self.read_cov_range[1] - self.read_cov_range[0])
        self.multiple_sep_int = False
        self.multimer_int = False
        self.check_for_concatemers()
        self.trunc_constr_read_pos()
        self.seq = seq
        self.seq_wo_constr = self.gen_seq_wo_constr()
        self.len_seq_wo_constr = len(self.seq_wo_constr)
        self.constr_len_in_read = self.get_constr_len()
        self.len_wo_constr = self.slen - self.constr_len_in_read
        if abs(self.len_wo_constr-self.len_seq_wo_constr) > 5:
            ilog.vlprint(
                f"WARNING: Discrepancy between calculated readlength"
                f" without integration and length of constructed "
                f"sequence without int read:{self.RNAME}, len(seq): "
                f"{len(self.seq)}, len(seq_wo_constr): {self.len_seq_wo_constr}"
                f", calc_len-wo-constr: {self.len_wo_constr}. Possible reasons "
                f"include duplicate read-names/IDs in the read file/s or "
                f"fasta/q formatting errors.", 5
                )
        self.rtr = False

    def gen_seq_wo_constr(self):
        multi_read = any([self.multiple_sep_int, self.multimer_int])
        if not multi_read:
            if self.aligner == "BLAST":
                end_left = self.read_cov_range[0]-1
                start_right = self.read_cov_range[1]
            if self.aligner == "minimap2":
                end_left = self.read_cov_range[0]
                start_right = self.read_cov_range[1]+1
            return self.seq[:end_left] + self.seq[start_right:]
        else:
            if self.aligner == "BLAST":
                adjm = 0
                seq_wo = self.seq
                for i in self.orig_alignment:
                    seq_wo = seq_wo[:i[0]-adjm] + seq_wo[i[1]+1-adjm:]
                    adjm += i[1]-i[0]+1
            if self.aligner == "minimap2":
                # this will not accurately excise multiple separate integrations
                end_left = self.read_cov_range[0]
                start_right = self.read_cov_range[1]+1
            return seq_wo
    
    def get_constr_len(self):
        disquali = [self.trunc_constr, self.multiple_sep_int, self.multimer_int]
        if not any(disquali):
            return abs(self.read_cov_range[1] - self.read_cov_range[0])
        else:
            ranges = [r[1]-r[0] for r in self.orig_alignment]
            return sum(ranges) # /len(ranges)

    def trunc_constr_read_pos(self):
        # using <50 for truncated constructs at the end of the read assumes
        # that there are no adapters or barcodes left at the end of the reads;
        # if there were adpt or barcodes left, this value would have to be 
        # increased accordingly
        #
        # check if length of alignment to constr is shorter than constr
        # Beware: this approach leaves open a possibility for constructs
        # with repetetive sequenes to generate alignments that are counted
        # longer than the the read and thus be judged as full integrations,
        # although really being truncated.

        args = self.args
        if (self.al_length - self.constr_len) < 0:
                # check if the short construct is very close to 5'-end of read
                if self.trunc_constr == 35:
                    # if value previously set by check_for_concatemers
                    return "Contains truncated integration sequence at both ends"
                elif self.read_cov_range[0] < 50 + args.untrimmed:
                    self.trunc_constr = 5
                    return "Truncated integration sequence on 5'-end of read"
                # check if the short construct is very close to 3'-end of read
                elif self.read_cov_range[1] - self.slen < 50 + args.untrimmed:
                    self.trunc_constr = 3
                    return "Truncated integration sequence on 3'-end of read"
                # check if the short construct longer than 90% of construct-len
                elif (self.constr_len - self.al_length) > (self.constr_len/10):
                    self.trunc_constr = 1
                    return "Truncated integration shorter than 90% of search " \
                        "sequence in the middle of the read. " \
                        "Potential partial integration or false-positive"
                else:
                    return "Integration containing >90% of search sequence"
        else:
            # check if readlen is enough to harbor full int-construct, because
            # alignments can be longer than read e.g. for repetetive constructs
            if (self.slen - self.constr_len) > 0:
                if self.trunc_constr == 35:
                    return "Contains truncated integration sequence at both ends"
                else:
                    return "Complete Integration"
            else:
                if self.trunc_constr == 35:
                    return "Contains truncated integration sequence at both ends"
                else:
                    self.trunc_constr = 2
                    return "Truncated Integration. Long alignment, but read" \
                        " shorter than construct. Can happen e.g. for" \
                        "repetetive contructs."
                        

    def check_for_concatemers(self):
        if len(self.orig_alignment) > 1:
            pos_sep_ints = []
            psi_count = 1
            # check if significantly more than one construct-length is aligned to read
            if self.constr_len*1.2 < self.read_cov_len:
                for ir in range(len(self.orig_alignment)-1):
                    # check if individual alignments are separate ints 
                    if (self.orig_alignment[ir+1][0] - self.orig_alignment[ir][1]) > self.args.min_dist_ints:
                        psi_count += 1
                        pos_sep_ints.append(psi_count)
            elif len(self.orig_alignment)==2:
                # check if read could carry truncated construct on both ends
                if self.orig_alignment[0][0]<10 and (self.orig_alignment[1][1]-self.slen)<10:
                    self.trunc_constr = 35
                    psi_count += 1
                    pos_sep_ints.append(psi_count)
                    ilog.vlprint(
                        f"{self.RNAME}: Read contains truncated integration"
                        f"sequences at both ends", 2
                        )
            else:
                # only minimap2 alignment performed like for --intra and --poly
                pass
            if any(pos_sep_ints):
                self.multiple_sep_int = max(pos_sep_ints)
                ilog.vlprint(
                    f"{self.RNAME}: Read contains up to {max(pos_sep_ints)} "
                    f"separate insertions of the integrating DNA element", 2
                    )
                return 0
            else:
                self.multimer_int = True
                ilog.vlprint(
                    f"{self.RNAME}: Read contains an insertion of a "
                    f"multimer/concatemer of the integrating DNA element", 2
                    )
                return 0
        else:
            return 1
            
    def report_al_all(self):
        for attr in self.__dict__:
            print(attr, self.__dict__[attr])

    def report_al_sum(self):
        return [
                self.ID,
                self.RNAME,
                self.slen,
                self.read_cov_range,
                self.al_length,
                self.constr_cov_range,
                self.pident,
                self.mismatch,
                self.gapopen,
                self.evalue,
                self.bitscore,
                self.qcovhsp,
                self.qcovus,
                self.constr_len,
                ]

    def print_al_sum(self, sum_rep_file):
        al_sum = self.report_al_sum()
        sum_rep_file.write(" ".join([str(i) for i in al_sum]) + "\n")
        ilog.vlprint(" ".join([str(i) for i in al_sum]), 2, logging=False)
    
    def calc_qcov(self):
        qcov = self.al_length/self.qlen


class MmapCandidateRead(BlastCandidateRead):
    """Class for creation of minimap2 candidate read objects"""
    def __init__(
        self, ID, RNAME, slen, length, constr_cov_range, read_cov_range, pident,  
        mismatch, gapopen, evalue, bitscore, qcovhsp, qcovus, constr_len, seq,
        orig_alignment, strand, matches, mapq, al_type, minimizer, chain_score_prim,
        seq_div, len_rep_seeds, args
        ):
        super().__init__(
            ID, RNAME, slen, length, constr_cov_range, read_cov_range,   
            pident, mismatch, gapopen, evalue, bitscore, qcovhsp, qcovus, 
            constr_len, seq, orig_alignment, strand, args, aligner="minimap2"
        )
        
        # self.strand = strand
        self.mapq = mapq
        self.al_type = al_type
        self.minimizer = minimizer
        self.chain_score_prim = chain_score_prim
        self.seq_div = seq_div
        self.len_rep_seeds = len_rep_seeds
    

    def get_constr_len(self):
        # disquali = [self.trunc_constr, self.multiple_sep_int, self.multimer_int]
        # print(f"trunc_constr: {self.trunc_constr}")
        # if not any(disquali):
            # print("not any(disquali)")
        return abs(self.read_cov_range[1] - self.read_cov_range[0])
        # else:
        #     print("trunc or multi")
        #     ranges = [r[1]-r[0] for r in self.orig_alignment]
        #     return sum(ranges) # /len(ranges)


class LocalAlnm:
    def __init__(
        self, qaccver, qlen, saccver, slen, length, qstart, qend, sstart, send,
        pident, mismatch, gapopen, evalue, bitscore, qcovhsp, qcovus, score,
        gaps, nident, positive, sstrand
        ):
        self.read = qaccver
        self.read_len = int(qlen)
        self.chrom = saccver
        self.chrom_len = int(slen)
        self.al_len = int(length)
        self.read_start = int(qstart)
        self.read_end = int(qend)
        self.read_start_end = tuple(sorted([self.read_start, self.read_end]))
        self.chrom_start = int(sstart)
        self.chrom_end = int(send)
        self.chrom_start_end = tuple(sorted([self.chrom_start, self.chrom_end]))
        self.al_coord = (self.read_start_end, self.chrom_start_end)
        self.pident = float(pident)
        self.mismatch = int(mismatch)
        self.gapopen = int(gapopen)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)
        self.qcovhsp = int(qcovhsp)
        self.qcovus = int(qcovus)
        self.score = float(score)
        self.gaps = int(gaps)
        self.num_ident_match = int(nident)
        self.pos_scor_matches = int(positive)
        self.sstrand = sstrand
    
    def report_all(self):
        print(self.__dict__)
    
    def to_json(self):
        return self.__dict__


class ContAlnm:
    """LocalAlnm Version; Class for creation of contiguous alignment objects"""
    def __init__(self, readlen_wo_constr, chromosome, loc_alnm_objs_lst): 
        self.loc_alnm_objs = list(enumerate(loc_alnm_objs_lst))
        self.readlen_wo_constr = readlen_wo_constr
        self.chromosome = chromosome
        self.al_lens = {k:v.al_len for k,v in self.loc_alnm_objs}
        self.al_coord = {k:v.al_coord for k,v in self.loc_alnm_objs}
        self.idents = {k:v.pident for k,v in self.loc_alnm_objs}
        self.mismatches = {k:v.mismatch for k,v in self.loc_alnm_objs}
        self.gaps = {k:v.gapopen for k,v in self.loc_alnm_objs}
        self.evals = {k:v.evalue for k,v in self.loc_alnm_objs}
        self.bitscores = {k:v.bitscore for k,v in self.loc_alnm_objs}
        self.qcovhsps = {k:v.qcovhsp for k,v in self.loc_alnm_objs}
        self.qcovuss = max(set([v.qcovus for k,v in self.loc_alnm_objs]))
        self.score = {k:v.score for k,v in self.loc_alnm_objs}
        self.gaps = {k:v.gaps for k,v in self.loc_alnm_objs}
        self.matches = {k:v.num_ident_match for k,v in self.loc_alnm_objs}
        self.sstrand = {k:v.sstrand for k,v in self.loc_alnm_objs}

    def _check_al_contiguity(self):
        """method to select and keep only the most contiguous group of
        individual local alignments, with respect to chromosome coordinates"""
        coords = sorted(self.al_coord.values(), key=lambda x: x[1])
        prev = 0
        curr = 0
        nxt = 0
        discontinuities = 0
        keep_separate = dict()
        for coord_pair in coords:
            nxt = coord_pair[1]
            if prev == 0:
                keep_separate[nxt] = [nxt]
                curr = nxt
            else:
                if nxt[0] - prev[1] > self.readlen_wo_constr / 4:
                    discontinuities += 1
                    keep_separate[nxt] = [nxt]
                    curr = nxt
                else:
                    keep_separate[curr].append(nxt)
            prev = nxt
        if discontinuities == 0:
            pass
        else:
            keepers = max(
                keep_separate,
                key=lambda x: sum(
                        [(abs(start-end)) for start,end in keep_separate[x]]
                        )
                    )
            keepers = keep_separate[keepers]
            self.al_coord = {
                k:v for k,v in self.al_coord.items() if v[1] in keepers
                }

    def _filter_redundant_als(self):
        """Filter out redundant alignmets of repetetive sequences, which lie 
        within the range of the size of the int-region and therefore are not 
        excluded through distance alone"""

        def check_overlap(range_x: tuple, range_y: tuple):
            lower = min(range_x[1], range_y[1])
            upper = max(range_x[0], range_y[0])
            overlap = lower - upper + 1
            if overlap <= 0:
                return False
            else:
                return overlap
        
        def chk_ovl_size(region_range, overlap):
            reg_len = region_range[1] - region_range[0] + 1
            ovl_proportion = overlap/reg_len * 100
            if ovl_proportion > 10:
                return True
            else:
                return False

        overlapping = defaultdict(list)
        for coopair in self.al_coord.values():
            others = [
                cpair for cpair in self.al_coord.values() if cpair != coopair
                ]
            for other in others:
                # check overlap between aligned range on query
                overlap = check_overlap(coopair[0], other[0])
                if overlap:
                    relevant = chk_ovl_size(coopair[0], overlap)
                    if relevant:
                        overlapping[coopair].append(other)
        
        compare_ovl_goups = [
            tuple(sorted([key, *overlapping[key]])) for key in overlapping
            ]
        compare_ovl_goups = set(compare_ovl_goups)
        for group in compare_ovl_goups:
            gap_sums = {}
            for redundant_al in group:
                tmp_exclude = [cop for cop in group if cop != redundant_al]
                test_coord_set = sorted(
                    [cop for cop in self.al_coord.values() if cop not in tmp_exclude],
                    key=lambda x: x[1][0]
                    )
                gaps_len = 0
                if len(test_coord_set) < 2:
                    pass
                else:
                    for idx in range(len(test_coord_set)-1):
                        gaps_len += test_coord_set[idx+1][1][0] - test_coord_set[idx][1][1]
                    gap_sums[redundant_al] = gaps_len
            if len(gap_sums) == 0:
                try:
                    best_al_obj = max(
                        [(k,o) for k,o in self.loc_alnm_objs if o.al_coord in group],
                        key=lambda x: x[1].bitscore
                        )
                    ba_key, best_al_coord = best_al_obj[0], best_al_obj[1].al_coord
                    exclude_keys = [k for k, v in self.loc_alnm_objs if k != ba_key]
                    self._discard_excluded(exclude_keys)
                except ValueError as e:
                    ilog.vlprint(
                        f"Unable to pick best_al_obj for {self.loc_alnm_objs[0][1].read}"
                        f" by max bitscore. {e}", 2
                        )
            else:
                try:
                    min_gap_al = min(
                        [(k,v) for k,v in gap_sums.items()], key=lambda x: x[1]
                        )
                    exclude_coord = [cop for cop in group if cop != min_gap_al]
                    exclude_keys = [k for k, v in self.loc_alnm_objs if v.al_coord in exclude_coord]
                    self._discard_excluded(exclude_keys)
                except ValueError as e:
                    ilog.vlprint(
                        f"Unable to pick min_gap_al for {self.loc_alnm_objs[0][1].read}. {e}", 1
                        )
    
    def report_all(self):
        print(self.__dict__)
    
    def to_json(self):
        return self.__dict__
    
    def _discard_excluded(self, excluded_keys: list):
        for attr in self.__dict__:
            if type(self.__dict__[attr]) == dict:
                 self.__dict__[attr] = {
                     k:v for k,v in self.__dict__[attr].items() if k not in excluded_keys
                     }
        self.loc_alnm_objs = [
            local for local in self.loc_alnm_objs if local[0] not in excluded_keys
            ]

    def _discard_range_excluded(self):
        excluded_keys = []
        for local in self.loc_alnm_objs:
            if local[1].al_coord not in self.al_coord.values():
                excluded_keys.append(local[0])
        self._discard_excluded(excluded_keys)

    def _listify(self):
        for attr in self.__dict__:
            if type(self.__dict__[attr]) == dict:
                self.__dict__[attr] = list(self.__dict__[attr].values())
    
    def refine_self(self):
        self._check_al_contiguity()
        self._discard_range_excluded()
        # best to _filter_redundant_als after range-exclusion for efficiency 
        self._filter_redundant_als()
        self._listify()


class ReadToRefAl:
    """Class for creation of read to reference alignment objects"""
    def __init__(self, READ, read_len, al_dict, cand_read_obj, args):
        self.cand_read_obj = cand_read_obj
        self.READ = READ
        self.trunc_constr = self.cand_read_obj.trunc_constr
        self.read_len = read_len
        self.readlen_wo_constr = self.cand_read_obj.len_wo_constr
        self.constr_proportion = (1-(self.readlen_wo_constr/self.read_len))*100
        self.all_alignments = [
            ContAlnm(self.readlen_wo_constr, chrom_key, al_dict[chrom_key])
            for chrom_key in al_dict
            ]
        self.refine_alignments()
        self.best_al = self.best_alignment(self.all_alignments)
        self.sstrand = tuple(set(self.best_al.sstrand))
        self.CHR = self.best_al.chromosome[0]
        self.chr_len = int(self.best_al.chromosome[1])
        self.als_lengths = sum(self.best_al.al_lens)
        self.read_start_end_chr_start_end = sorted(self.best_al.al_coord, key=lambda x: x[1])
        self.read_start_end, self.chr_start_end = self.get_start_end()
        self.genomic_flanking = len(self.best_al.loc_alnm_objs)
        self.flank_side = False
        self.pident = sum(self.best_al.idents)/len(self.best_al.idents)
        self.mismatch = self.best_al.mismatches
        self.gapopen = self.best_al.gaps
        self.evalue = self.best_al.evals
        self.bitscore = self.best_al.bitscores
        self.qcovhsp = self.best_al.qcovhsps
        self.qcovus = self.best_al.qcovuss
        self.cor_qcovus = self.best_al.qcovuss + self.constr_proportion
        self.gaps = self.best_al.gaps
        self.score = self.best_al.score
        self.matches = self.best_al.matches
        self.int_coord, self.ins_gap, self.coord_evid, self.split_coord = self.chrom_coord(args)
        # self.lead_coord = False
        self.assoc_nucl_sites = []
        # multi-int-read info from il_evaluate: mult_ind_int_read_info()
        self.multi_info = {"segments":[], "ins_chrom_pos":{}, "delta_sites": {}}
        # segments ordered list of all constr-to-read and read-to-ref alignments
        # ins_chrom_pos: segments[idx] of ins as keys and rtr-borders as vals
        # delta_sites: {(segments[idx], site-pos): abs(rtr-borders - site)}

    def refine_alignments(self):
        """Run refine_self method of all ContAlnm objects to remove the
        non/least contiguos and redundant individual local alignments"""
        for al in self.all_alignments:
            al.refine_self()   

    def best_alignment(self, al_lst):
        """Select best contiguous alignment location for each read"""
        # This method is made obsolete for integration_locator module by 
        # refine_alignments as long as culling_limit == 1.
        # However, it is still put to good use in il_evaluate, so keep around!

        if len(al_lst) > 0:
            # changed selection criterion for best_alignment for the sake
            # of il_evaluate from sum of al-length to sum of bitscores
            # return max(al_lst, key=lambda al_obj: sum(al_obj.al_lens))
            return max(al_lst, key=lambda al_obj: sum(al_obj.bitscores))
        else:
            return False

    def get_start_end(self):
        read_coord = []
        chrom_coord = []
        for tup_pair in self.read_start_end_chr_start_end:
            read_coord.extend(tup_pair[0])
            chrom_coord.extend(tup_pair[1])
        read_start_end = (min(read_coord), max(read_coord))
        chrom_start_end = (min(chrom_coord), max(chrom_coord))
        return read_start_end, chrom_start_end
    
    # def get_flank_sides(self):
    #     if self.genomic_flanking>1:
    #         return "both"
    #     elif self.genomic_flanking==1:
    #         start, end = self.chr_start_end
    #         if self.trunc_constr == 5:

    #     else:
    #         return False

    def chrom_coord(self, args):
        read_north_dict = {}
        read_south_dict = {}
        for n in range(len(self.read_start_end_chr_start_end)):
            tmp_south = []
            south_end1 = self.cand_read_obj.read_cov_range[0] - self.read_start_end_chr_start_end[n][0][0]
            south_end2 = self.cand_read_obj.read_cov_range[0] - self.read_start_end_chr_start_end[n][0][1]
            tmp_south.append(abs(south_end1))
            tmp_south.append(abs(south_end2))
            south_end = min(tmp_south)
            read_south_dict[n] = south_end
            tmp_north = []
            north_end1 = self.cand_read_obj.read_cov_range[1] - self.read_start_end_chr_start_end[n][0][0]
            north_end2 = self.cand_read_obj.read_cov_range[1] - self.read_start_end_chr_start_end[n][0][1]
            tmp_north.append(abs(north_end1))
            tmp_north.append(abs(north_end2))
            north_end = min(tmp_north)
            read_north_dict[n] = north_end

        chr_close_readsouth = self.read_start_end_chr_start_end[min(read_south_dict,
                                                                    key=lambda x: read_south_dict[x])][1]
        chr_close_readnorth = self.read_start_end_chr_start_end[min(read_north_dict,
                                                                    key=lambda x: read_north_dict[x])][1]
        
        coord_lst = [chr_close_readsouth, chr_close_readnorth]

        # These are the values between which the integration sits in chrom
        if self.sstrand[0]=="minus" and len(set(coord_lst))>1:
            int_coord = [coord_lst[0][0], coord_lst[1][1]]
        else:
            int_coord = [coord_lst[0][1], coord_lst[1][0]]

        ins_gap = abs(int_coord[0] - int_coord[1])

        if self.trunc_constr:
            ins_gap = None
            if self.trunc_constr == 5:
                if self.sstrand[0] == "minus":
                    self.flank_side = "left"
                    int_coord[1] = int_coord[0]
                else:
                    self.flank_side = "right"
                    int_coord[0] = int_coord[1]   
            elif self.trunc_constr == 3:
                if self.sstrand[0] == "minus":
                    self.flank_side = "right"
                    int_coord[0] = int_coord[1]
                else:
                    self.flank_side = "left"
                    int_coord[1] = int_coord[0]
            else:
                pass
        else:
            # if there is a truncated or full-length construct right at one or 
            # both read end/s or if for any other reason the int-seq is not 
            # flanked by genomic sequences on both side
            if self.genomic_flanking < 2:
                ins_gap = None
                read_cov = self.cand_read_obj.read_cov_range
                five_pri_dst = abs(0-read_cov[0])
                three_pr_dst = abs(read_cov[1]-self.cand_read_obj.slen)
                if five_pri_dst < three_pr_dst:
                    if self.sstrand[0] == "minus":
                        self.flank_side = "left"
                    else:
                        self.flank_side = "right"
                else:
                    if self.sstrand[0] == "minus":
                        self.flank_side = "right"
                    else:
                        self.flank_side = "left"
            else:
                pass

               

            
            # if there is a truncated or full-length construct right at one or 
            # both read end/s
            # cro = self.cand_read_obj
            # if cro.read_cov_range[0] <= 25:
            #     ins_gap = None
            # elif cro.slen-cro.read_cov_range[1] <= 25:
            #     ins_gap = None
            # else:
            #     pass

        # evid_qual = quality of coordinate evidence, the lower the value the 
        # higher the quality; RTR-alignments will be sorted according to this
        # values for int-site nucleation in gen_integration_sites.
        # 1 = "double_crisp" (two-sided, edges less than 5bp apart)
        # 2 = "double_ok" (two-sided, edges less than min_dist_ints apart)
        # 3 = "double_unsharp" (two-sided, edges more than min_dist_ints apart)
        # 4 = "one_sided" (int-read with only one genome/int border because read
        # ends with int-seq [usually truncated, or more rarely full-len right at
        # the read terminus])
        coord_evid = None
        if ins_gap == 0:
            coord_evid = 1 # "double_crisp"
        elif ins_gap == None:
            coord_evid = 4 # "one_sided"
        else:
            if ins_gap <= 5:
                coord_evid = 1 # "double_crisp"
            elif ins_gap > args.min_dist_ints:
                coord_evid = 3 # "double_unsharp"
            else:
                coord_evid = 2 # "double_ok"

        return tuple(int_coord), ins_gap, coord_evid, coord_lst
    
    def to_json(self):
        return self.best_al


class NucleationSite():
    """Class for nucleation of read-to-reference alignment seeds"""
    def __init__(self, chrom, rtr_al_seed, args):
        self.args = args
        self.chrom = chrom
        self.prim_seed = rtr_al_seed
        self.seed_type = self.prim_seed.coord_evid
        self.seed_coord_a = self.prim_seed.int_coord[0]
        self.seed_coord_b = self.prim_seed.int_coord[1]
        self.id = self.gen_id()
        self.linked_als = []
        self.linked_types = []
        self.linked_coords = []
        self.seed_coord_dists = {"A": [], "B": []}
        self.avg_dist = None
        self.avg_coord_evid = None
        self.all_als = []
        self.all_types = []
        self.all_coords = []


    def gen_id(self):
        transl_dct = {
            1: "double_crisp",
            2: "double_ok",
            3: "double_unsharp",
            4: "one_sided"
            }
        id = str(f"{self.chrom}_{self.seed_coord_a};{self.seed_coord_b}_"
            f"seed_{transl_dct[self.seed_type]}")
        return id
    
    def compare_incoming_rtr(self, rtr_al):
        dists_a = {}
        dists_b = {}
        for coord in rtr_al.int_coord:
            dists_a[abs(self.seed_coord_a - coord)] = coord
            dists_b[abs(self.seed_coord_b - coord)] = coord
        min_dist_a = min(dists_a.keys())
        min_dist_b = min(dists_b.keys())
        min_dists = [min_dist_a, min_dist_b]
        matching = [True for dist in min_dists if dist<self.args.min_dist_ints]
        if any(matching):
            coord_link = dists_a[min_dist_a] if min_dist_a<min_dist_b else dists_b[min_dist_b]
            self.add_rtr_alnm(rtr_al, coord_link, min_dist_a, min_dist_b)
            return True
        else:
            return False

    def add_rtr_alnm(self, rtr_al, coord_link, min_dist_a, min_dist_b):
        self.linked_als.append(rtr_al)
        self.linked_coords.append(coord_link)
        self.seed_coord_dists["A"].append(min_dist_a)
        self.seed_coord_dists["B"].append(min_dist_b)
        self.linked_types.append(rtr_al.coord_evid)
        # also attach the Nucleation-Site to the rtr_al to keep track of how 
        # many nucleation sites a rtr-al is associated with
        rtr_al.assoc_nucl_sites.append(self.id)
    
    def update_alnm_incl_seed_vals(self):
        self.all_als = self.linked_als
        self.all_als.append(self.prim_seed)
        self.all_types = self.linked_types
        self.all_types.append(self.seed_type)
        # The problem with adding the seed coordinates to an common list with
        # the linked coords is that from the linked ones only one coord is added 
        # (the better one); if the seed coordinates would be unsharp, the total 
        # location quality might suffer if both borders were included (especially 
        # if only the better border had attracted linked coord before). Possibilities
        # for handling this issue: 1. only add the coordinate that has the smallest delta 
        # to the average linked coordinates value; 2. Add both values to all_coords 
        # and than in later steps check which has the smaller standard deviation 
        # (all_coords or linked coords) and use the one with the smaller deviation
        # self.all_coords = self.linked_types
        incl_a = list(self.linked_coords)
        incl_a.append(self.seed_coord_a)
        incl_b = list(self.linked_coords)
        incl_b.append(self.seed_coord_b)
        incl_ab = list(self.linked_coords)
        incl_ab.append(self.seed_coord_a)
        incl_ab.append(self.seed_coord_b)
        incl_lsts = {"incl_a": incl_a, "incl_b": incl_b, "incl_ab": incl_ab}
        # the problem with the stdevs is when the list are shorter then three
        incl_lsts = {k:v for k,v in incl_lsts.items() if len(v)>2}
        incl_sds = {k:stats.stdev(v) for k,v in incl_lsts.items()}

        if len(incl_sds) > 0:
            minsd_coord_incl = min(incl_sds.items(), key=lambda x: x[1])[0]
            self.all_coords = incl_lsts[minsd_coord_incl]
        else:
            self.all_coords = self.linked_coords
            self.all_coords.extend(incl_ab)


    def calc_metrics(self):
        if len(self.seed_coord_dists["A"]) > 0:
            avg_dists_a = sum(self.seed_coord_dists["A"])/len(self.seed_coord_dists["A"])
        else:
            avg_dists_a = 1000000000
        if len(self.seed_coord_dists["B"]) > 0:
            avg_dists_b = sum(self.seed_coord_dists["B"])/len(self.seed_coord_dists["B"])
        else:
            avg_dists_b = 1000000000
        self.avg_dist = min([avg_dists_a, avg_dists_b])
        self.update_alnm_incl_seed_vals()
        if len(self.linked_als) > 0:
            # self.update_alnm_incl_seed_vals()
            self.avg_coord_evid = sum(self.all_types)/len(self.all_types)
        else:
            self.avg_coord_evid = self.seed_type
        # all_types = self.linked_types
        # all_types.append(self.seed_type)
        # self.avg_coord_evid = sum(self.linked_types)/len(self.linked_types)
    
    def report_nucl_site(self):
        ilog.vlprint(self.id, 5, logging=False)
        ilog.vlprint(f"avg_dist {self.avg_dist}", 5, logging=False)
        ilog.vlprint(f"coord_evidence {self.avg_coord_evid}", 5, logging=False)
        ilog.vlprint(f"supporting alnms {len(self.all_als)}", 5, logging=False)
        ilog.vlprint("##################", 5, logging=False)


class PotIntSite():
    """Class for creation of potential integration site objects"""
    def __init__(self, chrom, chrom_coord, supp_als, blast_params, args):
        self.chrom = chrom
        self.chrom_coord = chrom_coord
        self.supp_als = supp_als
        self.supp_absorbed_fp = []
        self.sample = args.sample
        self.args = args
        self.prune_alignments()
        self.blast_params = blast_params
        self.read_names = self.get_blast_readnames()
        self.app_loc_sd = "n.a."
        self.blast_loc = False
        self.mm2_loc = False
        self.mm2_chrom = False
        self.mm_only_cand_reads = {}
        self.approx_loc = self.approx_prec_loc()
        self.del_at_junction = False
        self.name = self.chrom + "_" + str(self.approx_loc)
        self.old_name = False
        self.merge_with = set()
        self.absorbed = False
        self.av_supp_qual = None
        self.achieved_ident_matches = self.calc_achieved_ident_matches()
        self.max_aims = self.calc_max_achievable_ident_matches()
        self.ident_match_ratio = self.calc_imr()
        self.ssbr = "n.a."
        self.minimap_supp = {}
        self.minimap_supp_sr = {}
        self.int_als = {}
        self.pot_non_int_als = {}
        self.conf_non_int_als = []
        self.proxhom_non_int_als = {}
        self.sum_non_int_als = int()
        self.total_cov = int()
        self.spanned_region = {"blast": None, "minimap2": None}
        self.spanned_seq = ""
        # the integration is part of a site which contains concatemers/multimers
        # or multiple closely spaced separate integrations sites
        self.multi_site = False
        # an integration exists in a similar, but non-identical, location on
        # a homologous chromosome
        self.prox_homolog = []
        # mrli_supp = additional support inferred from multi-read linkage 
        self.mrli_supp = {}
        # multread_linked = other PotIntSites linked by multireads
        self.multread_linked = {}
        self.multi_sep_supp = set()
        self.multimer_supp = set()
        self.total_support = self.calc_total_support()
        self.seq_loss_info = False
        self.seq_gain_info = False

    def get_blast_readnames(self):
        return [al.READ for al in self.supp_als]
    
    def get_spanned_region_blast(self):
        """Determine genomic region spanned by mapped candidate reads"""
        corner_vals = [val for rtr in self.supp_als for val in rtr.split_coord]
        corner_vals = sorted([val for tup in corner_vals for val in tup])
        max_span = (corner_vals[0], corner_vals[-1])
        self.spanned_region["blast"] = max_span

    def calc_total_support(self):
        total_support = []
        total_support.extend(self.get_blast_readnames())
        total_support.extend(list(self.minimap_supp.keys()))
        total_support.extend(list(self.mrli_supp.keys()))
        return set(total_support)

    def calc_total_cov(self, mode="potential"):
        # int_reads = max([len(self.supp_als), len(self.non_int_als)])
        pot_non = self.pot_non_int_als
        conf_non = self.conf_non_int_als
        ph_non = self.proxhom_non_int_als
        if mode=="potential":
            self.total_cov = len(self.int_als) + len(pot_non)
            self.sum_non_int_als = len(pot_non)
        elif mode=="confirmed":
            self.total_cov = len(self.int_als) + len(conf_non)
            self.sum_non_int_als = len(conf_non)
        elif mode=="all":
            self.total_cov = len(self.int_als) + len(conf_non) + len(ph_non)
            self.sum_non_int_als = len(conf_non) + len(ph_non)

    def calc_achieved_ident_matches(self):
        return sum([sum(rtr.matches) for rtr in self.supp_als])

    def calc_max_achievable_ident_matches(self):
        "Determine cumulative length of supporting reads without integration"
        max_scores = []
        for rtr_al in self.supp_als:
            cro = rtr_al.cand_read_obj
            eff_len = cro.len_wo_constr
            max_scores.append(eff_len)
        return sum(max_scores)
    
    def calc_imr(self):
        "Calculate identical matches ratio"
        try:
            imr = self.achieved_ident_matches / self.max_aims
            return imr
        except ZeroDivisionError:
            return "n.a." # "no RTRs"

    def prune_alignments(self):
        if self.sample == "clonal":
            self.supp_als = [al for al in self.supp_als if al.qcovus >= 25]
        elif self.sample == "polyclonal":
            self.supp_als = [al for al in self.supp_als if al.qcovus >= 1]
        else:
            pass

    def eval_supp_qual(self):
        "Evaluate quality of supporting alignments"
        cor_qcovuss = []
        for aln in self.supp_als:
            cor_qcovuss.append(aln.cor_qcovus)
        if len(cor_qcovuss) > 0:
            av_qcov = sum(cor_qcovuss)/len(cor_qcovuss)
            return av_qcov
        else:
            return 0

    def dump_reads(self):
        dumped_reads = []
        with open(self.name + "_supporting_reads.fasta", "w") as out:
            for al in self.supp_als:
                dumped_reads.append(">" + al.cand_read_obj.RNAME + "\n")
                out.write(">" + al.cand_read_obj.RNAME + "\n")
                out.write(al.cand_read_obj.seq)
                if not al.cand_read_obj.seq.endswith("\n"):
                    out.write("\n")
            for al in self.minimap_supp:
                if al not in dumped_reads:
                    dumped_reads.append(">" + al + "\n")
                    out.write(">" + al + "\n")
                    out.write(self.mm_only_cand_reads[al].seq + "\n")
                else:
                    pass
        return dumped_reads
    
    def eliminate_outlier(self, spl):
        diff_dic = {}
        med = stats.median(spl)
        for idx, num in enumerate(spl):
            diff_dic[abs(num - med)] = idx
        outlier = max(diff_dic)
        spl.pop(diff_dic[outlier])
        return spl

    def condense_coord(self, al_egdes_lst):
        """Condense al_edge data to data points lying closest together keeping 
        at least 1/3 of the data"""
        # Alignment edges converge on the actual integration-site from both 
        # sides. The distances of the alignment-edge data to the actual 
        # integration-site cluster in two groups, one close around the median 
        # corresponding to alignments that continue up to the actual integration
        # -site and one formed by alignments that break off significantly  
        # farther away, due to sequencing and/or alignment errors. Although for 
        # high coverage data using the median of all edges already generates
        # accurate results in most cases, with lower coverages these outliers  
        # can distort the approximation of the integration-site location 
        # substantially. Based on this observation the data is condensed to the 
        # third of data clustering closest around the median, for a first
        # approximation of integration-site location.
        spl = sorted(al_egdes_lst)
        orig_len = len(spl)
        curr_med = stats.median(spl)
        max_dist = max([abs(val - curr_med) for val in spl])
        while len(spl) >= orig_len*0.3:
            if len(spl) > 3 and max_dist > 25:
                self.eliminate_outlier(spl)
            else:
                break 
        return spl
    
    def calc_blast_loc(self): # incl_multi=False
        both_borders = []
        double_trunc = []
        for alignment in self.supp_als:
            if alignment.trunc_constr==5 or alignment.trunc_constr==3:
                both_borders.append(alignment.int_coord[0])
            elif alignment.trunc_constr==35:
                double_trunc.append(alignment.int_coord)
            else:
                both_borders.append(alignment.int_coord[0])
                both_borders.append(alignment.int_coord[1])

        # if incl_multi:
        #     for rtr in self.mrli_supp:
        #         rtro, idx = self.mrli_supp[rtr]
        #         both_borders.append(rtro.multi_info["ins_chrom_pos"][idx])
        #         print(
        #             self.name, f'Appended {rtro.multi_info["ins_chrom_pos"][idx]}'
        #             )

        if len(double_trunc) > 0:
            if len(both_borders) > 0:
                avg = sum(both_borders)/len(both_borders)
                for tup in double_trunc:
                    both_borders.append(min([abs(v-avg) for v in tup]))
            else:
                for tup in double_trunc:
                    both_borders.append(tup[0])
                    both_borders.append(tup[1])
       
        if len(self.supp_als) > 2:
            al_close_loc = self.condense_coord(both_borders)
            loc_median = round(stats.median(al_close_loc))
            self.app_loc_sd = round(stats.stdev(al_close_loc))
            self.blast_loc = loc_median
            return loc_median
        elif len(self.supp_als) == 0:
            return False
        else:
            loc_mean = round(stats.mean(both_borders))
            self.blast_loc = loc_mean
            return loc_mean
    
    def calc_mm2_loc(self):
        calc_pos = defaultdict(list)
        for cro in self.mm_only_cand_reads.values():
            pos_in_read_wo = cro.read_cov_range[0]
            chrom_name = self.minimap_supp[cro.RNAME]["tname"]
            tstart = self.minimap_supp[cro.RNAME]["tstart"]
            qstart = self.minimap_supp[cro.RNAME]["qstart"]
            pos_in_chrom = tstart + (pos_in_read_wo - qstart)
            calc_pos[chrom_name].append(pos_in_chrom)
        chrom, positions = max(calc_pos.items(), key=lambda x: len(x[1]))
        self.mm2_chrom = chrom
        if len(positions) > 2:
            al_close_loc = self.condense_coord(positions)
            loc_median = round(stats.median(al_close_loc))
            # self.app_loc_sd = round(stats.stdev(al_close_loc))
            self.mm2_loc, final_pos = loc_median, loc_median
        elif len(positions) == 0:
            return False
        else:
            loc_mean = round(stats.mean(positions))
            self.mm2_loc, final_pos = loc_mean, loc_mean
        # return chrom, final_pos
        return final_pos   

    def approx_prec_loc(self):
        """Approximate precise integration locations"""
        
        ilog.vlprint(
            'Approximate precise integration locations', 6, logging=False
            )
        if self.sample == "clonal":
            blast_loc = self.calc_blast_loc()
            return blast_loc
        elif self.sample == "polyclonal" or self.args.quick_mm2:
            if len(self.supp_als) > 0:
                self.blast_loc = self.calc_blast_loc()
            if len(self.mm_only_cand_reads) > 0:
                self.mm2_loc = self.calc_mm2_loc()
            if self.blast_loc:
                return self.blast_loc
            else:
                return self.mm2_loc
        else:
            blast_loc = self.calc_blast_loc()
            return blast_loc
    
    def calc_del_at_junction(self):
        """Determine (genomic) sequence loss at integration site"""
        int_gaps = [
            (rtr.READ, rtr.ins_gap) for rtr in self.supp_als if rtr.ins_gap!=None
            ]
        gap_coord = [rtr.int_coord for rtr in self.supp_als if rtr.ins_gap!=None]
        de_left = [t[0] for t in gap_coord]
        de_right = [t[1] for t in gap_coord]
        if len(int_gaps) > 0:
            avg_del = sum([t[1] for t in int_gaps])/len(int_gaps)
            med_del = stats.median([t[1] for t in int_gaps])
            min_del =  min([t[1] for t in int_gaps]) - 1
            min_del = min_del if min_del>-1 else 0
        else:
            avg_del = "n.a."
            min_del = "n.a."

        one_sided = [rtr for rtr in self.supp_als if rtr.genomic_flanking==1]
        lefties = [r.int_coord[0] for r in one_sided if r.flank_side=="left"]
        righties = [r.int_coord[0] for r in one_sided if r.flank_side=="right"]
        if all([len(lefties)>0, len(righties)>0]):
            singlegap = stats.median(righties)-stats.median(lefties)-1
            singlegap = singlegap if singlegap>-1 else 0
        else:
            singlegap = "n.a."
        
        if all([min_del=="n.a.", singlegap=="n.a."]):
            self.del_at_junction = "n.a."
        elif any([min_del=="n.a.", singlegap=="n.a."]):
            self.del_at_junction = [v for v in (min_del, singlegap) if v!="n.a."][0]
        else:
            self.del_at_junction = min((min_del, singlegap))
        
        all_left = list(lefties)
        all_left.extend(de_left)
        # sort_allft = sorted(all_left)
        all_right = list(righties)
        all_right.extend(de_right)
        # sort_alrgt = sorted(all_right)
        zip_all = list(zip(sorted(all_left), sorted(all_right)))
        delta_gaps = [abs(x-y)-1 for x,y in zip_all]
        delta_gaps = [v if v>=0 else 0 for v in delta_gaps]
        most_freq_gap = max(delta_gaps, key=delta_gaps.count)
        # print("delta_gaps", delta_gaps)
        # sd_gaps = stats.stdev([abs(x-y) for x,y in zip_all])
        # print("sd_gaps", sd_gaps)
        if all([len(all_left)>0, len(all_right)>0]):
            allgap = stats.median(all_right)-stats.median(all_left)-1
            allgap = allgap if allgap>-0.1 else 0
            if len(delta_gaps) >= 2:
                sd_gaps = stats.stdev(delta_gaps)
            else:
                sd_gaps = "n.a."
        else:
            allgap = "n.a."
        
        self.del_at_junction = 0 if self.del_at_junction<0 else self.del_at_junction
        self.seq_loss_info = {
            "double_ended": {
                "count": len(int_gaps),
                "avg_del": avg_del,
                "med_del": med_del,
                "min_del": min_del
                },
            "single_ended": {
                "count_all": len(one_sided),
                "count_left": len(lefties),
                "count_right": len(righties),
                "singlegap": singlegap, 
                },
            "min_loss": self.del_at_junction,
            "med_loss": allgap,
            "sd_loss": sd_gaps,
            "most_freq_gap": most_freq_gap
        }

    def calc_seq_gain(self):
        """
        Determine potential sequence gains (in addition to
         provided integrating-sequence) at integration site 
        """
        gains_left = []
        gains_right = []
        gain_info = {
            "left": {
                "min": "n.a",
                "median": "n.a",
                "most_freq": "n.a"
                },
            "right": {
                "min": "n.a",
                "median": "n.a",
                "most_freq": "n.a"
                }
        }

        def get_junction_delta(gen_cov_left, gen_cov_right, int_cov_rng, inv=False):
            if gen_cov_left:
                gain_left = int_cov_rng[0] - gen_cov_left[1] - 1
                if not inv:
                    gains_left.append(gain_left)
                else:
                    gains_right.append(gain_left)
            if gen_cov_right:
                gain_right = gen_cov_right[0] - int_cov_rng[1] - 1
                if not inv:
                    gains_right.append(gain_right)
                else:
                    gains_left.append(gain_right)

        read_ct = 0
        for rtr in self.supp_als:
            cro = rtr.cand_read_obj
            int_cov_rng = cro.read_cov_range
            gen_cov = sorted(
                    rtr.read_start_end_chr_start_end, key=lambda x: x[0][0]
                    )
            segments = rtr.multi_info["segments"]
            if len(rtr.multi_info["segments"]) > 0:
                read_ct += 1
                idx = -1
                fail_msg = f"Seq gain calc for {rtr.READ} in {self.name} failed"\
                            f" because int-range not flanked by gen-ranges. Lik"\
                            f"ely caused by closely spaced or multimeric ints"
                for seg in segments:
                    idx += 1
                    if type(seg)==list:
                        if idx==0:
                            int_rng, gen_right = segments[:2]
                            if type(gen_right)==tuple:
                                get_junction_delta(False, gen_right[0], int_rng)
                            else:
                                read_ct -= 1
                                ilog.vlprint(fail_msg, 5)
                        elif idx==len(segments)-1:
                            gen_left, int_rng = segments[idx-1:]
                            if type(gen_left)==tuple:
                                get_junction_delta(gen_left[0], False, int_rng)
                            else:
                                read_ct -= 1
                                ilog.vlprint(fail_msg, 5)
                        else:
                            gen_left, int_rng, gen_right = segments[idx-1:idx+2]
                            if all([type(gen_left)==tuple, type(gen_right)==tuple]):
                                get_junction_delta(gen_left[0], gen_right[0], int_rng)
                            else:
                                read_ct -= 1
                                ilog.vlprint(fail_msg, 5)
            elif len(rtr.read_start_end_chr_start_end)==2 and len(cro.orig_alignment)==1:
                # ilog.vlprint(f"{rtr.sstrand}, {rtr.cand_read_obj.slen}", 10)
                read_ct += 1
                gen_cov_left, gen_cov_right =  gen_cov[0][0], gen_cov[1][0]
                inv = rtr.sstrand[0]=="minus"
                get_junction_delta(gen_cov_left, gen_cov_right, int_cov_rng, inv=inv)
            elif len(rtr.read_start_end_chr_start_end)==1 and len(cro.orig_alignment)==1:
                read_ct += 1
                # ilog.vlprint(f"{rtr.sstrand}, {rtr.cand_read_obj.slen}", 10)
                gen_cov = gen_cov[0][0]
                if rtr.flank_side == "left":
                    if rtr.sstrand[0] == "minus":
                        get_junction_delta(False, gen_cov, int_cov_rng, inv=True)
                    else:
                        get_junction_delta(gen_cov, False, int_cov_rng)
                elif rtr.flank_side == "right":
                    if rtr.sstrand[0] == "minus":
                        get_junction_delta(gen_cov, False, int_cov_rng, inv=True)
                    else:
                        get_junction_delta(False, gen_cov, int_cov_rng)
                else:
                    read_ct -= 1
                    ilog.vlprint(
                        f"Flank side for Ref alignment of {rtr.READ} n.a.", 5
                        )
            else:
                ilog.vlprint(f"{rtr.READ} not categorized in calc_seq_gain", 7)
                pass

        gains = {"left": gains_left, "right": gains_right}
        for k in gains:
            gain = gains[k]
            if len(gain) > 1: # quantiles reuires at least 2 data points
                gain = [v if v>=0 else 0 for v in gain]
                quants = stats.quantiles(gain, n=5, method="inclusive")
                g_clip = [v for v in gain if v<quants[-1]]
                g_clip = g_clip if len(g_clip)/len(gain)>=0.66 else gain
                mf = max(Counter(gain).items(), key=lambda x: x[1])[0]
                gain_info[k]["most_freq"] = mf
                gain_info[k]["min"] = min(gain)
                gain_info[k]["median"] = informed_round(stats.median(g_clip), mf)
        ilog.vlprint(
            f"Processed {read_ct}/{len(self.supp_als)} supp. alignments of {self.name}", 8
            )
        self.seq_gain_info = gain_info

    def rename(self):
        blast_info = [self.chrom, self.blast_loc]
        mm2_info = [self.mm2_chrom, self.mm2_loc]
        all_info = [self.chrom, self.blast_loc, self.mm2_chrom, self.mm2_loc]
        if all(all_info):
            self.approx_prec_loc()
            self.name = self.chrom + "_" + str(self.blast_loc)
            self.approx_loc = self.blast_loc
        elif all(blast_info) and not all(mm2_info):
            self.approx_prec_loc()
            self.name = self.chrom + "_" + str(self.blast_loc)
            self.approx_loc = self.blast_loc
        elif all(mm2_info) and not all(blast_info):
            self.approx_prec_loc()
            self.name = self.mm2_chrom + "_" + str(self.mm2_loc)
            self.approx_loc = self.mm2_loc
        elif not any(all_info):
            self.name = "lost_all_support_info"
        else:
            pass

    def re_calc_attributes(self, rename=False):
        re_calc = {
                "approx_loc": self.approx_prec_loc,
                "achieved_ident_matches": self.calc_achieved_ident_matches,
                "max_aims": self.calc_max_achievable_ident_matches,
                "ident_match_ratio": self.calc_imr,
                "read_names": self.get_blast_readnames,
                "total_support": self.calc_total_support
                }
        for attr in re_calc:
             self.__dict__[attr] = re_calc[attr]()
        if rename:
            self.name = self.chrom + "_" + str(self.approx_loc)

    def merge_attributes(self, other):
        for attr in self.__dict__:
            if type(self.__dict__[attr]) == list:
                self.__dict__[attr].extend(other.__dict__[attr])
            elif type(self.__dict__[attr]) == dict:
                self.__dict__[attr].update(other.__dict__[attr])
            else:
                pass
        self.re_calc_attributes()

    def merge(self, others_set):
        if not self.absorbed:
            self_stats = f"PotIntSite {self.name} with {len(self.total_support)} supporting alignments"
            for other in others_set:
                other_stats = f"PotIntSite {other.name} with {len(other.total_support)} supp. alignments"
                # self.supp_als.extend(other.supp_als)
                # self.minimap_supp.update(other.minimap_supp)
                # self.read_names.extend(other.read_names)
                self.merge_attributes(other)
                other.absorbed = True
                ilog.vlprint(
                    f"Merged {other_stats} into {self.name}",
                    2, alt_log="site_specs.log"
                    )
            # self.approx_loc = self.approx_prec_loc()
            # self.name = self.chrom + "_" + str(self.approx_loc)
            self.rename()
            ilog.vlprint(
                f"{len(others_set)} site/s was/were merged into {self_stats} "
                f"generating new PotIntSite {self.name} with "
                f"{len(self.total_support)} supporting alignments", 
                2, alt_log="site_specs.log"
                )
    
    def to_json(self):
        return json.dumps(
            {self.name: self.supp_als},
            default=lambda obj: obj.to_json(),
            indent=2
            )
    
    def __repr__(self) -> str:
        return f"<il_classes.PotIntSite object {self.name}>"
