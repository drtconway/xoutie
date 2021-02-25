"""xoutie - detect outliers in count data

Usage:
    xoutie count [options] <reference> <annotation> <bam>...
    xoutie model [options] <counts-file>

Arguments:
    <annotation>                a gtf file with gene/transcript/exon annotations.
    <bam>                       a bam file to process to extract counts.
    <counts-file>               a tab separated file with raw counts for features.

Options:
    -h, --help                  show help message
    -o FILE, --output FILE      output filename [default: -]
    --context-name NAME         column name for the context (group).
    --category-name NAME        column name of the category of the count
    --first-sample-column N     the first column (counting from 1) for sample counts
                                (preceeding columns are assumed to be metadata).
"""
from docopt import docopt
from intervaltree import Interval, IntervalTree
from gtfparse import read_gtf
import pandas

import gzip
import subprocess

__version__ = "0.0.1"

def smart_open(filename, mode="rt"):
    if filename == "-":
        if 'r' in mode:
            return sys.stdin
        elif 'w' in mode:
            return sys.stdout
        else:
            raise SmartOpenError(f"unable to open '-' in mode '{mode}'")

    if filename.endswith(".gz"):
        return gzip.open(filename, mode)

    return open(filename, mode)

class atoms(object):
    def __init__(self, transform = None):
        self.index = {}
        self.transform = transform

    def __getitem__(self, thing):
        if thing not in self.index:
            if self.transform is None:
                self.index[thing] = thing
            else:
                self.index[thing] = self.transform(thing)
        return self.index[thing]

class pileup_decoder(object):
    def __init__(self):
        self.ref_chars = [0 for i in range(256)]
        self.ref_chars[b'.'[0]] = 1
        self.ref_chars[b','[0]] = 1

        self.alt_chars = [0 for i in range(256)]
        self.alt_chars[b'A'[0]] = 1
        self.alt_chars[b'a'[0]] = 1
        self.alt_chars[b'C'[0]] = 1
        self.alt_chars[b'c'[0]] = 1
        self.alt_chars[b'G'[0]] = 1
        self.alt_chars[b'g'[0]] = 1
        self.alt_chars[b'T'[0]] = 1
        self.alt_chars[b't'[0]] = 1

        self.skp_chars = [0 for i in range(256)]
        self.skp_chars[b'<'[0]] = 1
        self.skp_chars[b'>'[0]] = 1

        self.skips = [1 for i in range(256)]
        self.skips[b'^'[0]] = 2

    def decode(self, pile : str):
        ref_count = 0
        alt_count = 0
        skp_count = 0
        i = 0
        while i < len(pile):
            c = pile[i]
            ref_count += self.ref_chars[c]
            alt_count += self.alt_chars[c]
            skp_count += self.skp_chars[c]
            i += self.skips[c]
        return (ref_count, alt_count, skp_count)

def estimate_gamma(xs, biasCorrection = False):
    n = len(xs)
    sx = sum(xs)
    slx = sum([math.log(x) for x in xs])
    sxlx = sum([x*math.log(x) for x in xs])
    kHat = n*sx / (n*sxlx - slx*sx)
    thetaHat = (n*sxlx - slx*sx) / (n*n)
    if biasCorrection:
        thetaHat *= n / (n - 1.0)
        kHat1 = 1.0 + kHat
        kHat -= (3*kHat - 2.0/3.0 * (kHat/kHat1) - 4.0/5.0 * (kHat/(kHat1*kHat1)))/n
    return (kHat, thetaHat)

def counts(bam_filename : str, genome_fasta : str):
    p = subprocess.Popen(['samtools', 'mpileup', '-f', genome_fasta, bam_filename], stdout=subprocess.PIPE)
    chroms = atoms(lambda s: s.decode('ascii'))
    dec = pileup_decoder()
    for l in p.stdout:
        t = l.split()
        chrom = chroms[t[0]]
        pos = int(t[1])
        cts = dec.decode(t[4])
        yield (chrom, pos, cts)

class deep_dict(object):
    def __init__(self):
        self.idx = {}

    def __getitem__(self, ks):
        x = self.idx
        for k in ks:
            x = x[k]
        return x

    def __setitem__(self, ks, val):
        x = self.idx
        for k in ks[:-1]:
            if k not in x:
                x[k] = {}
            x = x[k]
        x[ks[-1]] = val

class segement_index(object):
    def __init__(self, annot):
        wanted = {'gene', 'transcript', 'exon'}

        self.gene_names = {}

        pos_sets = {}
        for ftr in annot.iterrows():
            ftr_type = ftr[1]['feature']
            if ftr_type not in wanted:
                continue
            chrom = ftr[1]['seqname']
            begin = ftr[1]['start']
            end = ftr[1]['end'] + 1

            if ftr_type == 'exon':
                gene_id = ftr[1]['gene_id']
                gene_nm = ftr[1]['gene_name']
                self.gene_names[gene_id] = (gene_nm, chrom)

                if chrom not in pos_sets:
                    pos_sets[chrom] = {}
                if gene_id not in pos_sets[chrom]:
                    pos_sets[chrom][gene_id] = IntervalTree()
                begMin = begin
                endMax = end
                ols = pos_sets[chrom][gene_id][begin:end]
                if len(ols) > 0:
                    for ivl in ols:
                        ivlBeg = ivl[0]
                        ivlEnd = ivl[1]
                        pos_sets[chrom][gene_id].remove(ivl)
                        begMin = min(begMin, ivlBeg)
                        endMax = max(endMax, ivlEnd)
                pos_sets[chrom][gene_id].add(Interval(begMin, endMax))

        self.segs = {}
        self.gene_segs = {}
        self.seg_idx = {}
        for chrom in pos_sets.keys():
            self.segs[chrom] = IntervalTree()
            for gene_id in pos_sets[chrom]:
                ps = set([])
                for ivl in pos_sets[chrom][gene_id]:
                    ps.add(ivl[0])
                    ps.add(ivl[1])
                ps = sorted(ps)
                self.gene_segs[gene_id] = len(ps) - 1
                self.seg_idx[gene_id] = []
                for i in range(1, len(ps)):
                    st = ps[i - 1]
                    en = ps[i]
                    self.segs[chrom][st:en] = (gene_id, i-1)
                    self.seg_idx[gene_id].append((st, en))

    def length(self, gene_id):
        return self.gene_segs[gene_id]

    def lookup(self, chrom, pos):
        if chrom not in self.segs:
            return set()
        return self.segs[chrom][pos]

    def gene_name(self, gene_id):
        return self.gene_names[gene_id]

    def gene_segments(self, gene_id):
        return self.seg_idx[gene_id]

def scanbam(opts):
    annot = read_gtf(opts["<annotation>"])
    idx = segement_index(annot)

    res = {}
    for bamName in opts["<bam>"]:
        print(bamName)
        print(opts["<reference>"])
        for t in counts(bamName, opts["<reference>"]):
            chrom = t[0]
            pos = t[1]
            cts = t[2]
            stab = idx.lookup(chrom, pos)
            if len(stab) == 0:
                continue

            cov = cts[0] + cts[1]

            for v in stab:
                w = v.data
                gene = w[0]
                seg = w[1]
                if gene not in res:
                    res[gene] = [{} for i in range(idx.length(gene))]
                if cov not in res[gene][seg]:
                    res[gene][seg][cov] = 0
                res[gene][seg][cov] += 1
    with smart_open(opts['--output'], "wt") as out:
        print(f"name\tid\tsegment\tchrom\tbegin\tend\tcoverage\tfrequency", file=out)
        for gene_id in sorted(res.keys()):
            (gene_nm, chrom) = idx.gene_name(gene_id)
            gene_segs = idx.gene_segments(gene_id)
            for i in range(len(res[gene_id])):
                seg = gene_segs[i]
                for (c,f) in sorted(res[gene_id][i].items()):
                    print(f"{gene_nm}\t{gene_id}\t{i}\t{chrom}\t{seg[0]}\t{seg[1]}\t{c}\t{f}", file=out)

def consmodel(opts):
    dat = pandas.read_table(opts['<counts-file>'])

def main():
    opts = docopt(__doc__, version=__version__)
    if opts["count"]:
        return scanbam(opts)
    if opts["model"]:
        return consmodel(opts)

if __name__ == "__main__":
    main()
