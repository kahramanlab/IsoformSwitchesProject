#!/usr/bin/env python3
"""
Calculates ClassIII peptide evidence for a cMDT
"""

import sys
import os
import argparse
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import re


class MDTAnalyzer:
    def __init__(self):
        self.verbose = False
        self.diann_sample_column_no = 1
        self.diann_gene_column_no = 5
        self.diann_peptide_abundance_column = 10
        self.diann_peptide_column_no = 14
        self.diann_qvalue_column_no = 17
        self.diann_pg_qvalue_column_no = 21
        
        self.cMDT_sample_column_no = 0
        self.cMDT_gene_column_no = 1

        self.max_qvalue = 0.1
        
        self.min_peptide_count = 5

        # File paths
        self.diann_file = None
        self.cmdt_file = None

    def parse_arguments(self):
        parser = argparse.ArgumentParser(description='Most Dominant Transcript Switch Analysis')

        parser.add_argument('--diann_file', required=True,
                            help='Peptide abundance file as analyzed by DIANN (e.g. AML/Main_Report.tsv)')
        parser.add_argument('--cmdt_file', required=True,
                            help='cMDT file (e.g. cMDT_peptides1+2_AML.tsv or SupplementaryFile1_AML.tsv)')
        parser.add_argument('--min_peptide_count', type=int, required=True,
                            help='minimum number of times a peptide needs to have been detected in order to check for significance')
        parser.add_argument('--max_qvalue', type=float, required=True,
                            help='maximum q-value')
        parser.add_argument('--verbose', action='store_true',
                            help='Print additional information')

        args = parser.parse_args()

        # Set attributes from arguments
        self.diann_file = args.diann_file
        self.cmdt_file = args.cmdt_file
        self.min_peptide_count = args.min_peptide_count
        self.max_qvalue = args.max_qvalue
        self.verbose = args.verbose

    def read_diann_file(self):
        """Read DIANN file"""

        diann = {}
        if self.verbose:
            print(f"Reading {self.diann_file} ... ", end="", file=sys.stderr)
        with open(self.diann_file, 'rt') as f:
            for line in f:
                if line.startswith('File.Name'):
                    continue
                line = line.strip()
                parts = line.split('\t')
                sample = re.sub("-.*", "", parts[self.diann_sample_column_no])
                gene = parts[self.diann_gene_column_no]
                peptide = parts[self.diann_peptide_column_no]
                qvalue = float(parts[self.diann_qvalue_column_no])
                pg_qvalue = float(parts[self.diann_pg_qvalue_column_no])
                if re.search("[0-9]", parts[self.diann_peptide_abundance_column]):
                    abundance = float(parts[self.diann_peptide_abundance_column])
                else:
                    if re.search("[0-9]", parts[self.diann_peptide_abundance_column-3]):
                        abundance = float(parts[self.diann_peptide_abundance_column-3])
                    else:
                        abundance = 0
                
                if(qvalue < self.max_qvalue and pg_qvalue < self.max_qvalue):
                    if gene not in diann:
                        diann[gene] = {}
                    if peptide not in diann[gene]:
                        diann[gene][peptide] = {}
                    if sample not in diann[gene][peptide]:
                        diann[gene][peptide][sample] = {}
                    diann[gene][peptide][sample] = abundance
        if self.verbose:
            print("Done.", file=sys.stderr)
        return diann

    def read_cMDT_file(self):
        """Read cMDT file"""

        cMDT = {}
        if self.verbose:
            print(f"Reading  {self.cmdt_file}  ... ", end="", file=sys.stderr)   
        with open(self.cmdt_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                sample = os.path.basename(parts[self.cMDT_sample_column_no])
                gene = parts[self.cMDT_gene_column_no]
                # only consider TuPro samples
                if re.search("^TP", parts[self.cMDT_sample_column_no]):
                    if gene not in cMDT:
                        cMDT[gene] = {}
                    cMDT[gene][sample] = line

        if self.verbose:
            print("Done.", file=sys.stderr)
        return cMDT 
    
    def run_analysis(self, diann, cMDT):
        if self.verbose:
            print(f"Running analysis ... ", end="", file=sys.stderr)

        for gene in cMDT.keys():
            # do this analysis only for cMDTs that have been detected at least min_peptide_count samples

            if gene in diann.keys():
                if len(cMDT[gene]) >= self.min_peptide_count:
                    pvalues = []
                    peptides = []
                    for peptide in diann[gene].keys():
                        if len(diann[gene][peptide]) >= self.min_peptide_count:
                            cMDT_abundances = []
                            xcMDT_abundances = []                    
                            for sample in diann[gene][peptide].keys():
                                abundance = diann[gene][peptide][sample]
                                if sample in cMDT[gene]:
                                    cMDT_abundances.append(abundance)
                                else:
                                    xcMDT_abundances.append(abundance)
                            
                            if len(cMDT_abundances) > 1 and len(xcMDT_abundances) > 1:
                                t_stat, pvalue = stats.ttest_ind(cMDT_abundances, xcMDT_abundances, equal_var=False) # use Welch's t-test to accomodate difference in cMDT vs xcMDT size differences and difference in abundance due to MDT
                                pvalues.append(pvalue)
                                peptides.append(peptide)
                    if len(pvalues)>0:
                        rejected, qvalues = fdrcorrection(pvalues, method='indep', alpha=self.max_qvalue)
                    else:
                        qvalues = []
                    hits = 0
                    hit_peptides = []
                    hit_qvalues = []
                    nonhit = 0
                    nonhit_peptides = []
                    nonhit_qvalues = []
                    i = -1
                    for q in qvalues:
                        i += 1
                        if q <= self.max_qvalue:
                            hits += 1
                            hit_peptides.append(peptides[i])
                            hit_qvalues.append(q)
                        else:
                            nonhit += 1
                            nonhit_peptides.append(peptides[i])
                            nonhit_qvalues.append(q)
            
                    for sample in cMDT[gene].keys():
                        if hits > 0 and nonhit > 0:
                            # print out peptide at the end of the line. Remove if needed the nan tag.
                            peptides = ""
                            for peptide in hit_peptides:
                                peptides += "," + peptide + "(ClassIII)"
                            
                            line = cMDT[gene][sample]
                            if re.search("nan$", line):
                                peptides = re.sub("^,", "", peptides)
                                print(re.sub("nan$", peptides, cMDT[gene][sample]))
                            else:
                                print(cMDT[gene][sample] + peptides)
                        else:
                            print(cMDT[gene][sample])
                else:
                    for sample in cMDT[gene].keys():
                        print(cMDT[gene][sample])

            else:
                for sample in cMDT[gene].keys():
                    print(cMDT[gene][sample])

        if self.verbose:
            print("Done.", file=sys.stderr)


def main():
    """Main function"""
    analyzer = MDTAnalyzer()
    analyzer.parse_arguments()
    cMDT_data = analyzer.read_cMDT_file()
    diann_data = analyzer.read_diann_file()
    analyzer.run_analysis(cMDT=cMDT_data, diann=diann_data)

if __name__ == "__main__":
    main()
