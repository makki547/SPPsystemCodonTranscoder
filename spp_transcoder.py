import argparse
import pandas as pd
import re
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
import sys
import os

CODON_USAGE_DB = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'ecoli_b_codon_usage.txt')

def parse_args():

    parser = argparse.ArgumentParser(description = 'Codon transcoder to make ACA-less sequence for SPP-system')


    parser.add_argument('fasta_input', type = str, help = 'FASTA file for the target')
    parser.add_argument('fasta_output', type = str, help = 'FASTA file for the output without ACA triplet.')
    parser.add_argument('-c', '--check_rarecodon', action = 'store_true', help = 'Check whether any rare-codons are introduced by transcoding')
    parser.add_argument('-d', '--codon_usage_db', type = str, default = None, help = 'Codon-usage DB data of Kazusa DNA Res. Inst., Format "A style like CodonFrequency output in GCG Wisconsin Package" is required.')
    parser.add_argument('-t', '--threshold', type = float, default = 9.0, help = 'Threshold for frequency of rare codon, 9.0 %% is default, considered from that of GGA codon that is enriched in Rosetta-strain.')

    return parser.parse_args()


    


def transcode_to_aca_less(seq):
    
    # Following trancoding rules are from Suzuki, Mao and Inouye, Nature Protocols, 2, 1802-1810(2007)
    def transcode_frame3_xxacax(double_frame):
        aa = str(Seq(double_frame[0]).translate())
        
        if aa in list('TRG'):
            return ('{}C'.format(double_frame[0][0:2]), 'CA{}'.format(double_frame[1][2:3]))

        else:
            return ('{}G'.format(double_frame[0][0:2]), 'CA{}'.format(double_frame[1][2:3]))
        
    
    frame1_rule = (re.compile(r'^ACA'), lambda x : ('ACC', x[1]) )
    
    frame2_rule = (re.compile(r'^[ATCG]ACA'), lambda x: ('{}AT'.format(x[0][0:1]), 'A{}'.format(x[1][1:3]) ) )
    
    frame3_rule1 = (re.compile(r'^ATACA'), lambda x: ('ATT', 'CA{}'.format(x[1][2:3])))
    frame3_rule2 = (re.compile(r'^AGACA'), lambda x: ('CGC', 'CA{}'.format(x[1][2:3])))
    frame3_rule3 = (re.compile(r'^TTACA'), lambda x: ('CTG', 'CA{}'.format(x[1][2:3])))
    frame3_rule4 = (re.compile(r'^T[GA]ACA'), lambda x: ('TAG', 'CA{}'.format(x[1][2:3])))
    frame3_rule5 = (re.compile(r'^[ATCG]{2}ACA'), transcode_frame3_xxacax)
    
    translate_rules = [
        frame1_rule,
        frame2_rule,
        frame3_rule1, frame3_rule2, frame3_rule3, frame3_rule4, frame3_rule5
    ]
    
    codons = [seq[i: i+3] for i in range(0, len(seq), 3)] + ['---']
    codons.reverse()
    ncodons = len(codons) - 1
    
    transcoded = ''
    for i in range(ncodons):
        
        codon1 = codons.pop()
        codon2 = codons[-1]
        
        double_frame = [codon1, codon2]
        
        subseq = codon1 + codon2
        for rule, proc in translate_rules:
            
            if rule.match(subseq):
                cc = proc(double_frame)
                transcoded += cc[0]
                codons[-1] = cc[1]
                break
        else:
            transcoded += codon1
    
    if Seq(transcoded).count('ACA') != 0:
        raise Exception('ACAs were found, transcoding logic should be wrong.')
        
    if Seq(transcoded).translate() != Seq(seq).translate():
        raise Exception ('Protein seq. from trancoded DNA seq is different from the original. Transcoding logic should be wrong.')
        
    return transcoded

def check_rarecodon(name, seq_original, seq_acaless, df_rarecodon):

    seq_orig = str(seq_original)
    seq_aca = str(seq_acaless)
    
    codons_orig = [seq_orig[i: i+3] for i in range(0, len(seq_orig), 3)]
    codons_acaless = [seq_aca[i: i+3] for i in range(0, len(seq_aca), 3)]

    rare_codons_orig = np.array(list(map(lambda codon: df_rarecodon.loc[codon ,'Rare'], codons_orig)))
    rare_codons_acaless = np.array(list(map(lambda codon: df_rarecodon.loc[codon ,'Rare'], codons_acaless)))
    rare_codons_xor = np.logical_xor(rare_codons_orig, rare_codons_acaless)

    if not np.any(np.logical_and(rare_codons_xor, rare_codons_acaless)):
        sys.stdout.write('Any rare codons are not introduced by the transcoding.\n\n')
        return
    
    sys.stdout.write('Rare codons of {} (following offsets are zero-based index):\n'.format(name))
    for i, (codon, rare) in enumerate(zip(codons_acaless, rare_codons_xor)):

        if rare and rare_codons_acaless[i]:
            
            sys.stdout.write('Offset: DNA seq. {:>4d} (AA seq., {:>3d}), AA: {}, Codon: {}, Frequency: {:>4.1f} %'.format(
                i*3, i, df_rarecodon.loc[codon ,'AA'], codon, df_rarecodon.loc[codon ,'Frequency']))

            sys.stdout.write('\n')

    sys.stdout.write('\n')



if __name__ == '__main__':

    parsed_args = parse_args()


    if parsed_args.check_rarecodon:

        if parsed_args.codon_usage_db is None:
            db_path = CODON_USAGE_DB
        else:
            db_path = parsed_args.codon_usage_db
            
        df_rarecodon = pd.read_csv(db_path, sep = '\s+')
        df_rarecodon = df_rarecodon.set_index('Codon')
        df_rarecodon['Rare'] = df_rarecodon['Frequency'] < parsed_args.threshold

    with open(parsed_args.fasta_output, 'w') as fout:
    
        for record in SeqIO.parse(parsed_args.fasta_input, 'fasta'):

            sys.stderr.write('{} is now being transcoded...\n'.format(record.name))
            transcoded = transcode_to_aca_less(str(record.seq))

            if parsed_args.check_rarecodon:

                
                check_rarecodon(record.name, record.seq, transcoded, df_rarecodon)
                
            if record.description != '':
                record.description += ', '

            record.description += 'ACA-less'
            record.seq = Seq(transcoded)


            SeqIO.write(record ,fout, 'fasta')


    sys.stderr.write('\nPlease cite this article:\n')
    sys.stderr.write('\tSuzuki, Mao, and Inouye, Nature Protocols, 2, 1802-1810(2007), DOI:10.1038/nprot.2007.252\n')
    if parsed_args.check_rarecodon:
        sys.stderr.write('\tNakamura, Gojobori, and Ikemura, Nucl. Acids Res. 28, 292(2000), DOI:10.1093/nar/28.1.292\n')
