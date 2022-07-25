import pandas as pd
from Bio import SeqIO
from scipy.stats import zscore

from collections import defaultdict
import gzip
import re


GTF_HEADER = ['seqname', 'source', 'feature', 'start', 'end', 'score',
              'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA = re.compile(r'\s*,\s*')
R_KEYVALUE = re.compile(r'(\s+|\s*=\s*)')


def dataframe(filename):
    """Open an optionally gzipped GTF file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value


df = pd.ExcelFile('https://drive.google.com/uc?export=download&id=1SXVX6hkj1DJYOiathpGU1CJ2ClI3g_XO')
df_human = pd.read_excel(df, "human", skiprows=[0, 1])
df_human['zscore'] = zscore(df_human['half-life (PC1)'])

with open(
        'D:/Documents/Bioinformatik/Master/ML/project_1/data/Homo_sapiens.GRCh38.83.chosenTranscript.3pUTRs.fa') as fasta_file:  # Will close handle cleanly
    UTR3_identifiers = []
    UTR3_seqs = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        UTR3_identifiers.append(seq_record.id)
        UTR3_seqs.append(seq_record.seq)

with open(
        'D:/Documents/Bioinformatik/Master/ML/project_1/data/Homo_sapiens.GRCh38.83.chosenTranscript.5pUTRs.fa') as fasta_file:  # Will close handle cleanly
    UTR5_identifiers = []
    UTR5_seqs = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        UTR5_identifiers.append(seq_record.id)
        UTR5_seqs.append(seq_record.seq)

with open(
        'D:/Documents/Bioinformatik/Master/ML/project_1/data/Homo_sapiens.GRCh38.83.chosenTranscript.ORFs.fa') as fasta_file:  # Will close handle cleanly
    ORF_identifiers = []
    ORF_seqs = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        ORF_identifiers.append(seq_record.id)
        ORF_seqs.append(seq_record.seq)

d = {'geneID': UTR5_identifiers, 'UTR5_seqs': UTR5_seqs}
UTR5 = pd.DataFrame(data=d)

d = {'geneID': ORF_identifiers, 'ORF_seqs': ORF_seqs}
ORF = pd.DataFrame(data=d)

d = {'geneID': UTR3_identifiers, 'UTR3_seqs': UTR3_seqs}
UTR3 = pd.DataFrame(data=d)

halflife = df_human[["Ensembl Gene Id", "zscore"]]

seqs = pd.merge(pd.merge(UTR5, ORF, on='geneID'), UTR3, on='geneID')
sequences = pd.merge(halflife, seqs, right_on='geneID', left_on='Ensembl Gene Id')
sequences = sequences.drop(columns=["geneID"])
sequences = sequences.rename(columns={"Ensembl Gene Id": "geneID"})

rubbish = [d, UTR3, ORF, UTR5, UTR5_seqs, UTR3_seqs, ORF_seqs, UTR3_identifiers, UTR5_identifiers, ORF_identifiers, df,
           df_human]
del rubbish

seqs = sequences['UTR5_seqs'] + sequences['ORF_seqs'] + sequences['UTR3_seqs']
seqs = seqs.dropna()

df = dataframe("D:/Documents/Bioinformatik/Master/ML/project_1/data/Homo_sapiens.GRCh38.83.chosenTranscript.gtf")
df[["transcriptID", "geneID"]] = df["INFO1"].str.split(':', expand=True)
df = df.drop(columns=["INFO1"])

f = open("D:/Documents/Bioinformatik/Master/ML/project_1/data/Kelley_et_al_exon_junctions.txt", "w")
f.write("GeneID\tUTR5_len\tExon_Junctions_In_Full_Sequence\n")

for index, row in sequences.iterrows():

    utr5_len = len(row["UTR5_seqs"])

    tmp_junction_indices = []
    pointer = utr5_len

    df_trans = df.loc[(df['geneID'] == row["geneID"]) & (df['feature'] == "transcript")]
    df_cds = df.loc[(df['geneID'] == row["geneID"]) & (df['feature'] == "CDS")]
    df_cds["start"] = df_cds["start"].astype(int)

    start = int(df_trans["start"])
    end = int(df_trans["end"])
    if "-" in df_trans["strand"]:
        df_cds.sort_values(by=["start"], ascending=False)
    else:
        df_cds.sort_values(by=["start"])

    for i, cds in df_cds.iterrows():
        cds_start = int(cds["start"])
        cds_stop = int(cds["end"])

        len_exon = abs(cds_start - cds_stop) + 1

        pointer += len_exon
        if not i == df_cds.index[-1]:
            tmp_junction_indices.append(pointer)

    f.write(row["geneID"] + "\t" + str(utr5_len) + "\t" + (",".join(str(ej) for ej in tmp_junction_indices) + "\n"))

f.close()
