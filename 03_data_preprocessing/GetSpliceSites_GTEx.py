from pyensembl import EnsemblRelease
import pandas as pd

ensemblv83 = EnsemblRelease(83)


def get_exon_junctions(transcript_id):
    try:
        exon_junctions = []
        transcript = ensemblv83.transcript_by_id(transcript_id)
        len_utr5 = len(transcript.five_prime_utr_sequence)
        cds_intervals = transcript.coding_sequence_position_ranges
        if transcript.on_positive_strand:
            cds_intervals.sort()
        else:
            cds_intervals.sort(reverse=True)
        pointer = len_utr5
        for cds in cds_intervals:
            cds_length = cds[1] - cds[0] + 1
            pointer += cds_length
            if not cds == cds_intervals[-1]:
                exon_junctions.append(pointer)
        return len_utr5, exon_junctions
    except ValueError:
        return None, None


df = pd.read_excel("D:/Documents/Bioinformatik/Master/ML/project_1/data/delta_hl_gtex.xlsx")

f = open("D:/Documents/Bioinformatik/Master/ML/project_1/data/GTEx_exon_junctions.txt", "w")
f.write("GeneID\tUTR5_len\tExon_Junctions_In_Full_Sequence\n")

for index, row in df.iterrows():
    transcript_id = str(row[0]).split(".")[0]
    len_utr5, exon_junctions = get_exon_junctions(transcript_id)
    if len_utr5 is not None:
        f.write(transcript_id + "\t" + str(len_utr5) + "\t" + (",".join(str(ej) for ej in exon_junctions) + "\n"))

f.close()
