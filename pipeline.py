import subprocess
import glob
import os
import ntpath
import statistics
from Bio import SeqIO
from Bio import AlignIO

############################################################################################################################################
def BuildingTrees(myInputBT):
    aligned_tmp_file = "/tmp/Aligned.fasta"
    trimmed_tmp_file = "/tmp/Trimmed.fasta"

    #Sequences length standard deviation
    list_sd = list()
    for rec in SeqIO.parse(myInputBT, 'fasta'):
        seqLen = len(rec)
        list_sd.append(int(seqLen))
    stdev = statistics.stdev(list_sd)

# 1. Alignment
    cmd = ['mafft', '--maxiterate', '1000', '--globalpair', myInputBT]
    with open(aligned_tmp_file, 'w+') as f:
        p1 = subprocess.Popen(cmd, stdout=f,
                              stderr=subprocess.DEVNULL)
        p1.communicate()
    
# 2. Trimming
    cmd = ['trimal', '-automated1', '-in', aligned_tmp_file]
    with open(trimmed_tmp_file, 'w+') as f:
        p2 = subprocess.Popen(cmd, stdout=f)
        p2.communicate()

    #Trimmed alignment length
    alignment = AlignIO.read(trimmed_tmp_file, 'fasta')
    trimmed_length = alignment.get_alignment_length()

# 3. Tree construction for multithreading add '-T', '8'
    cmd = ['raxml', '-p', '12345', '-m', 'PROTGAMMAWAG', '-#', '100', '-s', trimmed_tmp_file, '-f', 'a', '-x', '12345',
           '-n', ntpath.basename(myInputBT), '-o', 'Drosophila_melanogaster']
    p3 = subprocess.Popen(cmd)
    p3.communicate()

# 4. Tree certainty for multithreading add '-T', '8'
    cmd = ['raxml', '-b', '12345', '-m', 'PROTGAMMAWAG', '-#', '100', '-f', 'i',
           '-n', 'TC', '-z', 'RAxML_bootstrap.' + ntpath.basename(myInputBT), '-t', 'RAxML_BestTree.' + ntpath.basename(myInputBT), '-L', 'MR']
    p4 = subprocess.Popen(cmd)
    p4.communicate()

# 5. Parse TC file and extract relative tree certainty
    with open("RAxML_info.TC") as tc_file:
        for line in tc_file.readlines():
            if line.startswith("Relative tree certainty for this tree:"):
                tree_certainty = line.split(" ")[-1]
    

# 6. Files and output
    file_r = open(myInputBT + '.result','w')

    file_r.write("Standard deviation of protein length: ")
    file_r.write(str(stdev))
    file_r.write("\n")
    file_r.write("Trimmed alignment length: ")
    file_r.write(str(trimmed_length))
    file_r.close()

    os.remove(aligned_tmp_file)
    os.remove(trimmed_tmp_file)

    return(stdev, trimmed_length, tree_certainty)

############################################################################################################################################

def raw_cur_finder(folder):
    raw_found = False
    cur_found = False
    
    try:
        raw_path = glob.glob(folder + os.path.sep + folder + "_RAW" + os.path.sep + folder + "-*RAW*.fasta")
        if os.path.exists(raw_path[0]):
            raw_found = True   
    except:
        print("Raw file not found!")
    
    try:
        cur_path = glob.glob(folder + os.path.sep + folder + "_CUR" + os.path.sep + folder + "-*CUR*.fasta")
        if os.path.exists(cur_path[0]):
            cur_found = True 
    except:
        print("Curated file not found!")

    return(raw_path, raw_found, cur_path, cur_found)

############################################################################################################################################

def main(in_folder, gene_list):
    with open("metrics_out.txt", "w+") as metrics_out:
        metrics_out.write("gene" + "\t " + "file" + "\t " + "stdev" + "\t" + "trim_len" + "\t" + "relative_TC" + "\n")

        print("Genes to process:", gene_list)
        os.chdir(in_folder)

        for gene in gene_list:
            print(gene)

            raw_cur = raw_cur_finder(gene)
            if raw_cur[1] and raw_cur[3] == True:
                path_to_raw = raw_cur[0][0]
                path_to_cur = raw_cur[2][0]

                results_raw = BuildingTrees(path_to_raw)
                metrics_out.write("".join([str(gene), "\t raw \t", str(results_raw[0]), "\t", str(results_raw[1]), "\t", str(results_raw[2]), "\n"]))

                results_cur = BuildingTrees(path_to_cur)
                metrics_out.write("".join([str(gene), "\t curated \t", str(results_cur[0]), "\t", str(results_cur[1]), "\t", str(results_cur[2]), "\n"]))

            print("Gene: ", gene, "processed!")

    metrics_out.close()

main("data", ["STAT1", "UPD"])
