import sys


for row in sys.stdin:
    row_lst = row.strip("\n").split(".")
    lab = row_lst[1]
    K = lab.split("r")[0][1:]
    
    sys.stdout.write("{}\t{}\t{}\n".format(lab, K, row.strip("\n")))
