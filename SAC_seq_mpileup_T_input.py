import sys

for line in sys.stdin:
    info = line.split()
    counts = int(info[3])
    Forward_seq_mutation = [0, 0, 0, 0, 0, 0, 0, 0] # match, skip, mismatch_A, mismatch_G, mismatch_C, mismatch_T, INDEL, INSERT
    start = 0
    stop = 0
    number_str = 0
    jump = 0
    if info[2] in ["T","t"] and counts >= 20:
        for i in info[4]:
            number_str += 1
            if jump != 0:
                jump -= 1
                continue
            if i == '.':
                Forward_seq_mutation[0] += 1
                continue
            if i == ',':
                Forward_seq_mutation[0] += 1
                continue
            if i == '>':
                Forward_seq_mutation[1] += 1
                continue
            if i == '<':
                Forward_seq_mutation[1] += 1
                continue
            if i == 'A':
                Forward_seq_mutation[2] += 1
                continue
            if i == 'a':
                Forward_seq_mutation[2] += 1
                continue
            if i == 'G':
                Forward_seq_mutation[3] += 1
                continue
            if i == 'g':
                Forward_seq_mutation[3] += 1
                continue
            if i == 'C':
                Forward_seq_mutation[4] += 1
                continue
            if i == 'c':
                Forward_seq_mutation[4] += 1
                continue
            if i == 'T':
                Forward_seq_mutation[5] += 1
                continue
            if i == 't':
                Forward_seq_mutation[5] += 1
                continue
            if i == '-':
                jump = int(info[4][number_str]) + 1
                Forward_seq_mutation[6] += 1
                continue
            if i == '+':
                jump = int(info[4][number_str]) + 1
                Forward_seq_mutation[7] += 1
                continue
            if i == '^':
                jump = 1
                start += 1
                continue
            if i == '$':
                stop += 1
                continue
        total_counts = Forward_seq_mutation[0] + Forward_seq_mutation[2] + Forward_seq_mutation[3] + Forward_seq_mutation[4] \
                       + Forward_seq_mutation[5] + Forward_seq_mutation[6] + Forward_seq_mutation[7]
        mutation_counts = Forward_seq_mutation[2] + Forward_seq_mutation[3] \
                       + Forward_seq_mutation[4] + Forward_seq_mutation[6] + Forward_seq_mutation[7]
        try:
            mutation_ratio = round(mutation_counts / total_counts, 4)
            if total_counts >= 20 and mutation_ratio <= 0.025:
                print(f'{info[0]}\t{info[1]}\t-\t{mutation_ratio}\t{total_counts}\t{mutation_counts}\tmutation_T/G/C={Forward_seq_mutation[2]}:{Forward_seq_mutation[4]}:{Forward_seq_mutation[3]}\t{stop}')
        except:
            0
