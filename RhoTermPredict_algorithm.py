import numpy as np
import re
from Bio import SeqIO
from Bio.Seq import Seq
import csv


def palindrome_finder(sequence, file=None, gc_genome=None, strand=None, score=None, calc="no"):
    n_rut = 0
    list_scores = []
    sequence = Seq(sequence)
    sequence_r = str(sequence.reverse_complement())
    _pattern_pause_site1 = "GG\D{8}[C,T]G"
    _pattern_pause_site2 = "C[G,A]\D{8}CC" 
    for k in range(4, 8):
        for i in range(150):
            _x1 = i
            _y1 = i + k
            if _y1 > len(sequence):
                break
            else:
                window = sequence[_x1:_y1]
                window = str(window)
                gc_content = 100*(window.count("C")+window.count("G"))/len(window)
                score_p = score
                sub_sequence_r = sequence_r
                _start = 0
                while True:
                    if re.search(window, sub_sequence_r):
                        _research = re.search(window, sub_sequence_r)
                        _positions = _research.span()
                        _x2 = _positions[0]
                        _y2 = _positions[1]
                        sub_sequence_r = sub_sequence_r[_y2:]
                        _x2 += _start
                        _y2 += _start
                        _start += _y2
                        if 4 <= len(sequence)-_y2-_y1+1 <= 8:
                            if calc == "yes":
                                loop = len(sequence)-_y2-_y1
                                file.write("\nPalindromic sequences found at coordinates ")
                                file.write(str(_x1))
                                file.write("-")
                                file.write(str(_y1-1))
                                file.write(" e ")
                                file.write(str(len(sequence)-_y2))
                                file.write("-")
                                file.write(str(len(sequence)-_x2-1))
                                file.write("(Sequence:   ")
                                file.write(str(window))
                                file.write(")")
                                score_p += 3
                                if gc_content > gc_genome + 20:
                                    score_p += 2
                                elif gc_content > gc_genome + 10:
                                    score_p += 1
                                if len(window) > 4:
                                    score_p += 1
                                if loop < 6:
                                    score_p += 1
                                if strand == 1:
                                    a = _x1-5
                                    b = len(sequence)-_x2+5
                                    seq_pause = str(sequence[a:b])
                                    if re.search(_pattern_pause_site1, seq_pause):
                                        score_p += 3
                                        file.write(" (PAUSE-CONSENSUS PRESENT)")
                                if strand == -1:
                                    a = _x1-5
                                    b = len(sequence)-_x2+5
                                    seq_pause = str(sequence[a:b])
                                    if re.search(_pattern_pause_site2, seq_pause):
                                        score_p += 3
                                        file.write(" (PAUSE-CONSENSUS PRESENT)")
                                file.write(" (SCORE: ")
                                file.write(str(score_p))
                                file.write(")")
                                list_scores.append(score_p)
                            else:
                                n_rut += 1
                    else:
                        break
    if calc == "yes":
        if len(list_scores) > 0:
            return np.max(list_scores)
        else:
            return 0
    elif calc == "no":
        return n_rut

                        

sequences_to_analyze = {}

print("RhoTermPredict is a genome-wide predictor of transcription Rho-dependent terminators "
      "in bacterial genomes. It analyzes both the strands.\n\n")
print("Input: Genome sequences file")
print("Output: a xlsx file containing Rho-dependent terminators coordinates and a txt file "
      "containing informations about them")
print("Genome file must be in fasta format")
file_genome = input("Enter the input genome file name: ")

try:
    for seq_record in SeqIO.parse(file_genome, "fasta"):
        sequences_to_analyze[seq_record.id] = str(seq_record.seq)

    num_sequences = len(sequences_to_analyze)
    print(f"Sequences to analyze: {num_sequences}\n\n"
          f" RhoTermPredict is working, please wait...")

    for sq in sequences_to_analyze:
        genome = sequences_to_analyze[sq]
        genome = genome.upper()
        gc_whole_genome = 100*(genome.count("G")+genome.count("C"))/len(genome)
        pattern1 = "C\D{11,13}C\D{11,13}C\D{11,13}C\D{11,13}C\D{11,13}C"
        pattern2 = "G\D{11,13}G\D{11,13}G\D{11,13}G\D{11,13}G\D{11,13}G"
        pattern_pause_site1 = "GG\D{8}[C,T]G"
        pattern_pause_site2 = "C[G,A]\D{8}CC"
        predictions = 0
        cg_value = []
        scores = []
        num = 1
        cod = 1

        predictions_file = open(f'predictions_coordinates_{sq}.csv', mode='a')
        writer = csv.writer(predictions_file, delimiter='\t')
        writer.writerow(['Region', 'Start RUT', 'End RUT', 'Strand'])
        p = open(f"info_about_predictions_{sq}.txt", "a")
        p.write("Sequences of predicted Rho-dependent terminators")

        # positive strand

        scale = 0
        for j in range(len(genome)):
            x1 = scale + j
            x2 = scale + j + 78
            if x2 > len(genome):
                break
            else:
                w = genome[x1:x2]
                if w.count("G") > 0:
                    numG = w.count("G")
                    c_over_g = w.count("C")/numG
                else:
                    numG = 1
                    c_over_g = (w.count("C")+1)/numG
                if c_over_g > 1:
                    if re.search(pattern1, w):
                        data = []
                        for g in range(50):
                            if x2+g <= len(genome):
                                w = genome[x1+g:x2+g]
                                if w.count("G") > 0:
                                    numG = w.count("G")
                                    c_over_g = w.count("C")/numG
                                else:
                                    numG = 1
                                    c_over_g = (w.count("C")+1)/numG
                                if c_over_g > 1:
                                    if re.search(pattern1, w):
                                        data.append(c_over_g)
                                    else:
                                        data.append(0)
                                else:
                                    data.append(0)
                            else:
                                break
                        maxP = np.argmax(data)
                        c_over_g = np.max(data)
                        x1 = x1 + maxP
                        x2 = x2 + maxP
                        scale = x2 - j - 1
                        s = genome[x2:x2+150]
                        w = genome[x1:x2]
                        score = 3
                        ctrl = palindrome_finder(s)
                        if ctrl > 0 or re.search(pattern_pause_site1, s):
                            writer.writerow([f'T{cod}', x1, x2, 'plus'])
                            num += 1
                            predictions += 1
                            scale = x2 + 150 - j - 1
                            cg_value.append(c_over_g)
                            if c_over_g > 2:
                                score += 3
                            elif c_over_g > 1.50:
                                score += 2
                            elif c_over_g > 1.25:
                                score += 1
                            p.write("\n\n\nPREDICTED REGION NUMBER ")
                            p.write(str(cod))
                            p.write(" (STRAND POSITIVE) ")
                            p.write("\n\nGenomic sequence of putative RUT site (Coordinates:  ")
                            p.write(str(x1))
                            p.write("-")
                            p.write(str(x2))
                            p.write(", c/g = ")
                            p.write(str(c_over_g))
                            p.write(" ):      ")
                            p.write(str(w))
                            p.write("\n\nThe 150 nt long genomic region immediately downstream is ")
                            p.write(str(s))
                            p.write("\n")
                            cod += 1
                            final_score = palindrome_finder(s, p, gc_whole_genome, 1, score, calc="yes")
                            start = 0
                            while True:
                                if re.search(pattern_pause_site1, s):
                                    research = re.search(pattern_pause_site1, s)
                                    positions = research.span()
                                    x2 = positions[0]
                                    y2 = positions[1]
                                    s = s[y2:]
                                    p.write("\n\nPAUSE-CONSENSUS present at the coordinates ")
                                    p.write(str(x2))
                                    p.write("-")
                                    p.write(str(y2))
                                    x2 += start
                                    y2 += start
                                else:
                                    break
                            if final_score == 0:
                                final_score = score + 3
                            scores.append(final_score)

        # negative strand
        scale = 0
        for j in range(len(genome)):
            x1 = scale + j
            x2 = scale + j + 78
            if x2 > len(genome):
                break
            else:
                w = genome[x1:x2]
                if w.count("C") > 0:
                    numC = w.count("C")
                    c_over_g = w.count("G")/numC
                else:
                    numC = 1
                    c_over_g = (w.count("G")+1)/numC
                if c_over_g > 1:
                    if re.search(pattern2, w):
                        data = []
                        for g in range(50):
                            if x2+g <= len(genome):
                                w = genome[x1+g:x2+g]
                                if w.count("C") > 0:
                                    numC = w.count("C")
                                    c_over_g = w.count("G")/numC
                                else:
                                    numC = 1
                                    c_over_g = (w.count("G")+1)/numC
                                if c_over_g > 1:
                                    if re.search(pattern2, w):
                                        data.append(c_over_g)
                                    else:
                                        data.append(0)
                                else:
                                    data.append(0)
                            else:
                                break
                        maxP = np.argmax(data)
                        c_over_g = np.max(data)
                        x1 = x1 + maxP
                        x2 = x2 + maxP
                        scale = x2 - j - 1
                        s = genome[x1-150:x1]
                        w = genome[x1:x2]
                        score = 3
                        ctrl = palindrome_finder(s)
                        if ctrl > 0 or re.search(pattern_pause_site2, s):
                            writer.writerow([f'T{cod}', x1, x2, 'minus'])
                            num += 1
                            predictions += 1
                            scale = x2 + 150 - j - 1
                            cg_value.append(c_over_g)
                            if c_over_g > 2:
                                score += 3
                            elif c_over_g > 1.50:
                                score += 2
                            elif c_over_g > 1.25:
                                score += 1
                            p.write("\n\n\nPREDICTED REGION NUMBER ")
                            p.write(str(cod))
                            p.write(" (STRAND NEGATIVE)")
                            p.write("\n\nGenomic sequence of putative RUT site (Coordinates:  ")
                            p.write(str(x1))
                            p.write("-")
                            p.write(str(x2))
                            p.write(", c/g = ")
                            p.write(str(c_over_g))
                            p.write(" ):      ")
                            p.write(str(w))
                            p.write("\n\nThe 150 nt long genomic region immediately downstream is ")
                            p.write(str(s))
                            p.write("\n")
                            cod += 1
                            final_score = palindrome_finder(s, p, gc_whole_genome, -1, score, calc="yes")
                            start = 0
                            while True:
                                if re.search(pattern_pause_site2, s):
                                    research = re.search(pattern_pause_site2, s)
                                    positions = research.span()
                                    x2 = positions[0]
                                    y2 = positions[1]
                                    s = s[y2:]
                                    p.write("\n\nPAUSE-CONSENSUS present at the coordinates ")
                                    p.write(str(x2))
                                    p.write("-")
                                    p.write(str(y2))
                                    x2 += start
                                    y2 += start
                                else:
                                    break
                            if final_score == 0:
                                final_score = score + 3
                            scores.append(final_score)

        p.write("\n\n\n\n\n\nTotal number of predicted Rho-dependent terminators: ")
        p.write(str(predictions))
        p.write("\n\n\nMean C/G content of predicted terminators: ")
        p.write(str(np.mean(cg_value)))
        p.write("\nStandard deviation: ")
        p.write(str(np.std(cg_value)))

        p.close()
        predictions_file.close()

    print("Work finished, see output files in the current directory")
    
except IOError:   
    print(f"File {file_genome} not existent in the current directory!")
