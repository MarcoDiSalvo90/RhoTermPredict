import numpy as np
import matplotlib.pyplot as plt
import re
from Bio import SeqIO
from Bio.Seq import Seq
import openpyxl

def palindrome_control(sequenza,GC_genoma,strand):
    palindr = 0
    sequenza = Seq(sequenza)
    sequenza_r = sequenza.reverse_complement()
    for k in range (4,8):
        for i in range (150):
            x1 = i
            y1 = i + k
            if y1 > len(sequenza):
                break
            else:
                window = sequenza[x1:y1]
                window = str(window)
                gc_content = 100*(window.count("C")+window.count("G"))/len(window)
                sequenza_r = str(sequenza_r)
                sub_sequenza_r= sequenza_r
                par = 0
                while True:
                    if re.search(window, sub_sequenza_r):
                        ricerca = re.search(window, sub_sequenza_r)
                        positions = ricerca.span()
                        x2 = positions[0]
                        y2 = positions[1]
                        sub_sequenza_r=sub_sequenza_r[y2:]
                        x2 += par
                        y2 += par
                        par += y2
                        if 4 <= len(sequenza)-y2-y1+1 <= 8:
                            palindr +=1
                            break
                    else:
                        break
    return palindr


def palindrome_finder(sequenza,file,GC_genoma,strand,score):
    lista_scores = []
    sequenza = Seq(sequenza)
    sequenza_r = sequenza.reverse_complement()
    pattern_pause_site1 = "GG\D{8}[C,T]G"
    pattern_pause_site2 = "C[G,A]\D{8}CC" 
    for k in range (4,8):
        for i in range (150):
            x1 = i
            y1 = i + k
            if y1 > len(sequenza):
                break
            else:
                window = sequenza[x1:y1]
                window = str(window)
                gc_content = 100*(window.count("C")+window.count("G"))/len(window)
                sequenza_r = str(sequenza_r)
                score_p = score
                sub_sequenza_r= sequenza_r
                par = 0
                while True:
                    if re.search(window, sub_sequenza_r):
                        ricerca = re.search(window, sub_sequenza_r)
                        positions = ricerca.span()
                        x2 = positions[0]
                        y2 = positions[1]
                        sub_sequenza_r=sub_sequenza_r[y2:]
                        x2 += par
                        y2 += par
                        par += y2
                        if 4 <= len(sequenza)-y2-y1+1 <= 8:
                            loop = len(sequenza)-y2 - y1
                            file.write("\nPalindromic sequences found at coordinates ")
                            file.write(str(x1))
                            file.write("-")
                            file.write(str(y1-1))
                            file.write(" e ")
                            file.write(str(len(sequenza)-y2))
                            file.write("-")
                            file.write(str(len(sequenza)-x2-1))
                            file.write("(Sequence:   ")
                            file.write(str(window))
                            file.write(")")
                            score_p += 3
                            if gc_content > GC_genoma+20:
                                score_p += 2
                            elif gc_content > GC_genoma+10:
                                score_p += 1
                            if len(window) > 4:
                                score_p += 1
                            if loop < 6:
                                score_p += 1
                            if strand == 1:
                                a = x1-5
                                b = len(sequenza)-x2+5
                                seq_pause = sequenza[a:b]
                                seq_pause = str(seq_pause)
                                if re.search(pattern_pause_site1,seq_pause):
                                    score_p += 3
                                    file.write(" (PAUSE-CONSENSUS PRESENT)")
                            if strand == -1:
                                a = x1-5
                                b = len(sequenza)-x2+5
                                seq_pause = sequenza[a:b]
                                seq_pause = str(seq_pause)
                                if re.search(pattern_pause_site2,seq_pause):
                                    score_p += 3
                                    file.write(" (PAUSE-CONSENSUS PRESENT)")
                            file.write(" (SCORE: ")
                            file.write(str(score_p))
                            file.write(")")
                            lista_scores.append(score_p)
                    else:
                        break
    if len(lista_scores) > 0:
        return np.max(lista_scores)
    else:
        return 0
                        
    

MC = []

print("RhoTermPredict is a genome-wide predictor of transcription Rho-dependent terminators in bacterial genomes. It analyzes both the strands.\n\n")
print("Input: Genome sequences file")
print("Output: a xlsx file containing Rho-dependent terminators coordinates and a txt file containing informations about them")
print("Genome file must be in fasta format")
File_genome = input("Enter the input genome file name: ")

try:
    for seq_record in SeqIO.parse(File_genome, "fasta"):
        MC.append(str(seq_record.seq))
    print ("RhoTermPredict is working, please wait...")
    genome = MC[0] 
    genome = genome.upper()
    GC_whole_genome = 100*(genome.count("G")+genome.count("C"))/len(genome)
    pattern1 = "C\D{11,13}C\D{11,13}C\D{11,13}C\D{11,13}C\D{11,13}C"
    pattern2 = "G\D{11,13}G\D{11,13}G\D{11,13}G\D{11,13}G\D{11,13}G"
    pattern_pause_site1 = "GG\D{8}[C,T]G"
    pattern_pause_site2 = "C[G,A]\D{8}CC" 
    predictions = 0
    cg_value = []
    punteggi = []
    num = 1
    cod = 1
    
    nuovo_file = openpyxl.Workbook()
    sheet = nuovo_file.active
    sheet.title = "predictions"   #definisce un nome per il foglio di lavoro
    sheet["A1"] = "Region"
    sheet["B1"] = "Start RUT"
    sheet["C1"] = "End RUT"
    sheet["D1"] = "Strand"
    
    p = open("Informations about predictions.txt","a")
    p.write("Sequences of predicted Rho-dependent terminators")
    
    #positive strand
    
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
                rapporto = w.count("C")/numG
            else:
                numG = 1
                rapporto = (w.count("C")+1)/numG
            if rapporto > 1:
                if re.search(pattern1,w):
                    dati = []
                    for g in range(50):
                        if x2+g <= len(genome):
                            w = genome[x1+g:x2+g]
                            if w.count("G") > 0:
                                numG = w.count("G")
                                rapporto = w.count("C")/numG
                            else:
                                numG = 1
                                rapporto = (w.count("C")+1)/numG
                            if rapporto > 1:
                                if re.search(pattern1,w):
                                    dati.append(rapporto)
                                else:
                                    dati.append(0)
                            else:
                                dati.append(0)
                        else:
                            break
                    maxP = np.argmax(dati)
                    rapporto = np.max(dati)
                    x1 = x1 + maxP
                    x2 = x2 + maxP
                    scale = x2 - j - 1
                    s = genome[x2:x2+150]
                    w = genome[x1:x2]
                    score = 3
                    ctrl = palindrome_control(s,GC_whole_genome,1)
                    if ctrl > 0 or re.search(pattern_pause_site1,s):
                        sheet["A%d"%(num+1)] = "T%d"%(cod)
                        sheet["B%d"%(num+1)] = x1
                        sheet["C%d"%(num+1)] = x2
                        sheet["D%d"%(num+1)] = "plus"
                        num += 1
                        predictions += 1
                        scale = x2 + 150 - j - 1
                        cg_value.append(rapporto)
                        if rapporto > 2:
                            score += 3
                        elif rapporto > 1.50:
                            score += 2
                        elif rapporto > 1.25:
                            score += 1
                        p.write("\n\n\nPREDICTED REGION NUMBER ")
                        p.write(str(cod))
                        p.write(" (STRAND POSITIVE) ")
                        p.write("\n\nGenomic sequence of putative RUT site (Coordinates:  ")
                        p.write(str(x1))
                        p.write("-")
                        p.write(str(x2))
                        p.write(", c/g = ")
                        p.write(str(rapporto))
                        p.write(" ):      ")
                        p.write(str(w))
                        p.write("\n\nThe 150 nt long genomic region immediately downstream is ")
                        p.write(str(s))
                        p.write("\n")
                        cod += 1
                        final_score = 0
                        final_score = palindrome_finder(s,p,GC_whole_genome,1,score)
                        par = 0
                        while True:
                            if re.search(pattern_pause_site1,s):
                                ricerca = re.search(pattern_pause_site1,s)
                                positions = ricerca.span()
                                x2 = positions[0]
                                y2 = positions[1]
                                s = s[y2:]
                                p.write("\n\nPAUSE-CONSENSUS present at the coordinates ")
                                p.write(str(x2))
                                p.write("-")
                                p.write(str(y2))
                                x2 += par
                                y2 += par
                            else:
                                break
                        if final_score == 0:
                            final_score = score + 3
                        punteggi.append(final_score)    
                        

    #negative strand
    
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
                rapporto = w.count("G")/numC
            else:
                numC = 1
                rapporto = (w.count("G")+1)/numC
            if rapporto > 1:
                if re.search(pattern2,w):
                    dati = []
                    for g in range(50):
                        if x2+g <= len(genome):
                            w = genome[x1+g:x2+g]
                            if w.count("C") > 0:
                                numC = w.count("C")
                                rapporto = w.count("G")/numC
                            else:
                                numC = 1
                                rapporto = (w.count("G")+1)/numC
                            if rapporto > 1:
                                if re.search(pattern2,w):
                                    dati.append(rapporto)
                                else:
                                    dati.append(0)
                            else:
                                dati.append(0)
                        else:
                            break
                    maxP = np.argmax(dati)
                    rapporto = np.max(dati)
                    x1 = x1 + maxP
                    x2 = x2 + maxP
                    scale = x2 - j - 1
                    s = genome[x1-150:x1]
                    w = genome[x1:x2]
                    score = 3
                    ctrl = palindrome_control(s,GC_whole_genome,-1)
                    if ctrl > 0 or re.search(pattern_pause_site2,s):
                        sheet["A%d"%(num+1)] = "T%d"%(cod)
                        sheet["B%d"%(num+1)] = x1
                        sheet["C%d"%(num+1)] = x2
                        sheet["D%d"%(num+1)] = "minus"
                        num += 1
                        predictions += 1
                        scale = x2 + 150 - j - 1
                        cg_value.append(rapporto)
                        if rapporto > 2:
                            score += 3
                        elif rapporto > 1.50:
                            score += 2
                        elif rapporto > 1.25:
                            score += 1
                        p.write("\n\n\nPREDICTED REGION NUMBER ")
                        p.write(str(cod))
                        p.write(" (STRAND NEGATIVE)")
                        p.write("\n\nGenomic sequence of putative RUT site (Coordinates:  ")
                        p.write(str(x1))
                        p.write("-")
                        p.write(str(x2))
                        p.write(", c/g = ")
                        p.write(str(rapporto))
                        p.write(" ):      ")
                        p.write(str(w))
                        p.write("\n\nThe 150 nt long genomic region immediately downstream is ")
                        p.write(str(s))
                        p.write("\n")
                        cod += 1
                        final_score = 0
                        final_score = palindrome_finder(s,p,GC_whole_genome,-1,score)
                        par = 0
                        while True:
                            if re.search(pattern_pause_site2,s):
                                ricerca = re.search(pattern_pause_site2,s)
                                positions = ricerca.span()
                                x2 = positions[0]
                                y2 = positions[1]
                                s = s[y2:]
                                p.write("\n\nPAUSE-CONSENSUS present at the coordinates ")
                                p.write(str(x2))
                                p.write("-")
                                p.write(str(y2))
                                x2 += par
                                y2 += par
                            else:
                                break
                        if final_score == 0:
                            final_score = score + 3
                        punteggi.append(final_score)
    
   
    p.write("\n\n\n\n\n\nTotal number of predicted Rho-dependent terminators: ")
    p.write(str(predictions))
    p.write("\n\n\nMean C/G content of predicted terminators: ")
    p.write(str(np.mean(cg_value)))
    p.write("\nStandard deviation: ")
    p.write(str(np.std(cg_value)))
        
    p.close()
    nuovo_file.save('Rho-dependent terminators coordinates.xlsx')
    
    print("Work finished, see output files in the current directory")
    
except IOError:   
    print ("File %s inexistent in the current directory!" %(File_genome))
    
        
    