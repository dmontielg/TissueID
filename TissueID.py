#!/usr/bin/env python


import os
import sys
import time
import argparse
import subprocess

import pandas as pd
import numpy as np
import tensorflow as tf

from argparse import ArgumentParser
from sklearn.preprocessing import MinMaxScaler



# In[2]:


def get_frequency_table(mpileup):
    
    frequency_table = {}
    for i in mpileup.values:       
        fastadict = {"A":0,"T":0,"G":0,"C":0}                    
        sequence = i[4] #actual sequence                
        sequence = sequence.upper()         
        sequence = trimm_caret(sequence)                            
        sequence = sequence.replace("$", "")                   
        indel_pos = find_all_indels(sequence)
        ### Count number of indels
        indels = count_indels(sequence, indel_pos)        
        fastadict.update(indels)
        fastadict["-"] += sequence.count("*")
        ### Trimm Indels
        trimm_sequence = trimm_indels(sequence, indel_pos)        
        for seq in trimm_sequence:    
            if seq in fastadict:            
                fastadict[seq] +=1    
        frequency_table.update({i[1]:list(fastadict.values())})    
    df_frequency_table = pd.DataFrame.from_dict(frequency_table, orient='index')
    df_frequency_table.columns = bases  
    
    return df_frequency_table


def run_bwa_mem(path_fastq_file, output_file, threads, ref_genome):        
    cmd = "bwa mem -v 0 -t {} -B 1 -O 1 -E 1 -L 1 {} {} > {}".format(threads, ref_genome, path_fastq_file, output_file)                    
    try:        
        subprocess.call(cmd, shell=True)
        return True
    except OSError:
        return []


def generate_tmp_folder(tmp_folder):
    if os.path.isdir(tmp_folder):            
        cmd = 'rm -r {}'.format(tmp_folder)
        subprocess.call(cmd, shell=True)
        cmd = 'mkdir {}'.format(tmp_folder)
        subprocess.call(cmd, shell=True)
    else:
        cmd = 'mkdir {}'.format(tmp_folder)
        subprocess.call(cmd, shell=True)



def find_all_indels(s):
    find_all = lambda c,s: [x for x in range(c.find(s), len(c)) if c[x] == s]
    list_pos = []
    for i in find_all(s,"-"):
        list_pos.append(i)
    for i in find_all(s,"+"):
        list_pos.append(i)    
    return sorted(list_pos)



def count_indels(s, pos):    
    dict_indel = {"+":0,"-":0}    
    if pos == []:
        return dict_indel            
    if len(pos) > 0:
        for i in range(0,len(pos)): 
            try: # in case it is not a number but a base pair e.g. A
                dict_indel[s[pos[i]]] += int(s[pos[i]+1])                                                                        
            except ValueError:                
                dict_indel[s[pos[i]]] += 1
                continue                
    return dict_indel



def trimm_indels(s, pos):    
    ## Receives a sequence and trimms indels            
    if pos == []:
        return s
    u_sequence = ""  
    start =  pos[0]
    count = (start+1)    
    try: # in case it is not a number but a base pair e.g. A
        end = count+int(s[count])+1
    except ValueError:                
        end = start+1            
    u_sequence = s[:start]    
    if len(pos) > 1:
        for i in range(1,len(pos)):                      
            start = end                
            u_sequence += s[start:pos[i]]
            start = pos[i]
            count = (start+1)            
            try: # in case it is not a number but a base pair e.g. A
                end = count+int(s[count])+1  
            except ValueError:
                end = start+1                
            if pos[-1] == pos[i]:    
                #print(s[end:])
                u_sequence += s[end:]
    else:        
        u_sequence += s[end:]
        
    return u_sequence


def trimm_caret(s):       
    
    find_all = lambda c,s: [x for x in range(c.find(s), len(c)) if c[x] == s]
    list_pos = []
    for i in find_all(s,"^"):
        list_pos.append(i)    
    if list_pos == []:
        return s
    i = 0    
    start = 0
    end = 0
    sequence = ""
    while i<len(s):
        if s[i] == "^":        
            end = i
            sequence += (s[start:end])                    
            start = i+1
        elif i >= list_pos[-1]+1:
            sequence += (s[list_pos[-1]+1:])
            break
        i+=1
        
    return sequence



def process_frequencies(dirpath):
    
    dirpath = os.walk(dirpath)
    df_all = pd.DataFrame()
    name_files = []
    print("Reading files...")
    for dirpath, dirnames, filenames in dirpath:
        for filename in [f for f in filenames if f.endswith(".freq")]:
            path_freq_file = os.path.join(dirpath, filename)             
            name_freq_file = path_freq_file.split('/')[-1]                  
            name_files.append(name_freq_file)
            df_cpgs_file = pd.read_table(path_freq_file, sep="\t", index_col="position" )                                                            
            df_all = df_all.append(df_cpgs_file.T)                                                         
    df_all = df_all.fillna(0)
    df_all = df_all.T            
    print("Preprocessing Samples:...")
    print(len(name_files))
    X_train = preprocess_df(df_all, name_files)    
    print("Preprocessing Done!")
    
    return X_train, name_files



def transform_df(df_pos):
    df = pd.DataFrame()
    ### transoforms each sample to the sample scale row/all
    for i in range(len(df_pos)):
        row = df_pos.iloc[i]
        row = row/sum(row)
        df = df.append(pd.DataFrame(np.array(row).reshape(1,6), columns = df_pos.columns))
    df = df.fillna(0)
    #scaler = MinMaxScaler()
    #scaler.fit(df)
    #df_transformed = scaler.transform(df)    
    #df_transformed = pd.DataFrame(df_transformed, columns=df_pos.columns, index=df_pos.index)
    
    df_transformed = df
    df_transformed.index = df_pos.index
    return df_transformed



def preprocess_df(df_all, name_files):
        
    df = pd.DataFrame()        
    if df_all.shape[1]/6 > 1:                
        for index in range(len(df_all)):        
            
            df_tmp = pd.DataFrame()
            a = pd.DataFrame(df_all.iloc[index]["A"])
            a.index=name_files
            a.columns = ["A"]
            
            t = pd.DataFrame(df_all.iloc[index]["T"])    
            t.index=name_files
            t.columns = ["T"]    

            g = pd.DataFrame(df_all.iloc[index]["G"])
            g.index=name_files
            g.columns = ["G"]

            c = pd.DataFrame(df_all.iloc[index]["C"])
            c.index=name_files
            c.columns = ["C"]

            insertion = pd.DataFrame(df_all.iloc[index]["+"])
            insertion.index=name_files
            insertion.columns = ["+"]

            deletion = pd.DataFrame(df_all.iloc[index]["-"])
            deletion.index=name_files
            deletion.columns = ["-"]

            df_tmp = pd.concat([df_tmp, a, t, g, c, insertion, deletion], axis=1, sort=False)                
            df_scaled = transform_df(df_tmp)
            
            df = pd.concat([df, df_scaled], axis=1)
    else:
        df_all = transform_df(df_all)          
        for index in range(len(df_all)):
            row = df_all.iloc[index]
            df = pd.concat([df, row], axis=0, sort=False)                                
        df = df.T
        df.index = name_files        
    return df


def predict(X_test, dirpath):
    
    nets = 0
    list_models = []
    y_probs = np.array(0)
    dirpath = os.walk(dirpath)
    load_model = tf.keras.models.load_model        

    for dirpath, dirnames, filenames in dirpath:
        for filename in [f for f in filenames if f.endswith(".model")]:   
            path_model = os.path.join(dirpath, filename)    
            list_models.append(path_model)

    toolbar_width = len(list_models)
    
    # setup toolbar
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['
    
    for i in range(toolbar_width):                                                    
            path_model = os.path.join(dirpath, filename)             
            name_model = path_model.split('/')[-1]                   
            train_model = load_model(path_model)            
            #y_pred = train_model.predict_classes(X_test)
            y_probs = y_probs + train_model.predict(X_test)                                    
            nets += 1               
            time.sleep(0.01) # do real work here            
            sys.stdout.write("=")
            sys.stdout.flush()
    sys.stdout.write("\n")
    
    return y_probs/nets

def folder_exists(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

def dna():
    print("""

`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
   `=`,'=/     `=`,'=/     `=`,'=/     `=`,'=/
     y==/        y==/        y==/        y==/
   ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,-<=`.
,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
        """)

def get_arguments():    

    parser = ArgumentParser(description="\tERASMUS MC \n TissueID: Classification of different forensically \nrelevant human epithelial materials. ")    
    
    parser.add_argument("-fasta", "--FASTA",
        dest="fasta_dir", required=False, type=folder_exists,
        help="input path with fasta files", metavar="DIR")    
    
    parser.add_argument("-fastq", "--FASTQ",
            dest="fastq_dir", required=False, type=folder_exists,                        
            help="Input path with fastq files", metavar="DIR")
        
    parser.add_argument("-out", "--OUTPUT",
            dest="output_dir", required=True,                        
            help="Output folder path", metavar="DIR")
    
    parser.add_argument("-model", "--MODEL",
            dest="model_dir", required=True,                        
            help="Model folder path", metavar="DIR")
    
    parser.add_argument("-pos", "--POS_FILE",
            dest="pos_file", required=True,                        
            help=" Position in bed format", metavar="BED_file")
    
    parser.add_argument("-ref", "--REF",
            dest="ref_genome", required=True,                        
            help="path/Ecoli_K12_ref.fasta", metavar="REF_GENOME")
        
    parser.add_argument("-t", "--Threads", dest="threads",
            help="Set number of additional threads to use [CPUs]",
            type=int,
            default=2)
    
    args = parser.parse_args()    
    return args

def check_extension(filename,ext):
    flag = False
    for x in ext:
        if filename.endswith(x):
            flag = True
    return flag    

def check_if_folder(path,ext):
    
    list_files = []
    if os.path.isdir(path):
        dirpath = os.walk(path)
        for dirpath, dirnames, filenames in dirpath:
            for filename in [f for f in filenames if check_extension(f, ext)]:
                files = os.path.join(dirpath, filename)
                list_files.append(files)
        return list_files
    else:
        return [path]
    

def add_missed_positions(df_freq_table, pos_file):
        
    index_df = np.array(df_freq_table.index)        
    df_pos = pd.read_table(pos_file, sep="\t",header=None)        
    index_pos = df_pos[1].values        
    diff = (list(set(index_pos).difference(set(index_df))))
    
    if len(diff) > 0:
        for d in diff:                                    
            df_freq_table.loc[d] = np.zeros(6)            
        df_freq_table.sort_index(inplace=True)
        return df_freq_table
    else:
        return df_freq_table
    

if __name__ == "__main__":

    print("\tERASMUS MC Department of Genetic Identification \n\n\tTissueID: Classification of different forensically \n\trelevant human epithelial materials. ")    
    dna()
    bases = ["A","T","G","C","+","-"]
    
    args = get_arguments()  
    output_dir = args.output_dir
    model_dir = args.model_dir

    pos_file = args.pos_file
    ref_genome = args.ref_genome
    threads = args.threads    
    ext = []
    i = 1     
    if args.fasta_dir:
        ext.append(".fsa",".fasta",".fa")
        files = check_if_folder(args.fasta_dir,ext) 
    elif args.fastq_dir:
        ext.append(".fastq")
        files = check_if_folder(args.fastq_dir,ext) 
    
    alignment_dir   = "alignments"
    pileup_dir      = "pileups"
    frequency_dir   = "frequencies"

    generate_tmp_folder(output_dir) 
    generate_tmp_folder(output_dir+"/"+alignment_dir) 
    generate_tmp_folder(output_dir+"/"+pileup_dir) 
    generate_tmp_folder(output_dir+"/"+frequency_dir)     
     
    for path_fastq_file in files:            
        name_fastq_file = path_fastq_file.split('/')[-1]                                                
            
        sam_file = output_dir+"/"+alignment_dir+"/"+name_fastq_file+".sam"                    
        bam_file = output_dir+"/"+alignment_dir+"/"+name_fastq_file+".bam"                    
        mpileup = output_dir+"/"+pileup_dir+"/"+name_fastq_file+".mpileup"        
        frequency_file = output_dir+"/"+frequency_dir+"/"+name_fastq_file+".freq" 
            
        print(str(i)+".-Preprocessing file: "+name_fastq_file+" ...")                
        print("\tAligning with {}".format(ref_genome))        
        run_bwa_mem(path_fastq_file, sam_file, threads, ref_genome)                            
        cmd = "samtools view -@ {} -bS {} | samtools sort -@ {} -m 2G -o {}".format(threads, sam_file, threads, bam_file)                                
        subprocess.call(cmd, shell=True)                                
        cmd = "samtools index -@ {} {}".format(threads, bam_file)        
        subprocess.call(cmd, shell=True)                    
        print("\tGenerating pileup")        
        cmd = "samtools mpileup -d 100000 -l {} {} > {}".format(pos_file, bam_file, mpileup )                        
        subprocess.call(cmd, shell=True)        
        mpileup = pd.read_table(mpileup, names=["chr","pos","ref","reads","seq","qual"])
        print("\tGenerating frequency table ")
        df_freq_table = get_frequency_table(mpileup)                                                
        df_freq_table = add_missed_positions(df_freq_table, pos_file)                                
        df_freq_table.index.names = ['position']
        df_freq_table.to_csv(frequency_file, sep="\t", index=True)                
        i = i+1            
    X_test, name_files_test = process_frequencies(output_dir+"/"+frequency_dir)
    
    y_probs = predict(X_test.values, model_dir)
    df_probs = pd.DataFrame(y_probs,name_files_test)
    df_probs.columns = ["Skin", "Oral","Vagina"]
    df_probs.to_csv(output_dir+"/"+"predictions.csv")
    print(df_probs)
    print("Done!")    
    
