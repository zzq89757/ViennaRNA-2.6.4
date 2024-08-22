import re
def replace_str(raw_str:str, data_li:list, line_ct:int) -> str:
    # extract num in raw string and replaced by new data list
    # o_handle.write(raw_str)
    raw_str_num_li = re.findall(r"(\-?\d+\,\s+\-?\d+\,\s+\-?\d+\,\s+\-?\d+\,\s+\-?\d+)",raw_str)
    # o_handle.write(raw_str_num_li)
    # return
    idx = 0
    new_str:str
    new_str = raw_str.replace(raw_str_num_li[0],",   ".join(data_li[line_ct]))
    return new_str
    
    
def process_pair(rna_file:str, dna_file:str, output:str) -> None:
    r_handle = open(rna_file,'r')
    d_handle = open(dna_file,'r')
    o_handle = open(output, 'w')
    # storage data of dna file into list
    all_num_li = []
    for p_line in d_handle:
        # o_handle.write(p_line,end="")
        num_str, comment = p_line.split("/*")
        all_num_li.append(num_str.split()) 
    # replace rna data
    line_ct = 0
    for c_line in r_handle:
        if c_line.find("INF") > -1:
            o_handle.write(c_line)
            continue
        # data need to replace
        if c_line.find(",") > -1:
            replaced_str = replace_str(c_line, all_num_li, line_ct)
            o_handle.write(replaced_str)
            line_ct += 1
        else:
            o_handle.write(c_line)
            ...

def main():
    rna_par_file = "matrices/intl22.h"
    dna_par_file = "/home/wayne/Repository/ViennaRNA-2.6.4/src/bin/fold_pre_process/int22_dna.par"
    output_file = "intl22_dna.h"
    process_pair(rna_par_file, dna_par_file, output_file)
    

if __name__ == "__main__":
    main()