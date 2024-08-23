def parse_params_file(params_file):
    handle = open(params_file,'r')
    line_ct = 1
    for line in handle:
        data, comment = line.split("/*")
        comment = "/*" + comment
        data_li = data.split()
        inf_str = ",".join([" INF "] * 5)
        inf_para = r" ,{{{ " + inf_str + " }\n" + (r"  ,{  " + inf_str +"}\n") * 4 + "}\n"
        
        # print(data)
        if line_ct % 35 == 1:
            
            print(inf_para,end="")
        if line_ct % 5 == 1:
            print(r"  ,{{",end="")
        else:
            print(r" ,{ ",end="")
            # print(f"{line}")
        print(", ".join(data_li),end="")
        print(" }",end="   ")
        print(comment,end="")
        if line_ct % 5 == 0:
            print("}")
        if line_ct % 35 == 0:
            
            print(" }")
        line_ct += 1

def main():
    parse_params_file("/home/wayne/Repository/ViennaRNA-2.6.4/src/bin/fold_pre_process/intl21.par")


if __name__ == "__main__":
    main()
            