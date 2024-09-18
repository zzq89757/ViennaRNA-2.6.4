def parse_rnafold_params(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    data_blocks = {}
    current_block = None

    for line in lines:
        line = line.strip()
        if line.startswith("#"):
            # 每遇到一个 #，标识一个新的数据块
            current_block = line[1:].strip().replace(" ", "_")
            data_blocks[current_block] = []
        elif current_block and not line.startswith("/*") and line:
            # 解析数据行并添加到当前数据块中
            line = line.split("/*")[0].strip()
            try:
                row = list(map(int, line.split()))
                data_blocks[current_block].append(row)
            except ValueError:
                continue

    return data_blocks

def array_to_c_format(array, name):
    size = len(array)
    c_code = f"int {name}[{size}][{size}] = {{\n"
    for row in array:
        row_str = ", ".join(map(str, row))
        c_code += f"    {{ {row_str} }},\n"
    c_code = c_code.rstrip(",\n") + "\n};\n"
    return c_code

def main():
    # 输入文件路径
    input_file = '/home/wayne/Repository/ViennaRNA-2.6.4/misc/dna_mathews2004.par'

    # 解析并转换
    data_blocks = parse_rnafold_params(input_file)
    
    # 输出每个数据块的 C 代码
    for block_name, array in data_blocks.items():
        c_code = array_to_c_format(array, block_name)
        print(c_code)

if __name__ == "__main__":
    main()
