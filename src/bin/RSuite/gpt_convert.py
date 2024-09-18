import numpy as np

# 读取数据文件，假设文件中每行数据是以若干个空格分隔的
with open('/home/wayne/Repository/ViennaRNA-2.6.4/src/bin/fold_pre_process/intl21.par', 'r') as file:
    lines = file.readlines()

# 将每行数据拆分并转换为整数，然后存储在一个列表中
data = []
for line in lines:
    num_str = line.split("/*")[0]  # 去掉注释部分
    if num_str.strip():  # 确保不为空行
        data.append([int(x) for x in num_str.split()])

# 将列表转换为 numpy 数组
data = np.array(data)

# 验证数据形状是否正确（应为 1225 行 5 列）
if data.shape != (1225, 5):
    raise ValueError(f"数据形状错误，当前形状为 {data.shape}，预期形状为 (1225, 5)")

# 定义一个全为 INF 的 5*5*5 数组
inf_5_5_5 = np.full((5, 5, 5), np.inf)

# 定义一个全为 INF 的 8*5*5*5 数组
inf_8_5_5_5 = np.full((8, 5, 5, 5), np.inf)

# 初始化最终结果列表
final_data = []

# 填充数据
data_index = 0
for i in range(8):
    if i % 7 == 0 and i != 0:
        final_data.append(inf_8_5_5_5)  # 在每 7 个 8*5*5*5 前添加全 INF 数组
    inner_data = []
    for j in range(8):
        if j % 7 == 0 and j != 0:
            inner_data.append(inf_5_5_5)  # 在每 7 个 5*5*5 前添加全 INF 数组
        inner_data.append(data[data_index:data_index + 7].reshape(7, 5, 5, 5))
        data_index += 7
    final_data.append(np.array(inner_data))

# 将最终结果转换为 numpy 数组
final_data = np.array(final_data)

# 输出最终的 C 语言数组格式
c_array = str(final_data.tolist()).replace('[', '{').replace(']', '}')
print(c_array)
