from multitax import NcbiTx


def convert_multiple_files_to_krona_format(input_files, output_file):
    # 初始化 NcbiTx 对象
    tax = NcbiTx()

    # 存储每个 taxid 的计数
    taxid_count = {}

    # 处理每个输入文件
    for input_file in input_files:
        with open(input_file, 'r') as infile:
            for line in infile:
                # 按制表符分割
                parts = line.strip().split('\t')

                # 获取第二列中的 taxid 信息
                taxid_info = parts[1].split('\t')

                for info in taxid_info:
                    if info == "unclassified":
                        taxid_count["unclassified"] = taxid_count.get("unclassified", 0) + 1
                    else:
                        taxid = info.split(':')[0]  # 只取 taxid，忽略后面的数字
                        taxid_count[taxid] = taxid_count.get(taxid, 0) + 1

    # 写入合并后的输出文件
    with open(output_file + '.tsv', 'w') as outfile:
        for taxid, count in taxid_count.items():
            if taxid == "unclassified":
                outfile.write(f"{count + 1}\tunclassified\n")
            else:
                # 使用 NcbiTx 获取物种层级信息
                taxonomy_levels = tax.name_lineage(taxid)
                taxonomy_string = '\t'.join(taxonomy_levels)  # 用 \t 作为分隔符
                outfile.write(f"{count + 1}\t{taxonomy_string}\n")

    print(f"Conversion completed. The combined output file is saved as {output_file}.")
