from multitax import NcbiTx
import math


def calculate_shannon_index(taxon_dict):
    total_count = sum(taxon_dict.values())
    if total_count == 0:
        return 0
    shannon_index = 0
    for count in taxon_dict.values():
        p_i = count / total_count
        if p_i > 0:
            shannon_index -= p_i * math.log(p_i)
    return shannon_index


def calculate_simpson_index(taxon_dict):
    total_count = sum(taxon_dict.values())
    if total_count == 0:
        return 0
    simpson_index = 0
    for count in taxon_dict.values():
        p_i = count / total_count
        simpson_index += p_i ** 2
    return 1 - simpson_index


def process_file(input_files, output_file):
    # 初始化 NcbiTx 对象
    tax = NcbiTx()

    # 存储不同分类层级的计数，包括clade和unclassified
    count_by_level = {
        "superkingdom": {},  # 域
        "clade": {},  # 进化枝
        "phylum": {},  # 门
        "class": {},  # 纲
        "order": {},  # 目
        "family": {},  # 科
        "genus": {},  # 属
        "species": {},  # 种
        "strain": {},  # 菌株
    }

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
                        # 如果是未分类的，所有层级都加1
                        for level in count_by_level:
                            count_by_level[level]["unclassified"] = count_by_level[level].get("unclassified", 0) + 1
                    else:
                        taxid = info.split(':')[0]  # 只取 taxid，忽略后面的数字

                        # 使用 NcbiTx 获取物种的层级信息
                        try:
                            taxonomy_levels = tax.name_lineage(taxid)  # 获取名称层级
                            rank_levels = tax.rank_lineage(taxid)  # 获取层级 rank 信息

                            # 遍历 rank_levels，匹配到具体的分类等级
                            has_level = {level: False for level in count_by_level}  # 跟踪每个层级是否有数据

                            for name, rank in zip(taxonomy_levels, rank_levels):
                                if rank == "no rank":
                                    continue  # 跳过没有分类层级的条目
                                if rank == "clade":
                                    count_by_level["clade"][name] = count_by_level["clade"].get(name, 0) + 1
                                    has_level["clade"] = True
                                elif rank == "strain":
                                    count_by_level["strain"][name] = count_by_level["strain"].get(name, 0) + 1
                                    has_level["strain"] = True
                                elif rank in count_by_level:
                                    count_by_level[rank][name] = count_by_level[rank].get(name, 0) + 1
                                    has_level[rank] = True

                            # 对于没有分类信息的层级，记录为“unclassified”
                            for level in has_level:
                                if not has_level[level]:
                                    count_by_level[level]["unclassified"] = count_by_level[level].get("unclassified",
                                                                                                      0) + 1

                        except:
                            # 查询失败则归类为未分类
                            for level in count_by_level:
                                count_by_level[level]["unclassified"] = count_by_level[level].get("unclassified", 0) + 1

    # 计算多样性指标及丰度，并输出
    total_counts_by_level = {level: sum(taxon_dict.values()) for level, taxon_dict in count_by_level.items()}

    with open(output_file + '.tsv', 'w') as outfile:
        # 写入头部信息
        outfile.write("Taxon\tCount\tRelative Abundance (%)\tShannon Index\tSimpson Index\n")

        for level, taxon_dict in count_by_level.items():
            shannon_index = calculate_shannon_index(taxon_dict)
            simpson_index = calculate_simpson_index(taxon_dict)
            total_count = total_counts_by_level[level]

            # 输出层级标题
            outfile.write(f"\n## {level.capitalize()} Level ##\n")

            # 按照计数从大到小排序，并将unclassified放在最后
            sorted_items = []
            unclassified_item = None
            
            # 分离unclassified和其他项
            for taxon, count in taxon_dict.items():
                if taxon == "unclassified":
                    unclassified_item = (taxon, count)
                else:
                    sorted_items.append((taxon, count))
            
            # 对其他项按计数从大到小排序
            sorted_items.sort(key=lambda x: x[1], reverse=True)
            
            # 如果有unclassified项，添加到最后
            if unclassified_item:
                sorted_items.append(unclassified_item)
            
            # 输出排序后的结果
            for taxon, count in sorted_items:
                relative_abundance = (count / total_count) * 100 if total_count > 0 else 0
                # 输出每个分类单元的计数和相对丰度
                outfile.write(
                    f"{taxon}\t{count}\t{relative_abundance:.2f}\t{shannon_index:.4f}\t{simpson_index:.4f}\n")
