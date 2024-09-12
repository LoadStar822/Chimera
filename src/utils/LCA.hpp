#pragma once

#include <robin_hood.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>

class LCA
{
public:
	LCA() = default;

	void        addEdge(const std::string& father, const std::string& son);
	void        doEulerWalk(const std::string& root_node);
	int         getLCA(unsigned int u, unsigned int v) const;
	std::string getLCA(const std::vector< std::string >& taxIds) const;

private:
	void depthFirstSearch(const std::string& current, unsigned int depth);
	void preProcessRMQ();
	int  queryRMQ(unsigned int i, unsigned int j) const;

	static constexpr int first_appearance_init = -1;    ///< 初始值，表示未出现过。

	robin_hood::unordered_map< std::string, std::vector< std::string > > m_parents; ///< 存储父节点到子节点的映射。
	std::vector< int > m_euler;                                      ///< 存储欧拉遍历的节点编号。
	std::vector< int > m_depth;                                      ///< 存储每个节点的深度。
	std::vector< int > m_firstAppearance;                            ///< 存储每个节点第一次出现的位置。
	int                m_vertices = 0;                               ///< 节点总数。
	robin_hood::unordered_map< std::string, unsigned int > m_encode; ///< 节点名称到编号的映射。
	std::vector< std::string >                             m_decode; ///< 节点编号到名称的映射。
	std::vector< std::vector< int > >                      m_M;      ///< RMQ 预处理矩阵。
};

inline void LCA::addEdge(const std::string& father, const std::string& son)
{
	// 如果父节点还没有被编码
	if (m_encode.count(father) == 0)
	{
		// 将父节点插入编码映射，并分配新的节点编号
		m_encode.insert({ father, m_vertices });
		// 在解码向量中插入父节点名称，位置为新编号对应的索引
		m_decode.insert(m_decode.begin() + m_vertices, father);
		// 增加节点总数
		++m_vertices;
	}

	// 如果子节点还没有被编码
	if (m_encode.count(son) == 0)
	{
		// 将子节点插入编码映射，并分配新的节点编号
		m_encode.insert({ son, m_vertices });
		// 在解码向量中插入子节点名称，位置为新编号对应的索引
		m_decode.insert(m_decode.begin() + m_vertices, son);
		// 增加节点总数
		++m_vertices;
	}

	// 如果父节点还没有任何子节点，初始化其子节点列表
	if (m_parents.count(father) == 0)
	{
		m_parents[father] = { son };
	}
	else
	{
		// 否则，直接将子节点添加到父节点的子节点列表中
		m_parents[father].emplace_back(son);
	}
}

inline void LCA::depthFirstSearch(const std::string& current, unsigned int depth)
{
	// 获取当前节点的编码值
	const auto currentEncoded = m_encode[current];

	// 记录当前节点第一次出现的位置
	if (m_firstAppearance[currentEncoded] == first_appearance_init)
	{
		m_firstAppearance[currentEncoded] = m_euler.size();
	}

	// 将当前节点添加到欧拉序列中
	m_euler.push_back(currentEncoded);
	// 将当前节点的深度添加到深度序列中
	m_depth.push_back(depth);

	// 遍历当前节点的所有子节点
	for (const auto& node : m_parents[current])
	{
		// 递归地对每个子节点进行深度优先搜索，深度加1
		depthFirstSearch(node, depth + 1);
		// 再次将当前节点添加到欧拉序列中
		m_euler.push_back(currentEncoded);
		// 再次将当前节点的深度添加到深度序列中
		m_depth.push_back(depth);
	}
}

inline void LCA::doEulerWalk(const std::string& root_node)
{
	// 初始化 m_firstAppearance 数组，大小为节点总数，初始值为 first_appearance_init
	m_firstAppearance.resize(m_vertices, first_appearance_init);
	// 从根节点开始进行深度优先搜索，初始深度为0
	depthFirstSearch(root_node, 0);
	// 预处理区间最小值查询
	preProcessRMQ();
}

// <O(N logN) 预处理时间复杂度, O(1)查询时间复杂度>
inline void LCA::preProcessRMQ()
{
	// 获取深度数组的大小
	const auto size = m_depth.size();
	// 计算深度数组的对数（用于稀疏表）
	const int logDepth = std::ceil(std::log2(m_depth.size()));

	// 初始化 m_M 稀疏表，大小为 [size x logDepth]
	m_M.resize(size, std::vector< int >(logDepth));

	// 初始化稀疏表的第一列
	for (auto i = 0u; i < size; ++i)
	{
		m_M[i].front() = i;
	}

	// 按照从小区间到大区间的顺序计算稀疏表的值
	for (unsigned int j = 1; 1u << j <= size; j++)
	{
		for (unsigned int i = 0; i + (1 << j) - 1 < size; i++)
		{
			// 比较深度数组中两个区间的最小值，选择较小值对应的索引
			if (m_depth[m_M[i][j - 1]] < m_depth[m_M[i + (1 << (j - 1))][j - 1]])
			{
				m_M[i][j] = m_M[i][j - 1];
			}
			else
			{
				m_M[i][j] = m_M[i + (1 << (j - 1))][j - 1];
			}
		}
	}
}

inline int LCA::queryRMQ(unsigned int i, unsigned int j) const
{
	// 如果 i 大于 j，则交换 i 和 j 的值，以确保 i <= j
	if (i > j)
	{
		std::swap(i, j);
	}

	// 计算区间 [i, j] 的长度的对数值 k
	const auto k = static_cast<int>(std::log2(j - i + 1));
	// 获取区间 [i, i + 2^k - 1] 的最小值索引
	const auto term1 = m_M[i][k];
	// 获取区间 [j - 2^k + 1, j] 的最小值索引
	const auto term2 = m_M[j - (1 << k) + 1][k];

	// 返回两个区间中较小值的索引
	if (m_depth[term1] <= m_depth[term2])
	{
		return term1;
	}
	else
	{
		return term2;
	}
}

inline int LCA::getLCA(unsigned int u, unsigned int v) const
{
	// 确保 u 和 v 的索引在合法范围内
	assert(u < m_firstAppearance.size());
	assert(v < m_firstAppearance.size());

	// 简单情况，如果 u 和 v 相同，直接返回 u
	if (u == v)
	{
		return u;
	}

	// 确保 u 的首次出现位置小于 v 的首次出现位置，如果不是则交换 u 和 v
	if (m_firstAppearance[u] > m_firstAppearance[v])
	{
		std::swap(u, v);
	}

	// 在 [m_firstAppearance[u], m_firstAppearance[v]] 范围内执行 RMQ 查询，返回欧拉序列中的索引
	return m_euler[queryRMQ(m_firstAppearance[u], m_firstAppearance[v])];
}

inline std::string LCA::getLCA(const std::vector< std::string >& taxIds) const
{
	// 确保 taxIds 中至少有两个元素
	assert(taxIds.size() > 1);

	// 初始化 lca 为第一个和第二个节点的 LCA
	int lca = getLCA(m_encode.at(taxIds[0]), m_encode.at(taxIds[1]));
	// 依次计算 lca 与每个后续节点的 LCA
	for (unsigned int i = 2; i < taxIds.size(); ++i)
		lca = getLCA(lca, m_encode.at(taxIds[i]));

	// 返回最终 LCA 对应的解码值
	return m_decode.at(lca);
}