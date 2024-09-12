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

	static constexpr int first_appearance_init = -1;    ///< ��ʼֵ����ʾδ���ֹ���

	robin_hood::unordered_map< std::string, std::vector< std::string > > m_parents; ///< �洢���ڵ㵽�ӽڵ��ӳ�䡣
	std::vector< int > m_euler;                                      ///< �洢ŷ�������Ľڵ��š�
	std::vector< int > m_depth;                                      ///< �洢ÿ���ڵ����ȡ�
	std::vector< int > m_firstAppearance;                            ///< �洢ÿ���ڵ��һ�γ��ֵ�λ�á�
	int                m_vertices = 0;                               ///< �ڵ�������
	robin_hood::unordered_map< std::string, unsigned int > m_encode; ///< �ڵ����Ƶ���ŵ�ӳ�䡣
	std::vector< std::string >                             m_decode; ///< �ڵ��ŵ����Ƶ�ӳ�䡣
	std::vector< std::vector< int > >                      m_M;      ///< RMQ Ԥ�������
};

inline void LCA::addEdge(const std::string& father, const std::string& son)
{
	// ������ڵ㻹û�б�����
	if (m_encode.count(father) == 0)
	{
		// �����ڵ�������ӳ�䣬�������µĽڵ���
		m_encode.insert({ father, m_vertices });
		// �ڽ��������в��븸�ڵ����ƣ�λ��Ϊ�±�Ŷ�Ӧ������
		m_decode.insert(m_decode.begin() + m_vertices, father);
		// ���ӽڵ�����
		++m_vertices;
	}

	// ����ӽڵ㻹û�б�����
	if (m_encode.count(son) == 0)
	{
		// ���ӽڵ�������ӳ�䣬�������µĽڵ���
		m_encode.insert({ son, m_vertices });
		// �ڽ��������в����ӽڵ����ƣ�λ��Ϊ�±�Ŷ�Ӧ������
		m_decode.insert(m_decode.begin() + m_vertices, son);
		// ���ӽڵ�����
		++m_vertices;
	}

	// ������ڵ㻹û���κ��ӽڵ㣬��ʼ�����ӽڵ��б�
	if (m_parents.count(father) == 0)
	{
		m_parents[father] = { son };
	}
	else
	{
		// ����ֱ�ӽ��ӽڵ���ӵ����ڵ���ӽڵ��б���
		m_parents[father].emplace_back(son);
	}
}

inline void LCA::depthFirstSearch(const std::string& current, unsigned int depth)
{
	// ��ȡ��ǰ�ڵ�ı���ֵ
	const auto currentEncoded = m_encode[current];

	// ��¼��ǰ�ڵ��һ�γ��ֵ�λ��
	if (m_firstAppearance[currentEncoded] == first_appearance_init)
	{
		m_firstAppearance[currentEncoded] = m_euler.size();
	}

	// ����ǰ�ڵ���ӵ�ŷ��������
	m_euler.push_back(currentEncoded);
	// ����ǰ�ڵ�������ӵ����������
	m_depth.push_back(depth);

	// ������ǰ�ڵ�������ӽڵ�
	for (const auto& node : m_parents[current])
	{
		// �ݹ�ض�ÿ���ӽڵ�������������������ȼ�1
		depthFirstSearch(node, depth + 1);
		// �ٴν���ǰ�ڵ���ӵ�ŷ��������
		m_euler.push_back(currentEncoded);
		// �ٴν���ǰ�ڵ�������ӵ����������
		m_depth.push_back(depth);
	}
}

inline void LCA::doEulerWalk(const std::string& root_node)
{
	// ��ʼ�� m_firstAppearance ���飬��СΪ�ڵ���������ʼֵΪ first_appearance_init
	m_firstAppearance.resize(m_vertices, first_appearance_init);
	// �Ӹ��ڵ㿪ʼ�������������������ʼ���Ϊ0
	depthFirstSearch(root_node, 0);
	// Ԥ����������Сֵ��ѯ
	preProcessRMQ();
}

// <O(N logN) Ԥ����ʱ�临�Ӷ�, O(1)��ѯʱ�临�Ӷ�>
inline void LCA::preProcessRMQ()
{
	// ��ȡ�������Ĵ�С
	const auto size = m_depth.size();
	// �����������Ķ���������ϡ���
	const int logDepth = std::ceil(std::log2(m_depth.size()));

	// ��ʼ�� m_M ϡ�����СΪ [size x logDepth]
	m_M.resize(size, std::vector< int >(logDepth));

	// ��ʼ��ϡ���ĵ�һ��
	for (auto i = 0u; i < size; ++i)
	{
		m_M[i].front() = i;
	}

	// ���մ�С���䵽�������˳�����ϡ����ֵ
	for (unsigned int j = 1; 1u << j <= size; j++)
	{
		for (unsigned int i = 0; i + (1 << j) - 1 < size; i++)
		{
			// �Ƚ���������������������Сֵ��ѡ���Сֵ��Ӧ������
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
	// ��� i ���� j���򽻻� i �� j ��ֵ����ȷ�� i <= j
	if (i > j)
	{
		std::swap(i, j);
	}

	// �������� [i, j] �ĳ��ȵĶ���ֵ k
	const auto k = static_cast<int>(std::log2(j - i + 1));
	// ��ȡ���� [i, i + 2^k - 1] ����Сֵ����
	const auto term1 = m_M[i][k];
	// ��ȡ���� [j - 2^k + 1, j] ����Сֵ����
	const auto term2 = m_M[j - (1 << k) + 1][k];

	// �������������н�Сֵ������
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
	// ȷ�� u �� v �������ںϷ���Χ��
	assert(u < m_firstAppearance.size());
	assert(v < m_firstAppearance.size());

	// ���������� u �� v ��ͬ��ֱ�ӷ��� u
	if (u == v)
	{
		return u;
	}

	// ȷ�� u ���״γ���λ��С�� v ���״γ���λ�ã���������򽻻� u �� v
	if (m_firstAppearance[u] > m_firstAppearance[v])
	{
		std::swap(u, v);
	}

	// �� [m_firstAppearance[u], m_firstAppearance[v]] ��Χ��ִ�� RMQ ��ѯ������ŷ�������е�����
	return m_euler[queryRMQ(m_firstAppearance[u], m_firstAppearance[v])];
}

inline std::string LCA::getLCA(const std::vector< std::string >& taxIds) const
{
	// ȷ�� taxIds ������������Ԫ��
	assert(taxIds.size() > 1);

	// ��ʼ�� lca Ϊ��һ���͵ڶ����ڵ�� LCA
	int lca = getLCA(m_encode.at(taxIds[0]), m_encode.at(taxIds[1]));
	// ���μ��� lca ��ÿ�������ڵ�� LCA
	for (unsigned int i = 2; i < taxIds.size(); ++i)
		lca = getLCA(lca, m_encode.at(taxIds[i]));

	// �������� LCA ��Ӧ�Ľ���ֵ
	return m_decode.at(lca);
}