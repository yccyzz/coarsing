#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <regex>
#include <queue>
#include <cmath>
#include <filesystem>
#include <functional>

using namespace std;

// 实例（点）结构
struct Instance {
    string name;           // 实例名 (如 inst_19132)
    string type;           // 类型 (如 FDRE, LUT6)
    vector<string> nets;   // 连接的所有net
    double x, y;           // 坐标
    int die;               // 所属die (0,1,2,3)

    Instance() : x(0), y(0), die(-1) {}
};

// 图粗化类
class NetlistGraphCoarsening {
private:
    unordered_map<string, Instance> instances;  // inst_name -> Instance
    unordered_map<string, vector<string>> net_to_insts;  // net_name -> [inst_names]
    unordered_map<string, int> inst_to_idx;  // inst_name -> index
    vector<string> idx_to_inst;  // index -> inst_name

    // 邻接表（基于共享net的连接）
    unordered_map<int, unordered_set<int>> adj_list;

    // 边界节点信息
    unordered_set<string> boundary_node_names;
    unordered_map<string, int> node_boundary_id;  // 节点名 -> 边界ID (0,1,2)
    unordered_map<int, int> idx_boundary_id;      // 节点索引 -> 边界ID

    // die边界位置
    vector<double> die_boundaries = {120.0, 240.0, 360.0};  // die0|1, die1|2, die2|3的边界

public:
    string filename;  // 文件名

    NetlistGraphCoarsening() {}

    // 从坐标文件加载坐标 (.pl文件)
    bool readPlaceFile(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Cannot open placement file: " << filename << endl;
            return false;
        }

        string line;
        int loaded = 0;

        cout << "\nLoading placement file: " << filename << endl;

        while (getline(file, line)) {
            istringstream iss(line);
            string inst_name;
            double x, y, z;
            string fixed_flag;

            // 格式: inst_name x y z [FIXED]
            // 第三列是y坐标，用于判断die
            if (iss >> inst_name >> x >> y >> z) {
                // 可能有FIXED标记，忽略它
                iss >> fixed_flag;

                // 如果实例不存在，先创建它
                if (instances.find(inst_name) == instances.end()) {
                    Instance inst;
                    inst.name = inst_name;
                    instances[inst_name] = inst;
                }

                instances[inst_name].x = x;
                instances[inst_name].y = y;

                // 根据y坐标判断die
                if (y < 120) {
                    instances[inst_name].die = 0;
                } else if (y < 240) {
                    instances[inst_name].die = 1;
                } else if (y < 360) {
                    instances[inst_name].die = 2;
                } else {
                    instances[inst_name].die = 3;
                }

                loaded++;
            }
        }

        file.close();

        cout << "Placement loading completed! Loaded " << loaded << " instances" << endl;

        // 统计die分布
        vector<int> die_counts(4, 0);
        for (const auto& [name, inst] : instances) {
            if (inst.die >= 0 && inst.die < 4) {
                die_counts[inst.die]++;
            }
        }

        cout << "Die distribution:" << endl;
        for (int i = 0; i < 4; i++) {
            cout << "  Die " << i << ": " << die_counts[i] << " instances" << endl;
        }

        return true;
    }

    // 解析网表文件 (.osv或.v文件)
    bool readNetFile(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Cannot open netlist file: " << filename << endl;
            return false;
        }

        string line;
        string current_inst_name;
        string current_inst_type;
        vector<string> current_nets;

        // 正则表达式匹配
        regex inst_regex(R"(^(\w+)\s+(\w+)\s*\()");  // 匹配 "FDRE inst_19132 ("
        regex net_regex(R"(\.(\w+)\(([\w_]+)\))");    // 匹配 ".Q(net_19079)"

        cout << "\nParsing netlist file: " << filename << endl;

        while (getline(file, line)) {
            // 去除前后空格
            line.erase(0, line.find_first_not_of(" \t"));
            line.erase(line.find_last_not_of(" \t") + 1);

            if (line.empty() || line[0] == '#') continue;

            // 匹配实例声明行
            smatch match;
            if (regex_search(line, match, inst_regex)) {
                // 如果之前有实例，先保存
                if (!current_inst_name.empty()) {
                    saveInstance(current_inst_name, current_inst_type, current_nets);
                }

                // 开始新实例
                current_inst_type = match[1];
                current_inst_name = match[2];
                current_nets.clear();
            }

            // 提取net连接
            string::const_iterator search_start(line.cbegin());
            while (regex_search(search_start, line.cend(), match, net_regex)) {
                string net_name = match[2];
                current_nets.push_back(net_name);
                search_start = match.suffix().first;
            }

            // 如果遇到 ); 表示实例结束
            if (line.find(");") != string::npos && !current_inst_name.empty()) {
                saveInstance(current_inst_name, current_inst_type, current_nets);
                current_inst_name.clear();
                current_nets.clear();
            }
        }

        file.close();

        cout << "Parsing completed!" << endl;
        cout << "  Total instances: " << instances.size() << endl;
        cout << "  Total nets: " << net_to_insts.size() << endl;

        return true;
    }

    // 保存实例信息
    void saveInstance(const string& name, const string& type, const vector<string>& nets) {
        // 如果实例已存在（从坐标文件读取），更新类型和nets
        if (instances.find(name) != instances.end()) {
            instances[name].type = type;
            instances[name].nets = nets;
        } else {
            // 否则创建新实例
            Instance inst;
            inst.name = name;
            inst.type = type;
            inst.nets = nets;
            instances[name] = inst;
        }

        // 建立net到instance的反向映射
        for (const auto& net : nets) {
            net_to_insts[net].push_back(name);
        }
    }

    // 构建图（两个实例共享net则相连）
    void buildGraph(int max_net_size = 100) {
        cout << "\nBuilding graph..." << endl;
        cout << "  Max net size for full connectivity: " << max_net_size << endl;

        // 建立索引映射
        int idx = 0;
        for (const auto& [name, inst] : instances) {
            inst_to_idx[name] = idx;
            idx_to_inst.push_back(name);
            idx++;
        }

        // 基于共享net构建边
        int edge_count = 0;
        int skipped_nets = 0;
        int large_net_edges = 0;

        for (const auto& [net, inst_list] : net_to_insts) {
            size_t net_size = inst_list.size();

            // 对于小型net：全连接
            if (net_size <= max_net_size) {
                for (size_t i = 0; i < inst_list.size(); i++) {
                    for (size_t j = i + 1; j < inst_list.size(); j++) {
                        int idx_i = inst_to_idx[inst_list[i]];
                        int idx_j = inst_to_idx[inst_list[j]];

                        adj_list[idx_i].insert(idx_j);
                        adj_list[idx_j].insert(idx_i);
                        edge_count++;
                    }
                }
            }
                // 对于大型net：星形连接（连接到第一个节点）
            else {
                // 将所有节点连接到第一个节点，形成星形
                int idx_center = inst_to_idx[inst_list[0]];

                for (size_t i = 1; i < inst_list.size(); i++) {
                    int idx_i = inst_to_idx[inst_list[i]];

                    adj_list[idx_center].insert(idx_i);
                    adj_list[idx_i].insert(idx_center);
                    large_net_edges++;
                }

                skipped_nets++;
            }

            // 显示进度
            if ((edge_count + large_net_edges) % 100000 == 0) {
                cout << "  Progress: " << edge_count + large_net_edges << " edges processed\r" << flush;
            }
        }

        cout << endl;
        cout << "  Regular edges: " << edge_count << endl;
        cout << "  Large net edges (star topology): " << large_net_edges << endl;
        cout << "  Large nets (>" << max_net_size << " instances): " << skipped_nets << endl;
        cout << "  Total edges: " << edge_count + large_net_edges << endl;
    }

    // 计算点到最近die边界的距离
    double distanceToDieBoundary(double y) {
        double min_dist = 1e9;

        // 检查到每条die边界的距离
        for (double boundary : die_boundaries) {
            double dist = abs(y - boundary);
            min_dist = min(min_dist, dist);
        }

        return min_dist;
    }

    // 先筛选边界节点（在建图之前）
    vector<int> filterBoundaryNodes(double threshold = 15.0) {
        vector<int> boundary_nodes;
        unordered_set<string> boundary_node_names;
        unordered_map<string, int> node_boundary_id;  // 节点 -> 边界ID (0,1,2)

        // 第一遍：找出所有边界节点并标记所属边界
        for (const auto& [name, inst] : instances) {
            double y = inst.y;

            // 检查属于哪个边界
            int boundary_id = -1;
            if (abs(y - die_boundaries[0]) <= threshold) {
                boundary_id = 0;  // die0|1边界
            } else if (abs(y - die_boundaries[1]) <= threshold) {
                boundary_id = 1;  // die1|2边界
            } else if (abs(y - die_boundaries[2]) <= threshold) {
                boundary_id = 2;  // die2|3边界
            }

            if (boundary_id >= 0) {
                boundary_node_names.insert(name);
                node_boundary_id[name] = boundary_id;
            }
        }

        cout << "\nFiltering die boundary nodes (distance <= " << threshold << "):" << endl;
        cout << "  Boundary nodes: " << boundary_node_names.size()
             << " / " << instances.size() << endl;

        // 统计每个die边界附近的节点数
        vector<int> boundary_counts(3, 0);
        for (const auto& [name, bid] : node_boundary_id) {
            boundary_counts[bid]++;
        }

        cout << "  Boundary distribution:" << endl;
        cout << "    die0|1 boundary (y=" << die_boundaries[0] << "): " << boundary_counts[0] << " nodes" << endl;
        cout << "    die1|2 boundary (y=" << die_boundaries[1] << "): " << boundary_counts[1] << " nodes" << endl;
        cout << "    die2|3 boundary (y=" << die_boundaries[2] << "): " << boundary_counts[2] << " nodes" << endl;

        // 保存边界节点名称集合和边界ID映射
        this->boundary_node_names = boundary_node_names;
        this->node_boundary_id = node_boundary_id;

        return boundary_nodes;  // 这里返回空的，实际索引在buildGraph后创建
    }

    // 只为边界节点建图
    void buildGraphForBoundaryNodes(int max_net_size = 100) {
        cout << "\nBuilding graph for boundary nodes only..." << endl;
        cout << "  Max net size for full connectivity: " << max_net_size << endl;

        // 只为边界节点建立索引映射，同时记录边界ID
        int idx = 0;
        for (const auto& name : boundary_node_names) {
            inst_to_idx[name] = idx;
            idx_to_inst.push_back(name);
            idx_boundary_id[idx] = node_boundary_id[name];  // 记录边界ID
            idx++;
        }

        cout << "  Indexed " << idx_to_inst.size() << " boundary nodes" << endl;

        // 基于共享net构建边（只考虑边界节点之间的连接）
        int edge_count = 0;
        int skipped_nets = 0;
        int large_net_edges = 0;
        int nets_with_boundary = 0;  // 连接边界节点的net数量
        int nets_only_one_boundary = 0;  // 只连接1个边界节点的net

        for (const auto& [net, inst_list] : net_to_insts) {
            // 只保留边界节点
            vector<string> boundary_insts;
            for (const auto& inst_name : inst_list) {
                if (boundary_node_names.count(inst_name)) {
                    boundary_insts.push_back(inst_name);
                }
            }

            // 统计：有多少net连接边界节点
            if (boundary_insts.size() >= 1) {
                nets_with_boundary++;
            }

            // 如果这个net只连接1个边界节点，无法建边
            if (boundary_insts.size() == 1) {
                nets_only_one_boundary++;
            }

            // 如果这个net没有连接2个以上边界节点，跳过
            if (boundary_insts.size() < 2) continue;

            size_t net_size = boundary_insts.size();

            // 对于小型net：全连接
            if (net_size <= max_net_size) {
                for (size_t i = 0; i < boundary_insts.size(); i++) {
                    for (size_t j = i + 1; j < boundary_insts.size(); j++) {
                        int idx_i = inst_to_idx[boundary_insts[i]];
                        int idx_j = inst_to_idx[boundary_insts[j]];

                        adj_list[idx_i].insert(idx_j);
                        adj_list[idx_j].insert(idx_i);
                        edge_count++;
                    }
                }
            }
                // 对于大型net：星形连接
            else {
                int idx_center = inst_to_idx[boundary_insts[0]];

                for (size_t i = 1; i < boundary_insts.size(); i++) {
                    int idx_i = inst_to_idx[boundary_insts[i]];

                    adj_list[idx_center].insert(idx_i);
                    adj_list[idx_i].insert(idx_center);
                    large_net_edges++;
                }

                skipped_nets++;
            }

            // 显示进度
            if ((edge_count + large_net_edges) % 10000 == 0 && (edge_count + large_net_edges) > 0) {
                cout << "  Progress: " << edge_count + large_net_edges << " edges\r" << flush;
            }
        }

        cout << endl;
        cout << "  Graph statistics:" << endl;
        cout << "    Total nets: " << net_to_insts.size() << endl;
        cout << "    Nets connecting boundary nodes: " << nets_with_boundary << endl;
        cout << "    Nets with only 1 boundary node: " << nets_only_one_boundary << " (cannot form edges)" << endl;
        cout << "    Regular edges: " << edge_count << endl;
        cout << "    Large net edges (star topology): " << large_net_edges << endl;
        cout << "    Large nets (>" << max_net_size << " boundary nodes): " << skipped_nets << endl;
        cout << "    Total edges: " << edge_count + large_net_edges << endl;

        // 统计连通性
        int isolated_nodes = 0;
        for (int i = 0; i < idx_to_inst.size(); i++) {
            if (adj_list.find(i) == adj_list.end() || adj_list[i].empty()) {
                isolated_nodes++;
            }
        }
        cout << "    Isolated boundary nodes (no connections): " << isolated_nodes << endl;

        // 警告
        if (edge_count + large_net_edges == 0) {
            cout << "\n  ⚠️  WARNING: No edges between boundary nodes!" << endl;
            cout << "  This means boundary nodes are not connected to each other." << endl;
            cout << "  Clustering will not be effective." << endl;
            cout << "  Suggestions:" << endl;
            cout << "    - Check if boundary threshold is too strict" << endl;
            cout << "    - Verify netlist file is correct" << endl;
            cout << "    - Consider including non-boundary nodes in connectivity" << endl;
        }
    }

    // 获取一阶邻居
    unordered_set<int> getFirstOrderNeighbors(int node) {
        unordered_set<int> neighbors;

        if (adj_list.find(node) != adj_list.end()) {
            for (int n : adj_list[node]) {
                neighbors.insert(n);
            }
        }

        return neighbors;
    }

    // 计算两个节点（或超节点）的关联度
    // 使用所有点的连接信息，只考虑一阶邻居
    int computeRawWeight(const unordered_set<int>& cluster_u, const unordered_set<int>& cluster_v) {
        int weight = 0;

        // 计算两个集合之间的连接数
        for (int u : cluster_u) {
            if (adj_list.find(u) == adj_list.end()) continue;

            for (int v : cluster_v) {
                // 如果u和v直接相连
                if (adj_list[u].count(v)) {
                    weight++;
                }
            }
        }

        return weight;
    }

    // 改进的聚类：边界隔离 + 簇大小控制
    vector<unordered_set<int>> improvedClustering(int target_cluster_count = 1000, int max_cluster_size = 50) {
        cout << "\nStarting improved clustering..." << endl;

        int num_boundary_nodes = idx_to_inst.size();
        cout << "  Boundary nodes: " << num_boundary_nodes << endl;
        cout << "  Target cluster count: " << target_cluster_count << endl;
        cout << "  Max cluster size: " << max_cluster_size << endl;

        // 如果节点数已经小于目标，直接返回
        if (num_boundary_nodes <= target_cluster_count) {
            cout << "  ℹ️  Boundary nodes (" << num_boundary_nodes
                 << ") <= target (" << target_cluster_count << ")" << endl;
            cout << "  No clustering needed, each node forms its own cluster." << endl;

            vector<unordered_set<int>> result;
            for (int i = 0; i < num_boundary_nodes; i++) {
                unordered_set<int> cluster;
                cluster.insert(i);
                result.push_back(cluster);
            }
            return result;
        }

        // 步骤1: 收集边，只保留同边界的节点对
        cout << "  Step 1: Collecting edges (same boundary only)..." << endl;

        map<pair<int,int>, int> edge_weights;
        int filtered_out = 0;  // 统计被过滤掉的跨边界边

        for (const auto& [node_i, neighbors] : adj_list) {
            int boundary_i = idx_boundary_id[node_i];

            for (int node_j : neighbors) {
                if (node_i >= node_j) continue;  // 避免重复

                int boundary_j = idx_boundary_id[node_j];

                // 只保留同边界的边
                if (boundary_i == boundary_j) {
                    edge_weights[{node_i, node_j}]++;
                } else {
                    filtered_out++;
                }
            }
        }

        cout << "  Found " << edge_weights.size() << " same-boundary edges" << endl;
        cout << "  Filtered out " << filtered_out << " cross-boundary edges" << endl;

        if (edge_weights.empty()) {
            cout << "  ⚠️  WARNING: No same-boundary edges found!" << endl;
            vector<unordered_set<int>> result;
            for (int i = 0; i < num_boundary_nodes; i++) {
                unordered_set<int> cluster;
                cluster.insert(i);
                result.push_back(cluster);
            }
            return result;
        }

        // 步骤2: 按权重排序
        cout << "  Step 2: Sorting edges by weight..." << endl;

        vector<tuple<int, int, int>> sorted_edges;
        for (const auto& [nodes, weight] : edge_weights) {
            sorted_edges.push_back({nodes.first, nodes.second, weight});
        }

        sort(sorted_edges.begin(), sorted_edges.end(),
             [](const auto& a, const auto& b) {
                 return get<2>(a) > get<2>(b);
             });

        cout << "  Sorted " << sorted_edges.size() << " edges" << endl;

        // 步骤3: 带严格大小限制的并查集合并
        cout << "  Step 3: Merging with strict size limit..." << endl;

        // 并查集初始化
        vector<int> parent(num_boundary_nodes);
        vector<int> cluster_size(num_boundary_nodes, 1);  // 每个簇的实际大小

        for (int i = 0; i < num_boundary_nodes; i++) {
            parent[i] = i;
        }

        // 并查集查找（路径压缩）
        function<int(int)> find_root = [&](int x) {
            if (parent[x] != x) {
                parent[x] = find_root(parent[x]);
            }
            return parent[x];
        };

        // 获取簇的实际大小
        auto get_cluster_size = [&](int x) -> int {
            return cluster_size[find_root(x)];
        };

        // 按顺序合并，严格检查大小
        int num_clusters = num_boundary_nodes;
        int merged_count = 0;
        int rejected_size = 0;

        for (const auto& [node_i, node_j, weight] : sorted_edges) {
            // 优先检查：达到目标
            if (num_clusters <= target_cluster_count) {
                cout << "  ✓ Reached target cluster count. Stopping." << endl;
                break;
            }

            // 查找当前的根
            int root_i = find_root(node_i);
            int root_j = find_root(node_j);

            // 已在同簇
            if (root_i == root_j) continue;

            // 严格检查大小限制
            int size_i = cluster_size[root_i];
            int size_j = cluster_size[root_j];

            if (size_i + size_j > max_cluster_size) {
                rejected_size++;
                continue;
            }

            // 执行合并（小树挂在大树下，保持平衡）
            if (size_i < size_j) {
                parent[root_i] = root_j;
                cluster_size[root_j] += cluster_size[root_i];
                cluster_size[root_i] = 0;  // 清零已合并的
            } else {
                parent[root_j] = root_i;
                cluster_size[root_i] += cluster_size[root_j];
                cluster_size[root_j] = 0;  // 清零已合并的
            }

            num_clusters--;
            merged_count++;

            if (merged_count % 100 == 0) {
                cout << "    Merged " << merged_count
                     << ", clusters: " << num_clusters << "\r" << flush;
            }

            // 再次检查
            if (num_clusters <= target_cluster_count) {
                cout << endl;
                cout << "  ✓ Reached target cluster count. Stopping." << endl;
                break;
            }
        }

        cout << endl;
        cout << "  Merging completed:" << endl;
        cout << "    Successful merges: " << merged_count << endl;
        cout << "    Rejected (size limit): " << rejected_size << endl;
        cout << "    Final cluster count: " << num_clusters << endl;

        // 步骤4: 收集结果并验证
        cout << "  Step 4: Collecting results and verifying..." << endl;

        unordered_map<int, unordered_set<int>> clusters_map;

        for (int i = 0; i < num_boundary_nodes; i++) {
            int root = find_root(i);
            clusters_map[root].insert(i);
        }

        vector<unordered_set<int>> result;
        for (const auto& [root, cluster] : clusters_map) {
            result.push_back(cluster);
        }

        // 验证：检查是否有超过限制的簇
        int violations = 0;
        int max_actual_size = 0;

        for (const auto& cluster : result) {
            if (cluster.size() > max_cluster_size) {
                violations++;
                max_actual_size = max(max_actual_size, (int)cluster.size());
            }
        }

        if (violations > 0) {
            cout << "  ⚠️  WARNING: " << violations << " clusters exceed size limit!" << endl;
            cout << "      Largest cluster: " << max_actual_size << " nodes" << endl;
            cout << "      Splitting oversized clusters..." << endl;

            // 分割超大的簇
            vector<unordered_set<int>> fixed_result;

            for (const auto& cluster : result) {
                if (cluster.size() <= max_cluster_size) {
                    // 大小合适，直接保留
                    fixed_result.push_back(cluster);
                } else {
                    // 超大簇，需要分割
                    vector<int> nodes(cluster.begin(), cluster.end());

                    for (size_t i = 0; i < nodes.size(); i += max_cluster_size) {
                        unordered_set<int> sub_cluster;
                        for (size_t j = i; j < min(i + max_cluster_size, nodes.size()); j++) {
                            sub_cluster.insert(nodes[j]);
                        }
                        fixed_result.push_back(sub_cluster);
                    }
                }
            }

            result = fixed_result;
            cout << "      After splitting: " << result.size() << " clusters" << endl;
        } else {
            cout << "   All clusters within size limit." << endl;
        }

        // 统计
        vector<int> sizes;
        int single_node_clusters = 0;
        int max_size_clusters = 0;

        for (const auto& cluster : result) {
            sizes.push_back(cluster.size());
            if (cluster.size() == 1) single_node_clusters++;
            if (cluster.size() == max_cluster_size) max_size_clusters++;
        }
        sort(sizes.begin(), sizes.end());

        if (!sizes.empty()) {
            cout << "  Cluster size statistics:" << endl;
            cout << "    Min: " << sizes.front() << endl;
            cout << "    Max: " << sizes.back() << endl;
            cout << "    Median: " << sizes[sizes.size()/2] << endl;
            cout << "    Average: " << (num_boundary_nodes * 1.0 / result.size()) << endl;
            cout << "    Single-node clusters: " << single_node_clusters
                 << " (" << (single_node_clusters * 100.0 / result.size()) << "%)" << endl;
            cout << "    Max-size clusters: " << max_size_clusters
                 << " (" << (max_size_clusters * 100.0 / result.size()) << "%)" << endl;
        }

        // 按边界统计
        vector<int> boundary_cluster_counts(3, 0);
        for (const auto& cluster : result) {
            int boundary_id = idx_boundary_id[*cluster.begin()];
            boundary_cluster_counts[boundary_id]++;
        }

        cout << "  Clusters per boundary:" << endl;
        cout << "    die0|1: " << boundary_cluster_counts[0] << " clusters" << endl;
        cout << "    die1|2: " << boundary_cluster_counts[1] << " clusters" << endl;
        cout << "    die2|3: " << boundary_cluster_counts[2] << " clusters" << endl;

        return result;
    }

    // 输出超节点到文件
    void outputSupernodes(const vector<unordered_set<int>>& supernodes, const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Cannot create output file: " << filename << endl;
            return;
        }

        for (size_t i = 0; i < supernodes.size(); i++) {
            file << "Supernode_" << i << " (" << supernodes[i].size() << " nodes):" << endl;
            for (int idx : supernodes[i]) {
                const auto& inst = instances[idx_to_inst[idx]];
                file << "  " << inst.name
                     << " (x=" << inst.x << ", y=" << inst.y << ", die=" << inst.die << ")" << endl;
            }
            file << endl;
        }

        file.close();
        cout << "\nSupernodes written to file: " << filename << endl;
    }

    // 完整粗化流程
    void coarsen(double boundary_threshold = 15.0, int target_cluster_count = 1000,
                 int max_net_size = 100, int max_cluster_size = 50,
                 const string& output_file = "supernodes.txt") {
        cout << "\n========================================" << endl;
        cout << "Graph Coarsening - FPGA Netlist (Die Boundary)" << endl;
        cout << "========================================" << endl;

        // 1. 先筛选边界节点并标记所属边界
        filterBoundaryNodes(boundary_threshold);

        if (boundary_node_names.empty()) {
            cout << "Warning: No boundary nodes!" << endl;
            return;
        }

        // 2. 只为边界节点建图
        buildGraphForBoundaryNodes(max_net_size);

        // 3. 执行改进的聚类（边界隔离 + 大小控制）
        auto supernodes = improvedClustering(target_cluster_count, max_cluster_size);

        // 4. 输出部分结果到屏幕
        cout << "\nFirst 5 supernodes example:" << endl;
        for (size_t i = 0; i < min(size_t(5), supernodes.size()); i++) {
            cout << "  Supernode " << i << ": contains " << supernodes[i].size() << " instances" << endl;

            // 显示边界ID
            int boundary_id = idx_boundary_id[*supernodes[i].begin()];
            cout << "    Boundary: die" << boundary_id << "|" << (boundary_id+1)
                 << " (y=" << die_boundaries[boundary_id] << ")" << endl;

            // 统计这个超节点中的die分布
            vector<int> die_dist(4, 0);
            for (int idx : supernodes[i]) {
                int die = instances[idx_to_inst[idx]].die;
                if (die >= 0 && die < 4) die_dist[die]++;
            }

            cout << "    Die distribution: ";
            for (int d = 0; d < 4; d++) {
                if (die_dist[d] > 0) {
                    cout << "die" << d << "=" << die_dist[d] << " ";
                }
            }
            cout << endl;

            // 显示前3个实例
            cout << "    Example instances: ";
            int count = 0;
            for (int idx : supernodes[i]) {
                if (count++ >= 3) break;
                cout << idx_to_inst[idx] << " ";
            }
            cout << "..." << endl;
        }

        // 5. 输出到文件
        if (!output_file.empty()) {
            outputSupernodes(supernodes, output_file);
        }
    }
};

// 主函数
int main() {
    // 文件路径
    std::string placement_file = "C:\\Users\\Lenovo\\Desktop\\plresults\\FPGA12.pl";
    std::string netlist_file = "C:\\Users\\Lenovo\\Desktop\\ispd2016_flexshelf\\FPGA12\\design.osv";

    NetlistGraphCoarsening reader;

    // 提取文件名
    std::filesystem::path p(placement_file);
    reader.filename = p.filename().string();

    std::cout << "Processing file: " << reader.filename << std::endl;

    // 读取文件
    if (!reader.readPlaceFile(placement_file) || !reader.readNetFile(netlist_file)) {
        std::cerr << "File read error\n";
        return 1;
    }

    // 执行粗化（内部会自动筛选边界节点并建图）
    // boundary_threshold: 边界距离阈值 (<=15)
    // target_cluster_count: 目标簇数量 (<1000 停止合并)
    // max_net_size: net大小阈值，超过则用星形拓扑
    // max_cluster_size: 簇大小上限（控制均匀性）
    string output_file = "supernodes_" + reader.filename + ".txt";
    reader.coarsen(15.0, 1000, 100, 50, output_file);

    return 0;
}