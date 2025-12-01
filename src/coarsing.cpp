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
#include <random>
#include <chrono>

using namespace std;

struct Instance {
    string name;
    string type;
    vector<string> nets;
    double x, y;
    int die;
    bool is_fixed;
    Instance() : x(0), y(0), die(-1), is_fixed(false) {}
};

// 图粗化类
class NetlistGraphCoarsening {
private:
    unordered_map<string, Instance> instances;
    unordered_map<string, Instance> fixed_instances;
    unordered_map<string, vector<string>> net_to_insts;
    unordered_map<string, int> inst_to_idx;
    vector<string> idx_to_inst;

    unordered_map<int, unordered_set<int>> adj_list;

    unordered_set<string> boundary_node_names;
    unordered_map<string, int> node_boundary_id;
    unordered_map<int, int> idx_boundary_id;

    vector<double> die_boundaries = {120.0, 240.0, 360.0};

public:
    string filename;
    bool is_coarsing = false;
    NetlistGraphCoarsening() {}

    bool readPlaceFile(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Cannot open placement file: " << filename << endl;
            return false;
        }

        string line;
        int loaded = 0;
        int fixed_count = 0;

        while (getline(file, line)) {
            if (line.empty()) continue;

            istringstream iss(line);
            string inst_name, token2, token3, token4, token5;

            if (!(iss >> inst_name >> token2)) continue;

            Instance inst;
            inst.name = inst_name;
            bool success = false;

            if (iss >> token3) {
                if (!is_coarsing) {
                    try {
                        double y = stod(token3);
                        inst.y = y;

                        if (y < 120) {
                            inst.die = 0;
                        } else if (y < 240) {
                            inst.die = 1;
                        } else if (y < 360) {
                            inst.die = 2;
                        } else {
                            inst.die = 3;
                        }

                        if (iss >> token5 && token5 == "FIXED") {
                            inst.is_fixed = true;
                        }

                        success = true;
                    } catch (...) {
                        cout << "read place file failed" << endl;
                    }
                } else {
                    try {
                        if (token3 == "FIXED") {
                            continue;
                        }

                        int die = stoi(token3);
                        inst.die = die;

                        inst.y = die * 120.0 + 60.0;

                        if (iss >> token4 && token4 == "FIXED") {
                            inst.is_fixed = true;
                        }
                        success = true;
                    } catch (...) {
                        cout << "read place file failed" << endl;
                    }
                }
            }

            if (success) {
                if (inst.is_fixed) {
                    fixed_instances[inst_name] = inst;
                    fixed_count++;
                } else {
                    instances[inst_name] = inst;
                }
                loaded++;
            }
        }

        file.close();
        cout << "Placement loading completed! Loaded " << loaded << " instances" << endl;
        return true;
    }

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
        return true;
    }

    // 保存实例信息
    void saveInstance(const string& name, const string& type, const vector<string>& nets) {
        // 【修改】只保存非FIXED实例的类型和nets
        if (instances.find(name) != instances.end()) {
            instances[name].type = type;
            instances[name].nets = nets;

            // 建立net到instance的反向映射（只针对非FIXED实例）
            for (const auto& net : nets) {
                net_to_insts[net].push_back(name);
            }
        } else if (fixed_instances.find(name) != fixed_instances.end()) {
            // FIXED 实例也更新类型和nets，但不参与图构建
            fixed_instances[name].type = type;
            fixed_instances[name].nets = nets;
        } else {
            // 如果实例既不在 instances 也不在 fixed_instances 中
            // 创建新实例（默认非FIXED）
            Instance inst;
            inst.name = name;
            inst.type = type;
            inst.nets = nets;
            inst.is_fixed = false;
            instances[name] = inst;

            for (const auto& net : nets) {
                net_to_insts[net].push_back(name);
            }
        }
    }


    void filterBoundaryNodes(double threshold = 15.0) {

        boundary_node_names.clear();
        node_boundary_id.clear();

        vector<int> boundary_counts(3, 0);  // 每个边界的节点数

        for (const auto& [name, inst] : instances) {  // 只考虑非FIXED实例
            double y = inst.y;

            // 检查每个边界
            for (size_t i = 0; i < die_boundaries.size(); i++) {
                double boundary = die_boundaries[i];
                double dist = abs(y - boundary);

                if (dist <= threshold) {
                    boundary_node_names.insert(name);
                    node_boundary_id[name] = i;
                    boundary_counts[i]++;
                    break;  // 一个节点只属于一个边界
                }
            }
        }
    }

    // 只为边界节点构建图
    void buildGraphForBoundaryNodes(int max_net_size = 100) {
        cout << "\nBuilding graph for boundary nodes only..." << endl;

        // 清空现有的索引和邻接表
        inst_to_idx.clear();
        idx_to_inst.clear();
        adj_list.clear();
        idx_boundary_id.clear();

        // 只为边界节点建立索引
        int idx = 0;
        for (const string& name : boundary_node_names) {
            inst_to_idx[name] = idx;
            idx_to_inst.push_back(name);

            // 同时记录这个idx对应的边界ID
            idx_boundary_id[idx] = node_boundary_id[name];

            idx++;
        }


        int edge_count = 0;
        int large_net_edges = 0;

        for (const auto& [net, inst_list] : net_to_insts) {

            vector<string> boundary_insts;
            for (const auto& inst_name : inst_list) {
                if (boundary_node_names.count(inst_name)) {
                    boundary_insts.push_back(inst_name);
                }
            }

            size_t net_size = boundary_insts.size();
            if (net_size <= 1) continue;

            // 限制：只连接同一边界的节点（避免跨边界聚类）
            unordered_map<int, vector<string>> boundary_groups;
            for (const auto& inst_name : boundary_insts) {
                int boundary_id = node_boundary_id[inst_name];
                boundary_groups[boundary_id].push_back(inst_name);
            }

            for (const auto& [boundary_id, group] : boundary_groups) {
                size_t group_size = group.size();
                if (group_size <= 1) continue;

                if (group_size > (size_t)max_net_size) {

                    int hub_idx = inst_to_idx[group[0]];
                    for (size_t i = 1; i < group_size; i++) {
                        int idx = inst_to_idx[group[i]];
                        adj_list[hub_idx].insert(idx);
                        adj_list[idx].insert(hub_idx);
                        large_net_edges++;
                    }
                } else {
                    for (size_t i = 0; i < group_size; i++) {
                        for (size_t j = i + 1; j < group_size; j++) {
                            int idx_i = inst_to_idx[group[i]];
                            int idx_j = inst_to_idx[group[j]];

                            adj_list[idx_i].insert(idx_j);
                            adj_list[idx_j].insert(idx_i);
                            edge_count++;
                        }
                    }
                }
            }
        }

    }

    // 改进的聚类算法
    vector<unordered_set<int>> improvedClustering(int target_cluster_count = 1000,
                                                  int max_cluster_size = 50) {

        int num_boundary_nodes = idx_to_inst.size();
        if (num_boundary_nodes == 0) {
            cout << "Warning: No boundary nodes to cluster!" << endl;
            return vector<unordered_set<int>>();
        }

        vector<vector<int>> boundary_node_lists(3);
        for (int idx = 0; idx < num_boundary_nodes; idx++) {
            int boundary_id = idx_boundary_id[idx];
            boundary_node_lists[boundary_id].push_back(idx);
        }

        vector<int> target_clusters_per_boundary(3);
        for (int b = 0; b < 3; b++) {
            target_clusters_per_boundary[b] = max(1,
                                                  (int)(target_cluster_count * boundary_node_lists[b].size() / (double)num_boundary_nodes));
        }


        vector<unordered_set<int>> result;

        for (int boundary_id = 0; boundary_id < 3; boundary_id++) {
            const auto& nodes = boundary_node_lists[boundary_id];
            if (nodes.empty()) continue;

            int target = target_clusters_per_boundary[boundary_id];
            auto clusters = clusterBoundary(nodes, target, max_cluster_size);

            result.insert(result.end(), clusters.begin(), clusters.end());
        }

        bool need_split = false;
        for (const auto& cluster : result) {
            if ((int)cluster.size() > max_cluster_size) {
                need_split = true;
                break;
            }
        }

        if (need_split) {
            vector<unordered_set<int>> fixed_result;

            for (const auto& cluster : result) {
                if ((int)cluster.size() <= max_cluster_size) {
                    fixed_result.push_back(cluster);
                } else {
                    vector<int> nodes(cluster.begin(), cluster.end());
                    int num_splits = (nodes.size() + max_cluster_size - 1) / max_cluster_size;

                    for (int i = 0; i < num_splits; i++) {
                        unordered_set<int> small_cluster;
                        int start = i * max_cluster_size;
                        int end = min(start + max_cluster_size, (int)nodes.size());

                        for (int j = start; j < end; j++) {
                            small_cluster.insert(nodes[j]);
                        }
                        fixed_result.push_back(small_cluster);
                    }
                }
            }

            result = fixed_result;
            cout << "      After splitting: " << result.size() << " clusters" << endl;
        } else {
            cout << "   All clusters within size limit." << endl;
        }

        vector<int> sizes;
        int single_node_clusters = 0;
        int max_size_clusters = 0;

        for (const auto& cluster : result) {
            sizes.push_back(cluster.size());
            if (cluster.size() == 1) single_node_clusters++;
            if (cluster.size() == max_cluster_size) max_size_clusters++;
        }
        sort(sizes.begin(), sizes.end());

        return result;
    }

    vector<unordered_set<int>> clusterBoundary(const vector<int>& nodes,
                                               int target_count,
                                               int max_size) {
        if (nodes.empty()) return vector<unordered_set<int>>();

        int n = nodes.size();
        int target_cluster_size = max(1, n / target_count);

        if (target_count >= n) {
            vector<unordered_set<int>> result;
            for (int idx : nodes) {
                result.push_back({idx});
            }
            return result;
        }

        vector<bool> visited(n, false);
        vector<unordered_set<int>> clusters;

        vector<int> shuffled_nodes = nodes;
        unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
        shuffle(shuffled_nodes.begin(), shuffled_nodes.end(), default_random_engine(seed));
        // BFS聚类
        for (int start_idx : shuffled_nodes) {
            if (visited[inst_to_idx[idx_to_inst[start_idx]] - inst_to_idx[idx_to_inst[nodes[0]]]]) {
                continue;
            }

            unordered_set<int> cluster;
            queue<int> q;
            q.push(start_idx);

            while (!q.empty() && (int)cluster.size() < target_cluster_size) {
                int curr = q.front();
                q.pop();

                int curr_pos = -1;
                for (int i = 0; i < n; i++) {
                    if (nodes[i] == curr) {
                        curr_pos = i;
                        break;
                    }
                }

                if (curr_pos == -1 || visited[curr_pos]) continue;

                visited[curr_pos] = true;
                cluster.insert(curr);

                // 将邻居加入队列
                if (adj_list.count(curr)) {
                    for (int neighbor : adj_list[curr]) {
                        bool is_boundary_neighbor = false;
                        for (int node : nodes) {
                            if (node == neighbor) {
                                is_boundary_neighbor = true;
                                break;
                            }
                        }

                        if (is_boundary_neighbor) {
                            int neighbor_pos = -1;
                            for (int i = 0; i < n; i++) {
                                if (nodes[i] == neighbor) {
                                    neighbor_pos = i;
                                    break;
                                }
                            }

                            if (neighbor_pos != -1 && !visited[neighbor_pos]) {
                                q.push(neighbor);
                            }
                        }
                    }
                }
            }

            if (!cluster.empty()) {
                clusters.push_back(cluster);
            }
        }

        for (int i = 0; i < n; i++) {
            if (!visited[i]) {
                if (!clusters.empty()) {
                    auto min_it = min_element(clusters.begin(), clusters.end(),
                                              [](const unordered_set<int>& a, const unordered_set<int>& b) {
                                                  return a.size() < b.size();
                                              });

                    if ((int)min_it->size() < max_size) {
                        min_it->insert(nodes[i]);
                    } else {
                        clusters.push_back({nodes[i]});
                    }
                } else {
                    clusters.push_back({nodes[i]});
                }
            }
        }

        return clusters;
    }

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
                     << " (y=" << inst.y << ", die=" << inst.die << ")" << endl;
            }
            file << endl;
        }

        file.close();
        cout << "\nDetailed supernodes sucess "<< endl;
    }

    void outputForFPGAReader(const vector<unordered_set<int>>& supernodes,
                             const string& placement_file,
                             const string& netlist_file) {

        unordered_set<string> coarsened_inst_names;
        for (const auto& cluster : supernodes) {
            for (int idx : cluster) {
                coarsened_inst_names.insert(idx_to_inst[idx]);
            }
        }

        // 2. 找出“内部点”(Inner Nodes)：即非 Fixed 且 未被粗化 的点
        vector<string> inner_node_names;
        for (const auto& [name, inst] : instances) {
            if (coarsened_inst_names.find(name) == coarsened_inst_names.end()) {
                inner_node_names.push_back(name);
            }
        }
        sort(inner_node_names.begin(), inner_node_names.end());

        vector<string> fixed_node_names;
        for (const auto& [name, inst] : fixed_instances) {
            fixed_node_names.push_back(name);
        }
        sort(fixed_node_names.begin(), fixed_node_names.end());

        ofstream pl_file(placement_file);
        if (!pl_file.is_open()) {
            cerr << "Cannot create placement file: " << placement_file << endl;
            return;
        }

        for (const string& name : fixed_node_names) {
            if (instances.count(name) == 0) {
                cerr << "WARNING: Inner node instance '" << name << "' not found in 'instances' map. Skipping." << endl;
                continue;
            }
            const auto& inst = instances.at(name);
            pl_file << name << " " << (int)inst.y << " " << inst.die << " FIXED" << endl;
        }

        for (const string& name : inner_node_names) {
            const auto& inst = instances.at(name);
            if(inst.die == -1) continue;
            pl_file << name << " " << (int)inst.y << " " << inst.die << endl;
        }

        for (size_t i = 0; i < supernodes.size(); i++) {
            int die = -1;
            if (!supernodes[i].empty()) {
                die = instances[idx_to_inst[*supernodes[i].begin()]].die;
            }
            pl_file << "supernode_" << i << die * 60 + 120.0<< die << endl;
        }
        pl_file.close();
        ofstream net_file(netlist_file);
        if (!net_file.is_open()) {
            cerr << "Cannot create netlist file: " << netlist_file << endl;
            return;
        }

        unordered_map<string, string> inst_to_cluster_name;

        for (size_t i = 0; i < supernodes.size(); i++) {
            string s_name = "supernode_" + to_string(i);
            for (int idx : supernodes[i]) {
                inst_to_cluster_name[idx_to_inst[idx]] = s_name;
            }
        }
        for (const string& name : fixed_node_names) {
            inst_to_cluster_name[name] = name;
        }
        for (const string& name : inner_node_names) {
            inst_to_cluster_name[name] = name;
        }

        auto write_module = [&](const string& cluster_name,
                                const vector<string>& internal_insts,
                                const string& cluster_type) {
            unordered_set<string> boundary_nets;

            for (const string& inst_name : internal_insts) {
                // 获取实例引用 (可能在 instances 或 fixed_instances 中)
                const vector<string>* nets_ptr = nullptr;
                if (instances.count(inst_name)) nets_ptr = &instances.at(inst_name).nets;
                else if (fixed_instances.count(inst_name)) nets_ptr = &fixed_instances.at(inst_name).nets;

                if (!nets_ptr) continue;

                for (const string& net : *nets_ptr) {
                    // 判断该 net 是否连接到当前 Cluster 外部
                    bool is_external = false;

                    if (net_to_insts.count(net)) {
                        for (const string& connected_inst : net_to_insts[net]) {
                            // 如果连接的实例属于不同的 Cluster (或未找到Cluster映射)，则为外部Net
                            auto it = inst_to_cluster_name.find(connected_inst);
                            if (it == inst_to_cluster_name.end() || it->second != cluster_name) {
                                is_external = true;
                                break;
                            }
                        }
                    }
                    if (is_external) boundary_nets.insert(net);
                }
            }

            // 【主要修改】在 cluster_name 前输出 cluster_type
            net_file << cluster_type << " " << cluster_name << " (" << endl;
            vector<string> sorted_nets(boundary_nets.begin(), boundary_nets.end());
            sort(sorted_nets.begin(), sorted_nets.end());

            for (size_t k = 0; k < sorted_nets.size(); k++) {
                net_file << "  .PIN" << k << "(" << sorted_nets[k] << ")";
                if (k < sorted_nets.size() - 1) net_file << ",";
                net_file << endl;
            }
            net_file << ");" << endl << endl;
        };
        for (const string& name : fixed_node_names) {
            const auto& inst = fixed_instances.at(name);
            write_module(name, {name}, inst.type); // 使用其自身的类型
        }

        for (const string& name : inner_node_names) {
            const auto& inst = instances.at(name);
            write_module(name, {name}, inst.type); // 使用其自身的类型
        }

        for (size_t i = 0; i < supernodes.size(); i++) {
            string s_name = "supernode_" + to_string(i);
            vector<string> members;
            for (int idx : supernodes[i]) {
                members.push_back(idx_to_inst[idx]);
            }

            string s_type = getClusterType(members);

            write_module(s_name, members, s_type);
        }

        net_file.close();
    }

    string getClusterType(const vector<string>& internal_insts) {
        if (internal_insts.empty()) return "MODULE";

        string first_inst_type = "";
        string lut_type = "";
        string ram_type = "";

        for (const string& inst_name : internal_insts) {
            const Instance* inst_ptr = nullptr;
            if (instances.count(inst_name)) inst_ptr = &instances.at(inst_name);
            else if (fixed_instances.count(inst_name)) inst_ptr = &fixed_instances.at(inst_name);

            if (!inst_ptr) continue;

            if (first_inst_type.empty()) {
                first_inst_type = inst_ptr->type;
            }

            if (inst_ptr->type.find("LUT") != string::npos) {
                lut_type = inst_ptr->type;
            } else if (inst_ptr->type.find("RAM") != string::npos) {
                ram_type = inst_ptr->type;
            }
        }

        if (!ram_type.empty()) return ram_type;
        if (!lut_type.empty()) return lut_type;
        if (!first_inst_type.empty()) return first_inst_type;

        return "MODULE";
    }


    vector<unordered_set<int>> coarsen(double boundary_threshold = 15.0, int target_cluster_count = 1000,
                                       int max_net_size = 100, int max_cluster_size = 50) {


        filterBoundaryNodes(boundary_threshold);

        if (boundary_node_names.empty()) {
            cout << "Warning: No boundary nodes!" << endl;
            return vector<unordered_set<int>>();
        }

        buildGraphForBoundaryNodes(max_net_size);

        auto supernodes =
                improvedClustering(target_cluster_count, max_cluster_size);
        return supernodes;
    }
};
