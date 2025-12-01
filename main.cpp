#include "src/coarsing.cpp"

// 主函数
int main() {
    std::string placement_file;
    std::string netlist_file;
    std::string output_dir = "../";

    placement_file = "C:\\Users\\Lenovo\\Desktop\\plresults\\FPGA01.pl";
    netlist_file = "C:\\Users\\Lenovo\\Desktop\\ispd2016_flexshelf\\FPGA01\\design-without.osv";


    // 确保输出目录存在
    std::filesystem::create_directories(output_dir);

    // ========== 创建 coarsening 对象 ==========
    NetlistGraphCoarsening reader;

    // 提取文件名
    std::filesystem::path p(placement_file);
    reader.filename = p.filename().string();

    reader.is_coarsing = false;
    if (!reader.readPlaceFile(placement_file)) {
        std::cerr << "ERROR: Failed to read placement file: " << placement_file << std::endl;
        return 1;
    }

    if (!reader.readNetFile(netlist_file)) {
        std::cerr << "ERROR: Failed to read netlist file: " << netlist_file << std::endl;
        return 1;
    }
    double boundary_threshold = 15.0;
    int target_cluster_count = 1000;
    int max_net_size = 100;
    int max_cluster_size = 50;

    auto supernodes = reader.coarsen(
            boundary_threshold,
            target_cluster_count,
            max_net_size,
            max_cluster_size
    );

    if (supernodes.empty()) {
        std::cerr << "ERROR: Coarsening failed - no supernodes generated" << std::endl;
        return 1;
    }


    // 构建输出文件名
    std::string base_name = reader.filename.substr(0, reader.filename.find_last_of('.'));
    std::string coarsened_pl = output_dir + "coarsened_" + base_name + ".pl";
    std::string coarsened_net = output_dir + "coarsened_" + base_name + ".net";
    std::string detailed_file = output_dir + "supernodes_" + base_name + ".txt";

    // 生成 fpga_reader 需要的文件（必需）
    // 注意：使用 outputForFPGAReader 而不是 outputForMoveProgram
    reader.outputForFPGAReader(supernodes, coarsened_pl, coarsened_net);

    // 可选：生成详细的超节点列表（用于查看和调试）
    bool output_detailed = false;  // 设置为 true 以生成详细列表
    if (output_detailed) {
        reader.outputSupernodes(supernodes, detailed_file);
        std::cout << "\n✓ Detailed supernodes list: " << detailed_file << std::endl;
    }


    if (output_detailed) {
        std::cout << "\n  3. " << detailed_file << std::endl;
        std::cout << "     Purpose: Detailed supernodes information (for debugging)" << std::endl;
    }

    return 0;
}