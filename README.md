项目简介
NetlistGraphCoarsening 是一个用于大规模 FPGA 网表粗化的预处理工具。它主要针对跨 Die 边界的节点
进行聚类（Clustering）， 将紧密连接的实例打包成“超节点（Supernode）”，从而降低后续布局或划分算法的复杂度。

输入文件格式
1. 布局文件 (.pl)
   工具根据 is_coarsing 参数的不同，支持两种输入格式：

        当 is_coarsing = true 时（读取已粗化过的文件）：
        格式：实例名 y坐标 Die索引 [FIXED]
        示例：
            Plaintext
            inst_1 65 0
            inst_2 130 1 FIXED
        注意：FIXED 标签表示该实例在粗化过程中不可移动且不会被合并。
        
        当 is_coarsing = false 时（读取原始布局文件）：
        格式：实例名 x坐标 y坐标 z坐标 [FIXED]
        示例：
            Plaintext
            inst_1 10 20 30
            inst_2 10 20 30 FIXED

2. 网表文件 (.osv)
   格式：Verilog 风格的实例声明。
支持：支持复杂的引脚定义，包括数组索引解析（如 .A[3:0](net_name) 或 .B[6](net)）。

核心参数说明 (主函数控制)
参数名                   类型     默认值          说明
is_coarsing             bool    false           判断传入的布局文件是否为已经粗化过的文件。
boundary_threshold      double   15.0           边界阈值：实例距离 
Die                     边界  （120/240/360）     在此范围内时，才会被视为“边界节点”并参与粗化。
target_cluster_count    int     1000            目标簇数量：期望生成的超节点总数。算法会根据各边界区域的节点密度自动分配聚类指标。
max_cluster_size        int     50              簇大小上限：单个超节点（Supernode）内允许包含的最大原始实例数量。
max_net_size            int     100             大网过滤阈值：连接数超过此值的 Net 将采用星形连接（Hub Model）而非全连接，以降低图构建时的内存和计算开销。

输出文件定义执行完成后，工具会生成以下三类文件：
粗化布局文件 (coarsened_*.pl)
    包含所有未合并的原始实例（Inner Nodes）和新生成的 supernode_i。
    超节点的坐标被重置到其所属 Die 的中心位置（例如 Die 0 的 y 坐标设为 60.0）。
    保留原始的 FIXED 属性。

粗化网表文件 (coarsened_*.net)
重新定义的模块化网表。
关键特性：Supernode 内部的连线（即该 Net 连接的所有点都在同一个 Supernode 内）
    会被自动隐藏/剔除，不再出现在端口列表中，确保后续 Cut 计算准确。

调试报告 (supernodes_*.txt)（可选）
详细列出每个 Supernode 包含的具体原始实例名称及其物理属性（坐标、Die ID）。