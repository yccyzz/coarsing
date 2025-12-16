1. 输入文件格式
工具输入文件：
    布局文件 (.pl)
        当参数is_coarsing = true时
            格式：每行包含 实例 y Die [FIXED]。
                示例：
                inst_1 65 0
                inst_2 130 1 FIXED
                注意：FIXED 标签表示该实例在粗化过程中不可移动且不会被合并。
        当参数is_coarsing = false时
            格式：每行包含 实例 x y z die [FIXED]。
               示例：
               inst_1 10 20 30 
               inst_2 10 20 30  FIXED
    网表文件 (.osv)
       格式：Verilog 风格的实例声明。
       支持：支持复杂的引脚定义，包括数组索引（如 .A[3:0](net_name)）。
2. 主函数中控制的参数
    is_coarsing,bool,false,
     判单传入文件是否为粗化后的文件
    boundary_threshold,double,15.0,
     边界阈值：实例距离 Die 边界（120/240/360）在此范围内时，才会被视为“边界节点”参与粗化。
    target_cluster_count,int,1000,
     目标簇数量：期望生成的超节点总数。算法会根据各边界节点密度自动分配指标。
    max_cluster_size,int,50,
     簇大小上限：单个超节点内允许包含的最大原始实例数量。
    max_net_size,int,100,
     大网过滤阈值：连接数超过此值的 Net 将采用星形连接（Hub Model）而非全连接，以降低图构建开销。
3. 输出文件定义
执行完成后，工具会生成以下三类文件：
    粗化布局文件 (coarsened_*.pl)：包含所有未合并的原始实例和新生成的 supernode_i。超节点坐标被重置到其所属 Die 的中心位置（例如 Die 0 的 $y$ 坐标设为 60.0）。
    粗化网表文件 (coarsened_*.net)： 重新定义的模块化网表。
    调试报告 (supernodes_*.txt)（可选）：详细列出每个超节点包含的具体原始实例名称及其物理属性。