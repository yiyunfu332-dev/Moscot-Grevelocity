# Zebrafish hematopoiesis workflow summary / 斑马鱼造血分析流程总结

## 1. Objective / 研究目标

**English.** The workflow estimates developmental transport, proliferation-associated growth, and directional cell-state transitions in zebrafish hematopoietic cells. Zebrafish genes are scored in zebrafish expression data; reviewed human orthologs are used only for cross-species interpretation.

**中文。** 本流程用于估计斑马鱼造血细胞发育过程中的细胞运输、增殖相关生长以及具有方向性的细胞状态转变。增殖评分使用斑马鱼表达数据中的斑马鱼基因；经过审核的人类直系同源基因仅用于跨物种解释。

## 2. Corrected species logic / 已修正的物种逻辑

**English.** The former `HUMAN_PROLIFERATION_MARKERS` block was incorrect and was removed. `all_genes.txt` contains 16,940 zebrafish features, but it is a candidate universe rather than a proliferation signature. The final proliferation panel contains only zebrafish genes supported by cell-cycle biology, presence in the full expression matrix, and curated zebrafish-to-human orthology.

**中文。** 原先的 `HUMAN_PROLIFERATION_MARKERS` 模块属于错误的物种标注，现已删除。`all_genes.txt` 含有 16,940 个斑马鱼基因特征，但它是候选基因全集，而不是增殖特征集。最终的增殖面板只保留满足以下条件的斑马鱼基因：具有细胞周期证据、存在于完整表达矩阵中，并且具有经过人工审核的斑马鱼到人类直系同源关系。

## 3. Two aligned expression objects / 两个按细胞 ID 对齐的表达对象

**English.** `zebrahub_full_velocity.h5ad` contains 120,800 cells but only 3,000 velocity-selected genes. It is used for Dynamo and GraphVelo. `zf_atlas_hematopoetic_endothelial_v4_release.h5ad` contains 5,441 hematopoietic/endothelial cells and 27,435 genes. Its cell IDs exactly overlap the relevant cells in the velocity and lineage objects, so it is used for proliferation and apoptosis scoring. This fixes the earlier result in which only `ube2c` was found.

**中文。** `zebrahub_full_velocity.h5ad` 含有 120,800 个细胞，但仅保留了 3,000 个用于速度分析的基因，因此用于 Dynamo 和 GraphVelo。`zf_atlas_hematopoetic_endothelial_v4_release.h5ad` 含有 5,441 个造血/内皮细胞和 27,435 个基因。它与速度对象及谱系对象中的相关细胞具有完全一致的细胞 ID，因此用于增殖和凋亡评分。该修正解决了先前只能匹配到 `ube2c` 的问题。

## 4. Cell alignment and annotation / 细胞对齐与注释

**English.** Velocity, lineage, time point, fish identity, and official fine cell type are joined only by exact cell IDs. No row-position matching or invented cluster-to-cell-type dictionary is used. Rare official labels are filtered by a minimum-cell threshold before modeling.

**中文。** 速度数据、谱系数据、时间点、个体编号和官方精细细胞类型仅通过完全一致的细胞 ID 进行连接。不使用按行位置匹配，也不人为构造聚类到细胞类型的映射字典。建模前根据最小细胞数阈值过滤样本量过少的官方标签。

## 5. Preprocessing and velocity / 预处理与 RNA 速度

**English.** Raw spliced and unspliced counts are reconstructed for aligned blood cells. Dynamo preprocessing, PCA, a 30-neighbor graph, moments, and stochastic RNA velocity are calculated once. Velocity genes are retained only after finite-value and quality-control checks. GraphVelo learns a manifold-constrained velocity and projects it into the shared PCA state space.

**中文。** 对齐后的血细胞重新构建 spliced 和 unspliced 原始计数。随后统一执行一次 Dynamo 预处理、PCA、30 近邻图、矩估计和随机 RNA 速度计算。速度基因需通过有限值和质量控制检查后才被保留。GraphVelo 学习受流形约束的速度，并将其投影到共享的 PCA 状态空间。

## 6. Reviewed marker mapping / 经审核的标志基因映射

**English.** The reviewed CSV stores `zebrafish_gene`, `human_gene`, mapping direction, database, and review status. Every retained row agrees with the local ZFIN manually curated human–zebrafish orthology table. Unsupported `bmf2 -> BMF` and expression-absent `POLD4` were removed. The current scoring table contains 38 proliferation and 15 apoptosis mappings. Teleost paralogs remain explicit, for example `baxa/baxb -> BAX`.

**中文。** 经审核的 CSV 保存 `zebrafish_gene`、`human_gene`、映射方向、数据库和审核状态。所有保留的记录均与本地 ZFIN 人工整理的人类—斑马鱼直系同源表一致。不受 ZFIN 支持的 `bmf2 -> BMF` 以及在完整评分矩阵中缺失的 `POLD4` 已被删除。当前评分表包含 38 个增殖映射和 15 个凋亡映射。硬骨鱼基因组复制产生的旁系同源基因被明确保留，例如 `baxa/baxb -> BAX`。

## 7. Growth score and real time / 生长评分与真实时间

Time is converted from hours post fertilization to days:

$$
t_{\mathrm{day}}=\frac{t_{\mathrm{hpf}}}{24}.
$$

For cell $i$, proliferation score $p_i$ and apoptosis score $a_i$ define the prior growth rate:

$$
g_i^{\mathrm{prior}}=
\exp\left(\frac{(p_i-a_i)(t_1-t_0)}{s}\right),
$$

where $s$ is the configured growth-scaling constant.

**中文。** 时间由受精后小时转换为天。细胞 $i$ 的增殖评分 $p_i$ 与凋亡评分 $a_i$ 共同定义先验生长率；$s$ 为生长缩放常数。

## 8. Transport models / 运输模型

**English.** M0 is a uniform-marginal moscot baseline. M1 introduces growth-informed unbalanced transport. If $P^{\mathrm{raw}}$ is the M1 coupling, source-cell transported mass is

$$
g_i=\sum_j P^{\mathrm{raw}}_{ij},
$$

and conditional fate is

$$
P^{\mathrm{cond}}_{ij}=
\frac{P^{\mathrm{raw}}_{ij}}{g_i}.
$$

**中文。** M0 是具有均匀边际的 moscot 基线。M1 加入生长信息并采用非平衡最优传输。原始耦合矩阵的行和表示源细胞运输质量，行归一化矩阵表示在发生运输条件下的命运概率。

## 9. GraphVelo reweighting / GraphVelo 方向重加权

Directional compatibility between source velocity $V_i$ and displacement to target $j$ is

$$
s_{ij}=\cos\left(V_i, X_j-X_i\right).
$$

M2 combines conditional transport and velocity direction:

$$
Q^{\mathrm{cond}}_{ij}\propto
\left(P^{\mathrm{cond}}_{ij}\right)^{\alpha}
\exp\left(\kappa s_{ij}\right).
$$

After row normalization, growth mass is restored:

$$
Q^{\mathrm{raw}}_{ij}=g_iQ^{\mathrm{cond}}_{ij}.
$$

**中文。** $s_{ij}$ 衡量源细胞速度与指向目标细胞位移之间的余弦相似度。M2 使用参数 $\alpha$ 和 $\kappa$ 将条件运输概率与速度方向相结合；行归一化后再恢复 M1 的生长质量。

## 10. Optional cross-species prior / 可选的跨物种先验

**English.** M3 can add a reviewed cross-species lineage-topology prior as a soft penalty. It is disabled by default because orthology and cell-type similarity should not override the observed zebrafish expression and velocity evidence.

**中文。** M3 可将经过审核的跨物种谱系拓扑作为软惩罚加入模型。默认关闭该功能，因为直系同源关系和细胞类型相似性不应覆盖斑马鱼自身的表达和速度证据。

## 11. Validation and interpretation / 验证与解释

**English.** Raw mass and conditional fate are aggregated separately. M0 versus M1 tests the effect of growth; M1 versus M2 tests directional reweighting. Sensitivity to $\kappa$ is reported without refitting moscot. Orthology supports a conserved human-equivalent proliferation program, but it does not prove identical regulation of every gene in both species.

**中文。** 原始运输质量与条件命运概率分别汇总。M0 与 M1 的比较用于评估生长信息的影响；M1 与 M2 的比较用于评估方向重加权的影响。在不重新拟合 moscot 的情况下报告 $\kappa$ 敏感性。直系同源证据支持“与人类相对应的保守增殖程序”，但不能证明每个基因在两个物种中具有完全相同的调控方式。

## 12. Correct rerun order / 正确的重新运行顺序

**English.** Reload the notebook from disk, restart the `vae_gpu` kernel, and run from configuration onward. Do not lower the minimum-marker thresholds to hide missing genes. The notebook now selects the correct full-expression object for scoring and fails with explicit diagnostics if cell IDs, markers, or mappings disagree.

**中文。** 从磁盘重新加载 notebook，重启 `vae_gpu` 内核，并从配置单元开始运行。不要通过降低最小标志基因阈值来掩盖基因缺失问题。当前 notebook 会选择正确的完整表达对象进行评分；若细胞 ID、标志基因或映射不一致，将给出明确诊断信息并停止。
