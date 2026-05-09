# EDTA 安装与转座子分析流程

## EDTA 安装

- GitHub 地址：https://github.com/oushujun/EDTA/releases

### 1. 创建 Conda 环境
```bash
conda env create -f EDTA2.2.yml
```

### 2. RepeatMasker 添加 RepBase 数据库
```bash
# 下载 RepBase
wget http://ftp.genek.cn:8888/Share/linux_software/RepBaseRepeatMaskerEdition-20181026.tar.gz

# 解压并拷贝到 RepeatMasker 目录下，再解压缩
# 然后执行：
share/RepeatMasker/addRepBase.pl -libdir Libraries
```

---

## 一、构建 genome.fai
```bash
samtools faidx Ath.genome.fa
```

---

## 二、简单粗暴：直接运行 RepeatMasker
```bash
RepeatMasker --species "Arabidopsis lyrata" -e rmblast -pa 10 -qq Aly.genome.fa
```

---

## 三、RepeatModeler + RepeatMasker 流程

### 1. 构建数据库
```bash
BuildDatabase -name genome Aly.genome.fa
```

### 2. 运行 RepeatModeler
```bash
RepeatModeler -database genome -pa 8 -LTRStruct
```
> **注意**：EDTA 下的 RepeatModeler 不能用 `-LTRStruct` 参数。

### 3. 运行 RepeatMasker
```bash
RepeatMasker -e rmblast -pa 8 -qq -lib genome-families.fa Aly.genome.fa
```

---

## 四、EDTA 流程

### 1. 运行 EDTA

| 参数 | 说明 |
|:---|:---|
| `--species` | 指定物种，用于 TIR 鉴定。可选：`Rice` \| `Maize` \| `others` |
| `--step` | 运行 EDTA 的哪一步：`all` \| `filter` \| `final` \| `anno` |
| `--sensitive [0\|1]` | 是否运行 RepeatModeler，`0` 否，`1` 是 |
| `--anno [0\|1]` | 是否进行注释（即运行 RepeatMasker） |
| `--overwrite [0\|1]` | 是否覆盖之前结果 |
| `--curatedlib <File>` | 提供自定义数据库 |
| `--cds <File>` | 进行 CDS 过滤 |
| `-u` | LTR-Retriever 的自然突变率，默认 `1.3e-8`（水稻） |

```bash
EDTA.pl --genome Aly.genome.fa \
        --species others \
        --sensitive 1 \
        --anno 1 \
        --threads 20
```

### 结果文件说明

#### （1）`genome.fa.mod.EDTA.raw`
- `LTR`、`Helitron`、`TIR` 目录：含各类型原始注释结果
- `.EDTA.intact.fa`：所有种类 repeat 的序列
- `.intact.gff3`：注释文件

#### （2）`genome.fa.mod.EDTA.raw/LTR`
- `.LTRlib.fa`：LTR 的 lib
- `.intact.fa`：所有完整 LTR 序列
- `.pass.list.gff3`：完整 LTR 的内部注释
- `.pass.list`：LTR-Retriever 运行结果

#### （3）`genome.fa.mod.EDTA.final`（final 步结果）
- `TElib.fa`：最终的 lib

#### （4）`genome.fa.mod.EDTA.anno`（anno 步结果）
- `TEanno.gff3`：最终注释结果（含 Identity 信息等）
- `TEanno.sum`：统计结果

### 补充：EDTA_raw.pl

可单独运行 LTR-Retriever、TIR-Learner、HelitronScanner（如计算 LTR 插入时间）。

| 参数 | 说明 |
|:---|:---|
| `--type` | 可选：`ltr` \| `tir` \| `helitron` \| `all` |

```bash
EDTA_raw.pl --genome genome.fa \
            --species others \
            -u 1.3e-8 \
            --type ltr \
            --threads 20
```

---

## 常见问题一：SINE 和 LINE 注释缺失

以下三种解决方法：

### 方法（1）：使用 Repbase 构建 curatedlib
```bash
famdb.py families -d Viridiplantae \
         --format fasta_name \
         --include-class-in-name > Viridiplantae.repeatmasker.fasta

# EDTA.pl 加 --curatedlib RepBase_rosids.fa
```

### 方法（2）：使用 `--sensitive 1`，运行 RepeatModeler

### 方法（3）：使用 SINE Base 数据库
- 数据库地址：https://sines.eimb.ru
- 点击 Downloads，下载 SINEBank 和 LINEBank
- ID 格式需处理为标准格式，如 `>AFC-2#SINE/AFC-2`；也可只分到 SINE 层次：`>AFC-2#SINE`

```bash
# 处理 SINE
awk '{print $1}' SINEs.bnk | \
    awk '{if(/>/) print $1"#SINE/"$1; else print $1}' | \
    sed 's/\/>/\//g' > SINE.fa

# 处理 LINE
awk '{print $1}' LINEs.bnk | \
    awk '{if(/>/) print $1"#LINE/"$1; else print $1}' | \
    sed 's/\/>/\//g' > LINE.fa

# 合并
cat SINE.fa LINE.fa > SINE_LINE.fa

# 运行 EDTA.pl 时加 --curatedlib SINE_LINE.fa，可以 --sensitive 0
```

---

## 五、对结果中的 unknown 进一步分类 — DeepTE

- GitHub 地址：https://github.com/LiLabAtVT/DeepTE

### 环境创建
```bash
conda install tensorflow-gpu=1.14.0 biopython keras=2.2.4 numpy=1.16.0 hmmer
```

### 1. 提取 LTR/unknown 序列
```bash
grep "LTR/unknown" genome.fa.mod.EDTA.TElib.fa | \
    sed 's/>//' | \
    seqtk subseq genome.fa.mod.EDTA.TElib.fa - > LTR_unknown.fa
```

### 2. 运行 DeepTE
> 需在 GitHub 下载 DeepTE 的 model（也可用基因课下载完的），然后 `tar -xvf` 解压。

| 参数 | 说明 |
|:---|:---|
| `-i` | 指定待鉴定序列 |
| `-sp` | 指定物种：`P`（植物）\| `M`（动物）\| `F`（真菌）\| `O`（其他） |
| `-m_dir` | 指定下载的 model |
| `-fam` | 指定 repeat 种类，也可不加（速度慢点） |

```bash
DeepTE.py -i LTR_unknown.fa \
          -sp P \
          -m_dir Plants_model \
          -fam LTR
```

#### 结果
- `opt_DeepTE.txt`：分类结果
- `opt_DeepTE.fasta`：将新分的类名加在 ID 后面

### 3. 合并序列

#### （1）提取已知分类的序列
```bash
grep -v "LTR/unknown" genome.fa.mod.EDTA.TElib.fa | \
    sed 's/>//' | \
    seqtk subseq genome.fa.mod.EDTA.TElib.fa - > LTR_known.fa
```

#### （2）处理 DeepTE 结果 ID，使其符合标准
```bash
sed 's/LTR\/unknown__ClassI_LTR_Copia/LTR\/Copia/' opt_DeepTE.fasta | \
    sed 's/LTR\/unknown__ClassI_LTR_Gypsy/LTR\/Gypsy/' | \
    sed 's/LTR\/unknown__ClassI_LTR/LTR\/unknown/' > LTR_unknown_DeepTE.fa
```

#### （3）合并
```bash
cat LTR_known.fa LTR_unknown_DeepTE.fa > genome.fa.mod.EDTA.TElib.fa
```

### 4. 重新运行 EDTA
```bash
# 将 TElib.fa 放在 final 目录下
cp genome.fa.mod.EDTA.TElib.fa genome.fa.mod.EDTA.final/

# 重新注释
EDTA.pl --genome genome.fa \
        --step anno \
        --overwrite 1 \
        --anno 1
```

---

## 补充：计算 LAI 指数

### 1. 输入文件
- `genome.fa.mod.EDTA.raw/LTR/genome.fa.mod.pass.list`
- `genome.fa.mod.EDTA.final/genome.fa.mod.out`

### 2. 运行 LAI
```bash
nohup LAI -t 24 \
          -genome genome.fa \
          -intact genome.fa.mod.pass.list \
          -all genome.fa.mod.out \
          &> lai.out &
```

| 参数 | 说明 |
|:---|:---|
| `-t 24` | BLAST 使用的线程数量 |
| `-genome` | 指定基因组 |
| `-intact` | LTR-Retriever 生成的非冗余 LTR-RT 文库列表 |
| `-all` | RepeatMasker 注释的所有 LTR 序列 |
| `-window 3000000` | 计算 LAI 的窗口大小，默认 3 Mb |
| `-step 300000` | 计算 LAI 的步长，默认 300 Kb |
| `-q` | 快速评估模式，建议大基因组使用，牺牲 0.5% 精度 |
| `-qq` | 超快速模式，不评估 LTR identity，只输出 raw_LAI（用于种间比较） |
| `-mono [file]` | 用文件提供序列名称，LAI 只计算指定序列。主要用于多倍体基因组，提供代表单倍体的序列 |

### 3. 结果
- `sample.fa.mod.out.LAI`：保存 raw_LAI 值和 LAI 值（全基因组和滑窗范围）
- `sample.fa.mod.out.LAI.LTR.ava.age`
- `sample.fa.mod.out.LAI.LTR.ava.out`