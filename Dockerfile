# 使用 Miniconda 基础镜像
FROM continuumio/miniconda3:latest

# 设置工作目录
WORKDIR /app

# 1. 替换 apt 为国内清华源并安装基础构建工具
RUN sed -i 's/deb.debian.org/mirrors.tuna.tsinghua.edu.cn/g' /etc/apt/sources.list.d/debian.sources && \
    apt-get update && \
    apt-get install -y build-essential && \
    rm -rf /var/lib/apt/lists/*

# 2. 配置 Conda 清华源
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ && \
    conda config --set show_channel_urls yes


RUN conda install -y python=3.10
 
# 3. 把当前项目代码拷贝进容器
COPY . /app

# 4. 赋予底层二进制工具最高权限 (在 Docker 里我们是 Root，随便赋权)
RUN chmod 777 /app/bin/*

# 5. 配置 pip 阿里源并安装 Python 依赖
RUN pip config set global.index-url https://mirrors.aliyun.com/pypi/simple/ && \
    pip install fastapi uvicorn python-multipart celery redis pandas && \
    pip install biopython==1.77 && \
    pip install -e .

# 6. 安装底层生物学依赖 (一锅端)
RUN conda install -y hmmer anarci pdbfixer libstdcxx-ng "dssp<4.0.0"

# 7. 修复 DSSP 命名问题
RUN ln -sf /opt/conda/bin/mkdssp /opt/conda/bin/dssp

# 暴露端口
EXPOSE 8001