# worker.py
import os
import subprocess
import json
import glob
from celery import Celery

# 配置 Celery 连接本机的 Redis
broker_url = os.environ.get('CELERY_BROKER_URL', 'redis://localhost:6379/0')
backend_url = os.environ.get('CELERY_RESULT_BACKEND', 'redis://localhost:6379/1')
celery_app = Celery('tnp_worker', broker=broker_url, backend=backend_url)

@celery_app.task(bind=True)
def run_tnp_pipeline(self, task_id, sequence=None, file_path=None):
    """
    后台计算核心节点：继承环境、注入补丁、执行模型
    """
    # 1. 继承系统变量并注入我们修复好的“神级补丁”
    my_env = os.environ.copy()
    my_env["PATH"] = f"/dev/shm/yiqiu_bin:{my_env.get('PATH', '')}"
    my_env["LD_LIBRARY_PATH"] = f"/mnt/ssd0/yiqiu/conda/envs/tnp_api_env/lib:{my_env.get('LD_LIBRARY_PATH', '')}"
    
    output_dir = f"./results/{task_id}"
    os.makedirs(output_dir, exist_ok=True)
    
    # 2. 动态拼装命令（适配单序列和文件上传）
    cmd = ["TNP", "--name", task_id, "--output", output_dir]
    if sequence:
        cmd.extend(["--seq", sequence])
    elif file_path:
        cmd.extend(["--file", file_path])
        
    # 3. 运行子进程
    process = subprocess.run(cmd, env=my_env, capture_output=True, text=True)
    
    if process.returncode != 0:
        # 把底层的报错直接抛出，FastAPI 那边会捕获为 FAILURE 状态
        raise Exception(f"TNP Error: {process.stderr}\n{process.stdout}")
        
    # 4. 智能寻找生成的 JSON 结果并返回
    # TNP 输出的 json 路径可能带有子目录，用 glob 递归寻找最稳妥
    json_files = glob.glob(f"{output_dir}/**/*.json*", recursive=True)
    
    if json_files:
        with open(json_files[0], 'r') as f:
            data = json.load(f)
        return data  # 直接把字典 return 出去，FastAPI 会自动封进 JSON 里
    else:
        raise Exception(f"JSON output not found in {output_dir}")