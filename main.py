# main.py
from fastapi import FastAPI, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from celery.result import AsyncResult
import shutil
import os
import uuid

# 导入 worker 中的 Celery 任务
from worker import run_tnp_pipeline 

app = FastAPI(title="TNP Descriptor API", description="供前端调用的抗体特征计算服务")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], 
    allow_methods=["*"],
    allow_headers=["*"],
)

UPLOAD_DIR = "./uploads"
os.makedirs(UPLOAD_DIR, exist_ok=True)

@app.post("/api/v1/submit_job", summary="提交计算任务")
async def submit_job(
    file: UploadFile = File(None, description="上传 FASTA 格式文件（可选）"),
    sequence: str = Form(None, description="直接输入单条氨基酸序列（可选）")
):
    """前端可以选择上传文件，或者直接传文本序列。二选一。"""
    task_id = str(uuid.uuid4())
    
    # 模式 1：前端直接传了纯文本序列 (对应 TNP 的 --seq)
    if sequence:
        task = run_tnp_pipeline.delay(task_id, sequence=sequence, file_path=None)
        return {"message": "Sequence task submitted", "task_id": task.id}
        
    # 模式 2：前端上传了文件 (对应 TNP 的 --file)
    elif file:
        file_path = os.path.join(UPLOAD_DIR, f"{task_id}_{file.filename}")
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        task = run_tnp_pipeline.delay(task_id, sequence=None, file_path=file_path)
        return {"message": "File task submitted", "task_id": task.id}
        
    else:
        return {"error": "Must provide either 'file' or 'sequence'"}

@app.get("/api/v1/job_status/{task_id}", summary="轮询任务状态与结果")
async def get_job_status(task_id: str):
    task_result = AsyncResult(task_id)
    
    if task_result.state in ['PENDING', 'STARTED']:
        return {"task_id": task_id, "status": "processing"}
    elif task_result.state == 'SUCCESS':
        return {"task_id": task_id, "status": "completed", "data": task_result.result}
    elif task_result.state == 'FAILURE':
        return {"task_id": task_id, "status": "failed", "error": str(task_result.info)}
    else:
        return {"task_id": task_id, "status": task_result.state}