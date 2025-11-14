#!/usr/bin/env python3
"""
主脚本：提交4个runoff任务到LSF队列
每个任务在单独节点上运行，处理一个runoff的所有站点
"""

import subprocess
import os
import time
from pathlib import Path

def submit_runoff_job(runoff_idx, forcing='FLUXNET-CH4', mode='no_spin_up',sample='standard'):
    """提交单个runoff的LSF任务"""
    runoff = f'runoff{runoff_idx}'
    
    # LSF任务脚本路径
    run_path = '/share/home/dq076/mode/ME/251030_r/run/'
    worker_script = f'/share/home/dq076/mode/ME/251030_r/scripts/worker_runoff.py'
    log_dir = f'{run_path}site/{forcing}/{mode}/{sample}/{runoff}/bsub_logs/'
    os.makedirs(log_dir, exist_ok=True)
    
    # LSF任务参数
    job_name = f'{runoff}'
    stdout_log = f'{log_dir}stdout.log'
    stderr_log = f'{log_dir}stderr.log'
    
    # 构建bsub命令
    bsub_cmd = [
        'bsub',
        '-J', job_name,                    # 任务名
        '-n', '48',                        # 48核
        '-q', 'normal',                    # 队列名（根据你的系统调整）
        '-o', stdout_log,                  # 标准输出
        '-e', stderr_log,                  # 标准错误
        '-R', 'span[ptile=48]',           # 所有核心在同一节点
        '-R', 'rusage[mem=150G]',           # 所有核心在同一节点
        'python3', worker_script,          # 执行的脚本
        str(runoff_idx),                   # 传递runoff索引
        forcing,
        mode,
        sample
    ]
    
    # 提交任务
    result = subprocess.run(bsub_cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        # 提取job ID
        job_id = result.stdout.strip().split('<')[1].split('>')[0] if '<' in result.stdout else 'unknown'
        print(f"✓ 提交 {runoff}: Job ID = {job_id}")
        return job_id
    else:
        print(f"✗ 提交 {runoff} 失败: {result.stderr}")
        return None

def check_job_status(job_ids):
    """检查任务状态"""
    result = subprocess.run(['bjobs', '-w'], capture_output=True, text=True)
    
    if result.returncode != 0:
        return {}
    
    # 解析bjobs输出
    lines = result.stdout.strip().split('\n')
    if len(lines) < 2:
        return {}
    
    status_dict = {}
    for line in lines[1:]:  # 跳过表头
        parts = line.split()
        if len(parts) >= 3:
            job_id = parts[0]
            status = parts[2]
            if job_id in job_ids:
                status_dict[job_id] = status
    
    return status_dict

def main():
    forcing = 'FLUXNET-CH4'
    mode = 'no_spin_up'
    sample = 'standard'
    num_runoffs = 4

    print(f"=" * 60)
    print(f"开始提交 {num_runoffs} 个runoff任务到LSF队列")
    print(f"每个任务使用1个节点48核")
    print(f"=" * 60)
    
    # 提交所有任务
    job_ids = []
    runoff_to_job = {}
    
    for i in range(num_runoffs):
        job_id = submit_runoff_job(i, forcing, mode, sample)
        if job_id:
            job_ids.append(job_id)
            runoff_to_job[job_id] = f'runoff{i}'
        time.sleep(0.1)  # 避免过快提交
    
    print(f"\n成功提交 {len(job_ids)} 个任务")
    print(f"失败: {num_runoffs - len(job_ids)} 个任务\n")
    
    if not job_ids:
        print("没有任务成功提交，退出。")
        return
    
    # 监控任务进度
    print("开始监控任务状态...")
    print("=" * 60)
    
    completed = set()
    last_status = {}
    
    while len(completed) < len(job_ids):
        status_dict = check_job_status(job_ids)
        
        # 检查状态变化
        for job_id in job_ids:
            if job_id in completed:
                continue
            
            current_status = status_dict.get(job_id, 'DONE')
            
            # 如果任务不在bjobs输出中，说明已完成
            if job_id not in status_dict:
                completed.add(job_id)
                runoff = runoff_to_job[job_id]
                print(f"✓ {runoff} (Job {job_id}) 已完成")
                continue
            
            # 状态变化时打印
            if current_status != last_status.get(job_id):
                runoff = runoff_to_job[job_id]
                print(f"  {runoff} (Job {job_id}): {last_status.get(job_id, 'PEND')} -> {current_status}")
                last_status[job_id] = current_status
        
        # 显示当前进度
        pend = sum(1 for s in status_dict.values() if s == 'PEND')
        run = sum(1 for s in status_dict.values() if s == 'RUN')
        done = len(completed)
        
        print(f"\r进度: 等待={pend}, 运行={run}, 完成={done}/{len(job_ids)}\n", end='', flush=True)
        
        time.sleep(60)  # 每60秒检查一次
    
    print(f"\n\n{'=' * 60}")
    print(f"所有任务完成！")
    print(f"{'=' * 60}")

if __name__ == "__main__":
    main()
