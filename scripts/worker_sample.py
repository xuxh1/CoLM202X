#!/usr/bin/env python3
"""
Worker脚本：在单个节点上处理一个sample的所有站点
使用节点的48个核心并行处理44个站点
"""

import subprocess
import os
import sys
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed

def load_environment(env_file):
    """加载环境变量"""
    cmd = f'bash -c "source {env_file} && env"'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    new_env = {}
    for line in result.stdout.strip().split('\n'):
        if '=' in line:
            key, value = line.split('=', 1)
            if not key.startswith('BASH_FUNC_'):
                new_env[key] = value
    
    os.environ.update(new_env)
    return os.environ.copy()

def run_colm(run_path, nml_path, log_path, nml_name, updated_env):
    """
    处理单个站点的CoLM模拟
    """
    nml_file = f'{nml_path}{nml_name}.nml'
    log_file = f'{log_path}{nml_name}.txt'
    
    try:
        with open(log_file, 'w', encoding='utf-8') as log:
            log.write(f"=== 处理 {nml_name}.nml ===\n")
            log.flush()
            
            commands = [
                [f'{run_path}mksrfdata.x', nml_file],
                [f'{run_path}mkinidata.x', nml_file],
                [f'{run_path}colm.x', nml_file]
            ]
            
            for cmd in commands:
                log.write(f"执行命令: {' '.join(cmd)}\n")
                log.flush()
                
                result = subprocess.run(
                    cmd, 
                    env=updated_env, 
                    stdout=log, 
                    stderr=subprocess.STDOUT, 
                    text=True,
                    timeout=3600  # 1小时超时
                )
                
                if result.returncode != 0:
                    log.write(f"\n命令失败，返回码: {result.returncode}\n")
                    return False
                
                log.write("\n" + "="*50 + "\n")
                log.flush()
            
            log.write(f"=== {nml_name} 处理完成 ===\n")
            log.flush()
        
        return True
    
    except subprocess.TimeoutExpired:
        with open(log_file, 'a', encoding='utf-8') as log:
            log.write(f"\n=== {nml_name} 超时 ===\n")
        return False
    
    except Exception as e:
        with open(log_file, 'a', encoding='utf-8') as log:
            log.write(f"\n=== {nml_name} 异常: {str(e)} ===\n")
        return False

def main(sample_idx, forcing, mode):
    """主函数：处理指定sample的所有站点"""
    sample = f'sample{sample_idx+1}'
    
    print(f"=" * 60)
    print(f"Worker启动: {sample}")
    print(f"节点: {os.environ.get('HOSTNAME', 'unknown')}")
    print(f"可用核心: {os.cpu_count()}")
    print(f"=" * 60)
    
    # 配置路径
    env_file = '/share/home/dq089/soft/gnu-env'
    run_path = '/share/home/dq076/mode/ME/251030_r/run/'
    nml_path = f'{run_path}site/{forcing}/{mode}/{sample}/'
    log_path = f'{nml_path}logs/'
    
    os.makedirs(log_path, exist_ok=True)
    
    # 加载环境
    print("加载环境变量...")
    updated_env = load_environment(env_file)
    
    # 获取所有nml文件
    nml_files = glob.glob(f'{nml_path}*.nml')
    nml_names = [
        os.path.splitext(os.path.basename(nml_file))[0] 
        for nml_file in nml_files 
        if 'HK-MPM' not in nml_file
    ]
    
    print(f"发现 {len(nml_names)} 个站点")
    print(f"站点列表: {nml_names[:5]}{'...' if len(nml_names) > 5 else ''}")
    
    # 并行处理所有站点
    max_workers = min(48, len(nml_names))
    print(f"使用 {max_workers} 个并行进程")
    print("=" * 60)
    
    success_count = 0
    fail_count = 0
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # 提交所有任务
        future_to_nml = {
            executor.submit(run_colm, run_path, nml_path, log_path, nml_name, updated_env): nml_name
            for nml_name in nml_names
        }
        
        # 监控完成情况
        for future in as_completed(future_to_nml):
            nml_name = future_to_nml[future]
            try:
                success = future.result()
                if success:
                    success_count += 1
                    status = "✓"
                else:
                    fail_count += 1
                    status = "✗"
                
                print(f"{status} {nml_name}: 完成 ({success_count + fail_count}/{len(nml_names)})")
                
            except Exception as exc:
                fail_count += 1
                print(f"✗ {nml_name}: 异常 - {exc}")
    
    # 汇总结果
    print("=" * 60)
    print(f"{sample} 处理完成")
    print(f"成功: {success_count}/{len(nml_names)}")
    print(f"失败: {fail_count}/{len(nml_names)}")
    print("=" * 60)
    
    # 返回适当的退出码
    return 0 if fail_count == 0 else 1

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"用法: {sys.argv[0]} <sample_idx> <forcing> <mode>")
        sys.exit(1)
    
    sample_idx = int(sys.argv[1])
    forcing = sys.argv[2]
    mode = sys.argv[3]
    
    exit_code = main(sample_idx, forcing, mode)
    sys.exit(exit_code)
