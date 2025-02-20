# memory_monitor.py
import time
import psutil
import pandas as pd
import sys
print('Check memory')
def monitor_process(pid, interval=1, duration=60):
    # Process 객체 생성
    p = psutil.Process(pid)
    memory_usage = []

    for _ in range(int(duration / interval)):
        try:
            # 메모리 사용량 체크
            memory_info = p.memory_info()
            memory_usage.append({'time': time.time(), 'memory': memory_info.rss})
            time.sleep(interval)
        except psutil.NoSuchProcess:
            break

    # 데이터 프레임으로 변환
    df = pd.DataFrame(memory_usage)
    df.head()
    df.to_csv('./memory_usage.csv', index=False)

if __name__ == "__main__":
    monitor_process(int(sys.argv[1]), interval=1, duration=60) # PID, interval, duration
