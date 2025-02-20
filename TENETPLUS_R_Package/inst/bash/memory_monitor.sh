#!/bin/bash


# 모니터링할 시간 간격(초)
INTERVAL=1
# 출력 파일
OUTPUT_FILE="system_memory_usage.csv"

# 종료 시 실행할 함수 정의
function finish {
    echo "Memory monitoring stopped. Data saved to $OUTPUT_FILE"
    exit
}

# 스크립트 종료 시 finish 함수 실행
trap finish EXIT INT

# 초기 메모리 사용량 측정
read initial_total initial_used initial_free <<< $(free | grep Mem | awk '{print $2, $3, $4}')

# CSV 파일 헤더
echo "time,total,used,free" > $OUTPUT_FILE

# 메모리 모니터링 루프
while true; do
    # 현재 메모리 사용량 측정
    read total used free <<< $(free | grep Mem | awk '{print $2, $3, $4}')

    # 측정값을 초기값과 비교하여 조정
    adjusted_total=$((total - initial_total))
    adjusted_used=$((used - initial_used))
    adjusted_free=$((free - initial_free))

    # 현재 시간과 조정된 메모리 사용량을 CSV 파일에 기록
    echo "$(date +%Y-%m-%dT%H:%M:%S), $adjusted_total, $adjusted_used, $adjusted_free" >> $OUTPUT_FILE

    # INTERVAL만큼 대기
    sleep $INTERVAL
done