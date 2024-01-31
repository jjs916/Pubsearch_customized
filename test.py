import pandas as pd

# 가상의 데이터프레임 생성
data = pd.DataFrame({
    'Column1': ['apple', 'orange', 'banana', 'grape', 'kiwi', 'apple tree', 'apple juice'],
    'Column2': ['red', 'orange', 'yellow', 'purple', 'green', 'brown', 'yellow']
})

# 'apple'이 포함되지 않는 행 선택
filtered_rows = data[~data['Column1'].str.contains('apple', case=False, na=False)]

# 결과 출력
print(filtered_rows)
