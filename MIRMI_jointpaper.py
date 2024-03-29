import pandas as pd
from collections import Counter

# 파일 불러오기
file_path = 'DBLPpub_post_processed_MIRMI.xlsx'
sheet_name = 'database'
data = pd.read_excel(file_path, sheet_name=sheet_name)

# 연도 범위와 특정 학술지들의 목록 설정
year_range = range(2019, 2024)
target_journals = [
    'ICRA', 'IROS', 'Robotics: Science and Systems',
    'CoRL', 'RO-MAN', 'Humanoids', 'ISRR', 'IEEE Robotics Autom. Lett.', 'IEEE Trans. Robotics',
    'Int. J. Robotics Res.', 'Science', 'Nature'
]

university_list = ['TUM']

# 특정 조건에 맞는 행 선택
selected_rows = data[
    (data.iloc[:, 1].isin(year_range)) &
    (data.iloc[:, 2].str.contains('|'.join(target_journals)))
]

# 중복된 행 찾기
duplicated_indices = selected_rows[selected_rows.duplicated(subset=selected_rows.columns[3], keep=False)].index
duplicated_rows = selected_rows.loc[duplicated_indices].sort_values(by=selected_rows.columns[3])

# 중복된 행의 4번째 열 값이 같은 행들끼리 그룹화하여 처리
grouped = duplicated_rows.groupby(duplicated_rows.columns[3])

# 각 그룹에서 5번째 열 값의 빈도수를 셈
result = {}
for name, group in grouped:
    university_counts = Counter(group.iloc[:, 4])
    result[name] = university_counts

# 테이블 초기화
university_journal_table = pd.DataFrame(0, index=university_list, columns=target_journals)

# 그룹의 3번째 열이 열과 일치하고, count가 2 이상인 경우에만 값을 할당
for name, frequencies in result.items():
    for university, count in frequencies.items():
        group_with_university = grouped.get_group(name)
        third_column_values = group_with_university.iloc[:, 2].unique()
        if count >= 2 and university in university_list:
            for third_value in third_column_values:
                # 수정된 부분: 각 target_journal에 대해 부분 일치 여부 확인
                for target_journal in target_journals:
                    if target_journal in third_value:
                        university_journal_table.loc[university, target_journal] += (count-1)

output_file_path = 'duplication.xlsx'
university_journal_table.to_excel(output_file_path)
print(f"Data has been saved to {output_file_path}")
